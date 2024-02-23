#include "Graph.h"
#include "Utility.h"
#include "KPlex_matrix.h"

using namespace std;

Graph::Graph(const int _K)
{
    kplex.clear();

    K = _K;

    n = m = pm = nm = 0;

    pstart = nullptr;
    pend = nullptr;
    edges = nullptr;

    p_pstart = nullptr;
    p_pend = nullptr;
    p_edges = nullptr;

    n_pstart = nullptr;
    n_pend = nullptr;
    n_edges = nullptr;

	degree = nullptr;
    p_degree = nullptr;
    n_degree = nullptr;
    tri_cnt = nullptr;

    s_n = 0;

    lb = 0, ub = 1e9;
}

Graph::~Graph()
{
	if(pstart != nullptr) {
		delete[] pstart;
		pstart = nullptr;
	}
	if(pend != nullptr) {
		delete[] pend;
		pend = nullptr;
	}
	if(edges != nullptr) {
		delete[] edges;
		edges = nullptr;
	}
	if(p_pstart != nullptr) {
		delete[] p_pstart;
		p_pstart = nullptr;
	}
	if(p_pend != nullptr) {
		delete[] p_pend;
		p_pend = nullptr;
	}
	if(p_edges != nullptr) {
		delete[] p_edges;
		p_edges = nullptr;
	}
	if(n_pstart != nullptr) {
		delete[] n_pstart;
		n_pstart = nullptr;
	}
	if(n_pend != nullptr) {
		delete[] n_pend;
		n_pend = nullptr;
	}
	if(n_edges != nullptr) {
		delete[] n_edges;
		n_edges = nullptr;
	}
	if(degree != nullptr) {
		delete[] degree;
		degree = nullptr;
	}
	if(p_degree != nullptr) {
		delete[] p_degree;
		p_degree = nullptr;
	}
	if(n_degree != nullptr) {
		delete[] n_degree;
		n_degree = nullptr;
	}
	if(tri_cnt != nullptr) {
		delete[] tri_cnt;
		tri_cnt = nullptr;
	}
}

void Graph::load_graph(string input_graph)
{
	Timer t;
	t.restart();
    string buffer;
    ifstream input_file(input_graph, ios::in);

    if (!input_file.is_open()) {
        cout<<"cannot open file : "<<input_graph<<endl;exit(1);
    }
    else {
		input_file >> n >> m;
		map<ui, int> * s_G = new map<ui, int>[n];
		ui u, v;
		int flag;
		while (input_file >> u >> v >> flag) {
			if (u == v) continue;
            assert(u >= 0 && u < n);
            assert(v >= 0 && v < n);
            assert(flag == 1 || flag == -1);
            s_G[u].insert(make_pair(v, flag));
            s_G[v].insert(make_pair(u, flag));
		}

		m = pm = nm = 0;
		for (ui i = 0; i < n; i++) {
			const map<ui, int> & nei = s_G[i];
			for (auto e : nei) {
				if (e.second == 1) ++pm;
				else ++nm;
			}
			m += nei.size();
		}
		assert(m%2 == 0);assert(pm%2 == 0);assert(nm%2 == 0);
		m /= 2; pm /= 2; nm /= 2;

        input_file.close();

		pstart = new ept[n+1];
		pend = new ept[n];
		edges = new ui[2*m];

		p_pstart = new ept[n+1];
        p_pend = new ept[n];
		p_edges = new ui[2*pm];

        n_pstart = new ept[n+1];
        n_pend = new ept[n];
        n_edges = new ui[2*nm];

        degree = new ui[n];
        p_degree = new ui[n];
        n_degree = new ui[n];

        tri_cnt = new ept[2*m];

        //construct positive edges
        p_pstart[0] = 0;
		for (ui i = 0; i < n; i++) {
			const map<ui, int> & nei = s_G[i];
			ui idx = p_pstart[i];
			for (auto e : nei) if (e.second == 1) {
				p_edges[idx++] = e.first;
			}
			p_pend[i] = p_pstart[i+1] = idx;
		}
		assert(p_pstart[n] == 2*pm);

        //construct negative edges
        n_pstart[0] = 0;
		for (ui i = 0; i < n; i++) {
			const map<ui, int> & nei = s_G[i];
			ui idx = n_pstart[i];
			for (auto e : nei) if (e.second == -1) {
				n_edges[idx++] = e.first;
			}
			n_pend[i] = n_pstart[i+1] = idx;
		}
		assert(n_pstart[n] == 2*nm);

        //construct edges
        pstart[0] = 0;
		for (ui i = 0; i < n; i++) {
			const map<ui, int> & nei = s_G[i];
			ui idx = pstart[i];
			for (auto e : nei) {
				edges[idx++] = e.first;
			}
			pend[i] = pstart[i+1] = idx;
		}
		assert(pstart[n] == 2*m);

        delete [] s_G;

		for(ui i = 0; i < n; i++) {
			sort(edges+pstart[i], edges+pend[i]);
			sort(p_edges+p_pstart[i], p_edges+p_pend[i]);
			sort(n_edges+n_pstart[i], n_edges+n_pend[i]);
		}
    }

	cout<<"\t load_graph: time cost = "<<integer_to_string(t.elapsed())<<endl;
	cout<<"\t G : n = "<<n<<", m = "<<m<<", pm = "<<pm<<", nm = "<<nm<<endl;
}

void Graph::get_k_core(int k)
{
	Timer t;
	t.restart();

	bool *v_del = new bool[n];
	memset(v_del, 0, sizeof(bool)*n);

    get_degree();
    queue<ui> q;
    for(ui i = 0; i < n; i++) if(degree[i] < k) q.push(i);
	while(!q.empty()) {
		ui u = q.front(); q.pop();
		v_del[u] = 1;
        for(ept i = pstart[u]; i < pend[u]; i++) {
			ui v = edges[i];
			if ((degree[v]--) == k) {
				q.push(v);
			}
        }
	}

    //rebuild
    rebuild_graph(v_del);

	delete [] v_del;

    cout<<"\t get_k_core, time cost = "<<integer_to_string(t.elapsed())<<",\t k = "<<k<<", n = "<<n<<", m = "<<m<<endl;
}

//get degree for all nodes of G
void Graph::get_degree()
{
	for (ui i = 0; i < n; i++) {
        p_degree[i] = p_pstart[i+1] - p_pstart[i];
        n_degree[i] = n_pstart[i+1] - n_pstart[i];
        degree[i] = p_degree[i] + n_degree[i];
	}
}

//get triangle count for all edges of G
void Graph::get_tricnt()
{
	ui *mark = new ui[n];
	memset(mark, 0, sizeof(ui)*n);
    memset(tri_cnt, 0, sizeof(ept)*m*2);

	for(ui u = 0; u < n; u++) {
		for(ept j = pstart[u]; j < pend[u]; j++) mark[edges[j]] = j+1;

		for(ept j = pstart[u]; j < pend[u]; j++) {
			ui v = edges[j];
			for(ept k = pstart[v]; k < pend[v]; k++) {
				ui w = edges[k];
				if(mark[w] && v < w) {
					tri_cnt[j]++;
					tri_cnt[k]++;
					tri_cnt[mark[w]-1]++;
				}
			}
		}

		for(ept j = pstart[u]; j < pend[u]; j++) mark[edges[j]] = 0;
	}

	delete [] mark;
}

void Graph::get_g(ui u, vector<pair<int,int> > &vp, vector<int> &sgn)
{
    Timer t;
	t.restart();

	bool *v_sel = new bool[n];
	memset(v_sel, 0, sizeof(bool)*n);
    v_sel[u] = 1;

	for(ept i = pstart[u]; i < pend[u]; i++) {
		ui v = edges[i];
		v_sel[v] = 1;
		for(ept j = pstart[v]; j < pend[v]; j++) {
			ui w = edges[j];
			v_sel[w] = 1;
		}
	}

    //build adjacent matrix
    ui *rid = new ui[n];
	ui cnt = 0;
	for(ui i = 0; i < n; i++) if(v_sel[i]) {
		rid[i] = cnt++;
	}
    s_n = cnt;
    
    for(ui u = 0; u < n; u++) if(v_sel[u]) {
        for(ept i = p_pstart[u]; i < p_pend[u]; i++) {
            ui v = p_edges[i];
            if(v_sel[v]) {
                if(u < v) {
                    vp.push_back(make_pair(rid[u], rid[v]));
					sgn.push_back(1);
                }
            }
        }
        for(ept i = n_pstart[u]; i < n_pend[u]; i++) {
            ui v = n_edges[i];
            if(v_sel[v]) {
                if(u < v) {
                    vp.push_back(make_pair(rid[u], rid[v]));
					sgn.push_back(-1);
                }
            }
        }
    }

	delete [] rid;
	delete [] v_sel;

#ifndef NDEBUG
	cout<<"\t get_g, T : "<<integer_to_string(t.elapsed())<<",\t n="<<s_n<<endl;
#endif
}

void Graph::rebuild_graph(bool *v_del)
{
	ui *rid = new ui[n];
	ui cnt = 0;
	for(ui i = 0; i < n; i++) if(!v_del[i]) {
		rid[i] = cnt++;
	}

    if(cnt != n) {
        cnt = 0;
		ept pos = 0, p_pos = 0, n_pos = 0;
        pstart[0] = p_pstart[0] = n_pstart[0] = 0;
		for(ui u = 0; u < n; u++) if(!v_del[u]) {
			for(ept i = pstart[u]; i < pend[u]; i++) {
				ui v = edges[i];
				if(!v_del[v]) {
					edges[pos++] = rid[v];
				}
			}
			pend[cnt] = pos;

			for(ept i = p_pstart[u]; i < p_pend[u]; i++) {
				ui v = p_edges[i];
				if(!v_del[v]) {
					p_edges[p_pos++] = rid[v];
				}
			}
			p_pend[cnt] = p_pos;

			for(ept i = n_pstart[u]; i < n_pend[u]; i++) {
				ui v = n_edges[i];
				if(!v_del[v]) {
					n_edges[n_pos++] = rid[v];
				}
			}
			n_pend[cnt] = n_pos;
			cnt++;
		}

		n = cnt;
		m = pos / 2;
		pm = p_pos / 2;
		nm = n_pos / 2;

		for(ui u = 1; u <= n; u++) {
			pstart[u] = pend[u-1];
			p_pstart[u] = p_pend[u-1];
			n_pstart[u] = n_pend[u-1];
		}
    }

    delete[] rid;
}

void Graph::CTCP(int del_v, bool lb_changed, int tv, int te)
{
    printf("\t CTCP: tv = %d, te = %d\n", tv, te);
	tv = max(0, tv); te = max(0, te);
    Timer t;
	t.restart();

	queue<ui> qv;
	queue<pair<ui, ept>> qe; //from, idx
	get_degree();
	get_tricnt();
	bool *v_del = new bool[n];
	bool *e_del = new bool[m*2];
    ui *mark = new ui[n];
	memset(v_del, 0, sizeof(bool)*n);
	memset(e_del, 0, sizeof(bool)*m*2);
    memset(mark, 0, sizeof(ui)*n);

    if(del_v != -1) qv.push((ui)del_v);
    if(lb_changed) {
        for(ui u = 0; u < n; u++) {
            if(degree[u] < tv) qv.push(u);
            for(ept i = pstart[u]; i < pend[u]; i++) {
                ui v = edges[i];
                if(u < v && tri_cnt[i] < te) {
                    // printf("%d %d\n", u, v);
                    qe.push(make_pair(u, i));
                }
            }
        }
    }

    while(!qv.empty() || !qe.empty()) {
        while(!qe.empty()) {
			auto ue = qe.front(); qe.pop();
			ui u = ue.first; ept id_uv = ue.second; ui v = edges[id_uv];
			ui id_vu = pstart[v] + find(edges + pstart[v], edges + pend[v], u);
			assert(e_del[id_uv] == e_del[id_vu]);
			if(v_del[u] || v_del[v] || e_del[id_uv]) continue;
            e_del[id_uv] = 1;
			e_del[id_vu] = 1;

            if((degree[u]--) == tv) qv.push(u);
            if((degree[v]--) == tv) qv.push(v);

            for(ept i = pstart[u]; i < pend[u]; i++) if(!e_del[i]) mark[edges[i]] = i+1;

            for(ept j = pstart[v]; j < pend[v]; j++) if(!e_del[j]) {
                ui w = edges[j];
                if(!v_del[w] && mark[w]) { // triangle count--
                    ui id_uw = mark[w] - 1;
                    ui id_wu = pstart[w] + find(edges + pstart[w], edges + pend[w], u);
                    if((tri_cnt[id_uw]--) == te) qe.push(make_pair(u, id_uw));
                    tri_cnt[id_wu]--;

                    ui id_vw = j;
                    ui id_wv = pstart[w] + find(edges + pstart[w], edges + pend[w], v);
                    if((tri_cnt[id_vw]--) == te) qe.push(make_pair(v, id_vw));
                    tri_cnt[id_wv]--;
                }
            }

            for(ept i = pstart[u]; i < pend[u]; i++) if(!e_del[i]) mark[edges[i]] = 0;
        }
        if(!qv.empty()) {
            ui u = qv.front(); qv.pop();
			// printf("%d\n", u);
            if(v_del[u]) continue;
            v_del[u] = 1;

			for(ept i = pstart[u]; i < pend[u]; i++) {
				ui v = edges[i];
				if(!e_del[i] && !v_del[v]) {
					if((degree[v]--) == tv) qv.push(v);
				}
				e_del[i] = 1;
			}

            for(ept i = pstart[u]; i < pend[u]; i++) if(!e_del[i]) mark[edges[i]] = i+1;

			for(ept i = pstart[u]; i < pend[u]; i++) {
				ui v = edges[i]; 
                if(e_del[i] || v_del[v]) continue;
                for(ept j = pstart[v]; j < pend[v]; j++) {
                    ui w = edges[j];
                    if(e_del[j] || v_del[w]) continue;
                    if(mark[w] && w > v) { // triangle count--
                        ui id_vw = j;
                        ui id_wv = pstart[w] + find(edges + pstart[w], edges + pend[w], v);
                        if((tri_cnt[id_vw]--) == te) qe.push(make_pair(u, id_vw));
                        tri_cnt[id_wv]--;
                    }
                }
			}

            for(ept i = pstart[u]; i < pend[u]; i++) if(!e_del[i]) mark[edges[i]] = 0;
        }
    }

    rebuild_graph(v_del);

    delete[] v_del;
    delete[] e_del;

// #ifndef NDEBUG
    cout<<"\t CTCP, T : "<<integer_to_string(t.elapsed())<<",\t n="<<n<<", m="<<m<<endl;
// #endif
}

// // degeneracy-based k-plex
// // return an upper bound of the maximum k-plex size
// // return dOrder
// void kplex_degen(ListLinearHeap *heap, int k, int *dOrder)
// {
// 	Timer t;
// 	int *peel_sequence = new int[G.n];
// 	int *vis = new int[G.n];
	
// 	for(int i = 0; i < G.n; i++) peel_sequence[i] = i;
// 	for(int i = 0; i < G.n; i++) vis[i] = 0;

// 	heap->init(G.n, G.n-1, peel_sequence, G.degree);
// 	for(int i = 0; i < G.n; i++) {
// 		int u, deg;
// 		heap->pop_min(u, deg);
// 		dOrder[i] = u;

// 		// if(deg+k >= G.n-i+1 && G.n-i+1 > lb) lb = G.n-i+1;
// 		ub = max(ub, min(deg+k, G.n-i+1));

// 		for(int j = G.pstart[u]; j < G.pend[u]; j++) {
// 			int v = G.edges[j];
// 			if(vis[v] == 0) heap->decrement(v, 1);
// 		}
// 		vis[u] = 1;
// 	}

// 	// for(int i = 0; i < G.n; i++) {
// 	// 	printf("%d ", dOrder[i]);
// 	// }
// 	// printf("\n");

// 	delete [] peel_sequence;
// 	delete [] vis;

// }

void Graph::heu_signed_kplex(int rounds, int k)
{
	Timer t; t.restart();
	if(rounds < 1) return;

	KPLEX_MATRIX *kplex_solver = new KPLEX_MATRIX();
	kplex_solver->allocateMemory(n, m/2);

	vector<pair<int,int> > vp; vp.reserve(m/2);
	vector<int> sgn; sgn.reserve(m/2);

	for(int round = 0; round < rounds && round < n; round++) {
		get_degree();
		ui u = 0;
		for(ui i = 1; i < n; i++) if(degree[u] > degree[i]) u = i;
		printf("%d %d\n", u, degree[u]);

		vp.clear();
		sgn.clear();
		get_g(u, vp, sgn);

        //heuristically find signed kplex
		kplex_solver->load_graph(s_n, vp, sgn);
		// kplex_solver->heu_kPlex(K, kplex);
// #ifndef NDEBUG
		// printf("\t heu_kplex = %d \n", (int)kplex.size());
		// vector<ui> KPLEX; KPLEX.clear();
		kplex_solver->kPlex(K, kplex, true);
		printf("\t KPLEX = %d \n", (int)kplex.size());
// #endif
		bool lb_changed = 0;
		if(kplex.size() > lb) {
			if(kplex.size() > lb) lb_changed = 1;
			lb = kplex.size();
		}
        // CTCP
		CTCP(u, lb_changed, lb+1-K, lb+1-2*K);
        // CTCP(-1, lb_changed, lb+1-K, lb+1-2*K);
	}
	cout<<"\t heu_find_kplex, T : "<<integer_to_string(t.elapsed())<<",\t heu_kplex_size = "<<kplex.size()<<endl;
}

void Graph::find_signed_kplex()
{
	Timer t;
	t.restart();

    lb = 0, ub = n;

	// find heuristic signed k-plex
	CTCP(-1, 1, lb+1-K, lb+1-2*K);
	heu_signed_kplex(10, K);
	// return;

    if((int)kplex.size() < ub) {
        lb = max((int)kplex.size(), 2*K-2);
        get_k_core(lb+1-K);
		get_degree();
		get_tricnt();
		CTCP(-1, 1, lb+1-K, lb+1-2*K);
		printf("---------------------------\n");

		KPLEX_MATRIX *kplex_solver = new KPLEX_MATRIX();
		kplex_solver->allocateMemory(n, m/2);

        vector<pair<int,int> > vp; vp.reserve(m/2);
		vector<int> sgn; sgn.reserve(m/2);

        while(n > 0) {
            //get u
			get_degree();
            ui u = 0;
            for(ui i = 1; i < n; i++) if(degree[u] > degree[i]) u = i;
            
            //get g
			vp.clear();
			sgn.clear();
			get_g(u, vp, sgn);

            //kplex
			kplex_solver->load_graph(s_n, vp, sgn);
			kplex_solver->kPlex(K, kplex, true);
			
			bool lb_changed = 0;

			if(kplex.size() > lb) {
				if(kplex.size() > lb) lb_changed = 1;
				lb = kplex.size();
			}
#ifndef NDEBUG
			printf("\t kplex = %d \n", (int)kplex.size());
#endif
            // CTCP
            CTCP(u, lb_changed, lb+1-K, lb+1-2*K);
        }
    }
    cout<<"\t find_signed_kplex, T : "<<integer_to_string(t.elapsed())<<",\t kplex_size="<<kplex.size()<<endl;
}