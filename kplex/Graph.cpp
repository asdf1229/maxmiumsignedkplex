#include "Graph.h"
#include "Utility.h"
#include "KPlex_matrix.h"

using namespace std;

Graph::Graph(const int _K)
{
    kplex.clear();
    K = _K;

    n = m = pm = nm = 0;

    pstart = NULL;
    pend = NULL;
    edges = NULL;

    p_pstart = NULL;
    p_pend = NULL;
    p_edges = NULL;

    n_pstart = NULL;
    n_pend = NULL;
    n_edges = NULL;

	degree = NULL;
    p_degree = NULL;
    n_degree = NULL;
    tri_cnt = NULL;

    lb = 0, ub = 1e9;
	s_n = 0;
}

Graph::~Graph()
{
	if(pstart != NULL) {
		delete[] pstart;
		pstart = NULL;
	}
	if(pend != NULL) {
		delete[] pend;
		pend = NULL;
	}
	if(edges != NULL) {
		delete[] edges;
		edges = NULL;
	}
	if(p_pstart != NULL) {
		delete[] p_pstart;
		p_pstart = NULL;
	}
	if(p_pend != NULL) {
		delete[] p_pend;
		p_pend = NULL;
	}
	if(p_edges != NULL) {
		delete[] p_edges;
		p_edges = NULL;
	}
	if(n_pstart != NULL) {
		delete[] n_pstart;
		n_pstart = NULL;
	}
	if(n_pend != NULL) {
		delete[] n_pend;
		n_pend = NULL;
	}
	if(n_edges != NULL) {
		delete[] n_edges;
		n_edges = NULL;
	}
	if(degree != NULL) {
		delete[] degree;
		degree = NULL;
	}
	if(p_degree != NULL) {
		delete[] p_degree;
		p_degree = NULL;
	}
	if(n_degree != NULL) {
		delete[] n_degree;
		n_degree = NULL;
	}
	if(tri_cnt != NULL) {
		delete[] tri_cnt;
		tri_cnt = NULL;
	}
}

//Read graph and remove self-loops and multiple edges.
void Graph::load_graph(string input_graph)
{
	Timer t;
	t.restart();
    string buffer;
    ifstream input_file(input_graph, ios::in);

    if(!input_file.is_open()) {
        cout<<"cannot open file : "<<input_graph<<endl; exit(1);
    }
    else {
		input_file >> n >> m;
		map<ui, int> *s_G = new map<ui, int>[n];
		ui u, v; int flag;
		while(input_file >> u >> v >> flag) {
			if(u == v) continue;
            assert(u >= 0 && u < n);
            assert(v >= 0 && v < n);
            assert(flag == 1 || flag == -1);
            s_G[u].insert(make_pair(v, flag));
            s_G[v].insert(make_pair(u, flag));
		}

		m = pm = nm = 0;
		for(ui i = 0; i < n; i++) {
			const map<ui, int> &nei = s_G[i];
			for(auto e : nei) {
				if(e.second == 1) ++pm;
				else ++nm;
			}
			m += nei.size();
		}
		assert(m%2 == 0); assert(pm%2 == 0); assert(nm%2 == 0);
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
		for(ui i = 0; i < n; i++) {
			const map<ui, int> &nei = s_G[i];
			ui idx = p_pstart[i];
			for(auto e : nei) if(e.second == 1) {
				p_edges[idx++] = e.first;
			}
			p_pend[i] = p_pstart[i+1] = idx;
		}
		assert(p_pstart[n] == 2*pm);

        //construct negative edges
        n_pstart[0] = 0;
		for(ui i = 0; i < n; i++) {
			const map<ui, int> &nei = s_G[i];
			ui idx = n_pstart[i];
			for(auto e : nei) if(e.second == -1) {
				n_edges[idx++] = e.first;
			}
			n_pend[i] = n_pstart[i+1] = idx;
		}
		assert(n_pstart[n] == 2*nm);

        //construct edges
        pstart[0] = 0;
		for(ui i = 0; i < n; i++) {
			const map<ui, int> &nei = s_G[i];
			ui idx = pstart[i];
			for(auto e : nei) {
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

//get degree for all nodes of G
void Graph::get_degree()
{
	for(ui i = 0; i < n; i++) {
        p_degree[i] = p_pend[i] - p_pstart[i];
        n_degree[i] = n_pend[i] - n_pstart[i];
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
		for(ept i = pstart[u]; i < pend[u]; i++) mark[edges[i]] = i+1;

		for(ept i = pstart[u]; i < pend[u]; i++) {
			ui v = edges[i];
			if(u < v) {
				for(ept j = pstart[v]; j < pend[v]; j++) {
					ui w = edges[j];
					if(mark[w] && v < w) {
						ept id_uv = i;
						ept id_vu = pstart[v] + find(edges + pstart[v], edges + pend[v], u);
						ept id_vw = j;
						ept id_wv = pstart[w] + find(edges + pstart[w], edges + pend[w], v);
						ept id_uw = mark[w]-1;
						ept id_wu = pstart[w] + find(edges + pstart[w], edges + pend[w], u);
#ifndef NDEBUG
						ept id_uv1 = pstart[u] + find(edges + pstart[u], edges + pend[u], v);
						ept id_vw1 = pstart[v] + find(edges + pstart[v], edges + pend[v], w);
						ept id_uw1 = pstart[u] + find(edges + pstart[u], edges + pend[u], w);
						assert(id_uv == id_uv1);
						assert(id_vw == id_vw1);
						assert(id_uw == id_uw1);
#endif
						tri_cnt[id_uv]++; tri_cnt[id_vu]++;
						tri_cnt[id_vw]++; tri_cnt[id_wv]++;
						tri_cnt[id_uw]++; tri_cnt[id_wu]++;
					}
				}		
			}
		}

		for(ept i = pstart[u]; i < pend[u]; i++) mark[edges[i]] = 0;
	}

	delete [] mark;
}

ui Graph::get_g(ui u, vector<pair<int,int> > &vp, vector<int> &sgn)
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
	ui rid_u = n;
    ui *rid = new ui[n];
	ui cnt = 0;
	for(ui i = 0; i < n; i++) if(v_sel[i]) {
		if(u == i) rid_u = cnt;
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
	assert(rid_u != n);
	return rid_u;
}

void Graph::rebuild_graph(bool *v_del, bool *e_del)
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
			ept p_pointer = p_pstart[u], n_pointer = n_pstart[u];
			for(ept i = pstart[u]; i < pend[u]; i++) if(!e_del[i]) {
				ui v = edges[i];
				if(!v_del[v]) {
					edges[pos++] = rid[v];

					while(p_pointer < p_pend[u] && p_edges[p_pointer] < v) p_pointer++;
					while(n_pointer < n_pend[u] && n_edges[n_pointer] < v) n_pointer++;

					if(p_pointer < p_pend[u] && p_edges[p_pointer] == v) {
						p_edges[p_pos++] = rid[v];
						p_pointer++;
					}
					else if(n_pointer < n_pend[u] && n_edges[n_pointer] == v) {
						n_edges[n_pos++] = rid[v];
						n_pointer++;
					}
					else {
						throw std::runtime_error("rebuild_graph error");
					}
				}
			}
			pend[cnt] = pos;
			p_pend[cnt] = p_pos;
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

//*
void Graph::CTCP(int del_v, int tv, int te)
{
	static int last_tv = 0;
    Timer t;
	t.restart();

    // printf("\t CTCP: tv = %d, te = %d\n", tv, te);
	tv = max(0, tv); te = max(0, te);
	assert(last_tv <= tv);

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
    if(last_tv < tv) {
        for(ui u = 0; u < n; u++) {
            if(degree[u] < tv) qv.push(u);
            for(ept i = pstart[u]; i < pend[u]; i++) {
                ui v = edges[i];
                if(u < v && tri_cnt[i] < te) {
                    qe.push(make_pair(u, i));
                }
            }
        }
    }
	last_tv = tv;

	while(!qv.empty() || !qe.empty()) {
		while(!qe.empty()) {
			auto ue = qe.front(); qe.pop();
			ui u = ue.first; ept id_uv = ue.second; ui v = edges[id_uv];
			ept id_vu = pstart[v] + find(edges + pstart[v], edges + pend[v], u);
			assert(e_del[id_uv] == e_del[id_vu]);
			if(v_del[u] || v_del[v] || e_del[id_uv]) continue;
			e_del[id_uv] = 1;
			e_del[id_vu] = 1;

			if((degree[u]--) == tv) qv.push(u);
			if((degree[v]--) == tv) qv.push(v);

			for(ept i = pstart[u]; i < pend[u]; i++) if(!e_del[i]) mark[edges[i]] = i+1;

			for(ept j = pstart[v]; j < pend[v]; j++) if(!e_del[j]) {
				ui w = edges[j];
				if(mark[w]) { // triangle count--
					ept id_uw = mark[w] - 1;
					ept id_wu = pstart[w] + find(edges + pstart[w], edges + pend[w], u);
					if((tri_cnt[id_uw]--) == te) qe.push(make_pair(u, id_uw));
					tri_cnt[id_wu]--;

					ept id_vw = j;
					ept id_wv = pstart[w] + find(edges + pstart[w], edges + pend[w], v);
					if((tri_cnt[id_vw]--) == te) qe.push(make_pair(v, id_vw));
					tri_cnt[id_wv]--;
				}
			}

			for(ept i = pstart[u]; i < pend[u]; i++) if(!e_del[i]) mark[edges[i]] = 0;
		}
		if(!qv.empty()) {
			ui u = qv.front(); qv.pop();
			if(v_del[u]) continue;
			v_del[u] = 1;

			for(ept i = pstart[u]; i < pend[u]; i++) if(!e_del[i]) mark[edges[i]] = i+1;

			for(ept i = pstart[u]; i < pend[u]; i++) if(!e_del[i]) {
				ui v = edges[i];
				for(ept j = pstart[v]; j < pend[v]; j++) if(!e_del[j]) {
					ui w = edges[j];
					if(mark[w] && v < w) { // triangle count--
						ept id_vw = j;
						ept id_wv = pstart[w] + find(edges + pstart[w], edges + pend[w], v);
						if((tri_cnt[id_vw]--) == te) qe.push(make_pair(v, id_vw));
						tri_cnt[id_wv]--;
					}
				}
			}

			for(ept i = pstart[u]; i < pend[u]; i++) if(!e_del[i]) mark[edges[i]] = 0;

			for(ept i = pstart[u]; i < pend[u]; i++) if(!e_del[i]) {
				ui v = edges[i];
				if((degree[v]--) == tv) qv.push(v);
				ept id_uv = i;
				ept id_vu = pstart[v] + find(edges + pstart[v], edges + pend[v], u);
				e_del[id_uv] = 1;
				e_del[id_vu] = 1;
			}
		}
	}

	rebuild_graph(v_del, e_del);

	delete[] v_del;
	delete[] e_del;

#ifndef NDEBUG
	get_degree();
	get_tricnt();
	for(ui u = 0; u < n; u++) {
		assert(degree[u] >= tv);
		for(ept i = pstart[u]; i < pend[u]; i++) {
			assert(tri_cnt[i] >= te);
		}
	}
#endif

	// #ifndef NDEBUG
	cout<<"\t CTCP, T : "<<integer_to_string(t.elapsed())<<",\t n = "<<n<<", m = "<<m<<endl;
	// #endif
}

//*
void Graph::heu_signed_kplex(int rounds, int k)
{
// 	Timer t; t.restart();
// 	if(rounds < 1) return;

// 	KPLEX_MATRIX *kplex_solver = new KPLEX_MATRIX();
// 	kplex_solver->allocateMemory(n, m/2);

// 	vector<pair<int,int> > vp; vp.reserve(m/2);
// 	vector<int> sgn; sgn.reserve(m/2);

// 	for(int round = 0; round < rounds && round < n; round++) {
// 		get_degree();
// 		ui u = 0;
// 		for(ui i = 1; i < n; i++) if(degree[u] > degree[i]) u = i;
// 		// printf("%d %d\n", u, degree[u]);

// 		vp.clear();
// 		sgn.clear();
// 		get_g(u, vp, sgn);

//         //heuristically find signed kplex
// 		kplex_solver->load_graph(s_n, vp, sgn);
// 		// kplex_solver->heu_kPlex(K, kplex);
// // #ifndef NDEBUG
// 		// printf("\t heu_kplex = %d \n", (int)kplex.size());
// 		// vector<ui> KPLEX; KPLEX.clear();
// 		kplex_solver->kPlex(K, kplex);
// 		printf("\t KPLEX = %d \n", (int)kplex.size());
// // #endif
// 		bool lb_changed = 0;
// 		if(kplex.size() > lb) {
// 			if(kplex.size() > lb) lb_changed = 1;
// 			lb = kplex.size();
// 		}
//         // CTCP
// 		CTCP(u, lb_changed, lb+1-K, lb+1-2*K);
//         // CTCP(-1, lb_changed, lb+1-K, lb+1-2*K);
// 	}
// 	cout<<"\t heu_find_kplex, T : "<<integer_to_string(t.elapsed())<<",\t heu_kplex_size = "<<kplex.size()<<endl;
}

void Graph::find_signed_kplex()
{
	Timer t;
	t.restart();

    lb = 2*K-2, ub = n;

	// find heuristic signed k-plex
	// CTCP(-1, 1, lb+1-K, lb+1-2*K);
	// heu_signed_kplex(10, K);
	// return;

	lb = max((int)kplex.size(), lb);
	CTCP(-1, lb+1-K, lb+1-2*K);
	printf("-----------------------------------\n");

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
		ui rid_u = get_g(u, vp, sgn);
		assert(0 <= rid_u && rid_u < n);

		//kplex
		kplex_solver->load_graph(s_n, vp, sgn);
		kplex_solver->kPlex(K, kplex, (s_n == n) ? -1 : rid_u);

		if(kplex.size() > lb) lb = kplex.size();
#ifndef NDEBUG
		printf("\t kplex = %d \n", (int)kplex.size());
#endif
		if(s_n == n) break;
		// CTCP
		CTCP(u, lb+1-K, lb+1-2*K);
	}
    cout<<"\t find_signed_kplex, T : "<<integer_to_string(t.elapsed())<<",\t kplex_size="<<kplex.size()<<endl;
}