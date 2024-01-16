#include "Graph.h"
#include "Utility.h"

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

    // s_n = 0;
    // s_matrix = nullptr;
    // s_solution_size = 0;
    // s_solution = nullptr;

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
    // for(int i = 0; i < max_core; i++) if(s_matrix[i] != nullptr){
    //     delete [] s_matrix[i];
    //     s_matrix[i] = nullptr;
    // }
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
		assert(p_pstart[n] == 2*nm);

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

//get G's k-core
void Graph::get_k_core(int k)
{
	printf("\t k = %d\n", k);
	Timer t;
	t.restart();

    get_degree();

	bool *v_del = new bool[n];
	memset(v_del, 0, sizeof(bool)*n);

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

    cout<<"\t get_k_core, T : "<<integer_to_string(t.elapsed())<<",\t n="<<n<<", m="<<m<<endl;
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

void Graph::get_g(ui u)
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
    
	for(ui i = 0; i < s_n; i++)
		for(ui j = 0; j < s_n; j++)
			s_matrix[i][j] = 0;
    
    for(ui u = 0; u < n; u++) if(v_sel[u]) {
        for(ept i = p_pstart[u]; i < p_pend[u]; i++) {
            ui v = p_edges[i];
            if(v_sel[v]) {
                s_matrix[rid[u]][rid[v]] = 1;
            }
        }
        for(ept i = n_pstart[u]; i < n_pend[u]; i++) {
            ui v = n_edges[i];
            if(v_sel[v]) {
                s_matrix[rid[u]][rid[v]] = -1;
            }
        }
    }

	// for(int i = 0; i < s_n; i++) {
	// 	for(int j = 0; j < s_n; j++) {
	// 		printf("%d ", Matrix[i][j]);
	// 	}
	// 	printf("\n");
	// }

	delete [] rid;
	delete [] v_sel;

	cout<<"\t get_g, T : "<<integer_to_string(t.elapsed())<<",\t n="<<n<<", m="<<m<<endl;
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
			pend[cnt] = pstart[cnt+1] = pos;

			for(ept i = p_pstart[u]; i < p_pend[u]; i++) {
				ui v = p_edges[i];
				if(!v_del[v]) {
					p_edges[p_pos++] = rid[v];
				}
			}
			p_pend[cnt] = p_pstart[cnt+1] = p_pos;

			for(ept i = n_pstart[u]; i < n_pend[u]; i++) {
				ui v = n_edges[i];
				if(!v_del[v]) {
					n_edges[n_pos++] = rid[v];
				}
			}
			n_pend[cnt] = n_pstart[cnt+1] = n_pos;
			cnt++;
		}

		n = cnt;
		m = pos / 2;
		pm = p_pos / 2;
		nm = n_pos / 2;
    }

    delete[] rid;
}

void Graph::CTCP(int del_v, bool lb_changed, ui tv, ui te)
{
    printf("\t CTCP: tv = %d, te = %d\n", tv, te);
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
                    printf("%d %d\n", u, v);
                    qe.push(make_pair(u, i));
                }
            }
        }
    }

    while(!qv.empty() || !qe.empty()) {
        while(!qe.empty()) {
			auto ue = qe.front(); qe.pop();
			ui u = ue.first; ept idx = ue.second; ui v = edges[idx];
			if(v_del[u] || v_del[v] || e_del[idx]) continue;
            e_del[idx] = 1;

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

    cout<<"\t CTCP, T : "<<integer_to_string(t.elapsed())<<",\t n="<<n<<", m="<<m<<endl;
}

void Graph::heu_signed_kplex(int rounds, int k)
{
	if(rounds < 1) return;
	priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> kset;
	for(int i = 0; i < rounds; i++) kset.push(make_pair(degree[i], i));

    for(int i = rounds; i < n; i++){
        if(degree[i] > kset.top().first){
            kset.pop();
            kset.push(make_pair(degree[i], i));
        }
    }
    vector<pair<int, int>> ordV(rounds);
    for(int i = 0; i < rounds; i++){
        ordV[i] = kset.top();
        kset.pop();
    }
    assert(kset.empty());
	sort(ordV.begin(), ordV.end(), greater<pair<int, int>>()); //decreasing order

    int * label = new int[n];
    int * vs_deg = new int[n];
	int * inR = new int[n];
	
	for(int round = 0; round < rounds && round < n; round++) {
        int u = ordV[round].second;

		get_g(u);
		
        memset(label, 0, sizeof(int)*n);
		memset(inR, 0, sizeof(int)*n);
        vector<ui> res;
        res.push_back(u);
		inR[u] = 1;
        vector<int> vsP, vsN;
        for(int i = p_pstart[u]; i < p_pend[u]; i++){
            int v = p_edges[i];
            vsP.push_back(v);
            label[v] = 1;
        }
        for(int i = n_pstart[u]; i < n_pend[u]; i++){
            int v = n_edges[i];
            vsN.push_back(v);
            label[v] = 2;
        }
        for(auto e : vsP) vs_deg[e] = 0;
        for(auto e : vsN) vs_deg[e] = 0;
        for(auto e : vsP){
            for(int i = p_pstart[e]; i < p_pend[e]; i++){
                int v = p_edges[i];
                if(label[v] == 1) ++ vs_deg[e];
            }
            for(int i = n_pstart[e]; i < n_pend[e]; i++){
                int v = n_edges[i];
                if(label[v] == 2) ++ vs_deg[e];
            }
        }
        for(auto e : vsN){
            for(int i = p_pstart[e]; i < p_pend[e]; i++){
                int v = p_edges[i];
                if(label[v] == 2) ++ vs_deg[e];
            }
            for(int i = n_pstart[e]; i < n_pend[e]; i++){
                int v = n_edges[i];
                if(label[v] == 1) ++ vs_deg[e];
            }
        }

		while (!vsP.empty() || !vsN.empty()) {
			if(!vsP.empty()) {
                int tmp_deg = 0;
                int next_v;
                for(int i = 0; i < vsP.size(); i++){
                    if(vs_deg[vsP[i]] >= tmp_deg){
                        tmp_deg = vs_deg[vsP[i]];
                        next_v = vsP[i];
                    }
                }
                res.push_back(next_v);
				inR[next_v] = 1;
                vector<int> new_vsP, new_vsN;
                assert(label[next_v] == 1);
                for(int i = p_pstart[next_v]; i < p_pend[next_v]; i++){
                    int v = p_edges[i];
                    if(label[v] == 1) new_vsP.push_back(v);
                }
                for(int i = n_pstart[next_v]; i < n_pend[next_v]; i++){
                    int v = n_edges[i];
                    if(label[v] == 2) new_vsN.push_back(v);
                }
                for(auto e : vsP) label[e] = 0;
                for(auto e : vsN) label[e] = 0;
                vsP = new_vsP;
                vsN = new_vsN;
                for(auto e : vsP) label[e] = 1;
                for(auto e : vsN) label[e] = 2;
                for(auto e : vsP) vs_deg[e] = 0;
                for(auto e : vsN) vs_deg[e] = 0;
                for(auto e : vsP){
                    for(int i = p_pstart[e]; i < p_pend[e]; i++){
                        int v = p_edges[i];
                        if(label[v] == 1) ++ vs_deg[e];
                    }
                    for(int i = n_pstart[e]; i < n_pend[e]; i++){
                        int v = n_edges[i];
                        if(label[v] == 2) ++ vs_deg[e];
                    }
                }
                for(auto e : vsN){
                    for(int i = p_pstart[e]; i < p_pend[e]; i++){
                        int v = p_edges[i];
                        if(label[v] == 2) ++ vs_deg[e];
                    }
                    for(int i = n_pstart[e]; i < n_pend[e]; i++){
                        int v = n_edges[i];
                        if(label[v] == 1) ++ vs_deg[e];
                    }
                }
			}
			else if(!vsN.empty()) {
                int tmp_deg = 0;
                int next_v;
                for(int i = 0; i < vsN.size(); i++){
                    if(vs_deg[vsN[i]] >= tmp_deg){
                        tmp_deg = vs_deg[vsN[i]];
                        next_v = vsN[i];
                    }
                }
                res.push_back(next_v);
				inR[next_v] = 1;
                vector<int> new_vsP, new_vsN;
                assert(label[next_v] == 2);
                for(int i = p_pstart[next_v]; i < p_pend[next_v]; i++){
                    int v = p_edges[i];
                    if(label[v] == 2) new_vsN.push_back(v);
                }
                for(int i = n_pstart[next_v]; i < n_pend[next_v]; i++){
                    int v = n_edges[i];
                    if(label[v] == 1) new_vsP.push_back(v);
                }
                for(auto e : vsP) label[e] = 0;
                for(auto e : vsN) label[e] = 0;
                vsP = new_vsP;
                vsN = new_vsN;
                for(auto e : vsP) label[e] = 1;
                for(auto e : vsN) label[e] = 2;
                for(auto e : vsP) vs_deg[e] = 0;
                for(auto e : vsN) vs_deg[e] = 0;
                for(auto e : vsP){
                    for(int i = p_pstart[e]; i < p_pend[e]; i++){
                        int v = p_edges[i];
                        if(label[v] == 1) ++ vs_deg[e];
                    }
                    for(int i = n_pstart[e]; i < n_pend[e]; i++){
                        int v = n_edges[i];
                        if(label[v] == 2) ++ vs_deg[e];
                    }
                }
                for(auto e : vsN){
                    for(int i = p_pstart[e]; i < p_pend[e]; i++){
                        int v = p_edges[i];
                        if(label[v] == 2) ++ vs_deg[e];
                    }
                    for(int i = n_pstart[e]; i < n_pend[e]; i++){
                        int v = n_edges[i];
                        if(label[v] == 1) ++ vs_deg[e];
                    }
                }
			}
		}

		// cout << "clique_size = " << res.size() << endl;
		if(res.size() > kplex.size()) {
			kplex = res;
		}

		ordV.clear();
		for(int i = 0; i < n; i++) {
			if(inR[i]) continue;
			ordV.push_back(make_pair(degree[i], i));
			sort(ordV.begin(), ordV.end(), greater<pair<int, int>>()); //decreasing order
		}
		for(auto pa : ordV) {
			int v = pa.second;
			// printf("%d ", v);
			res.push_back(v);
			if(!cal(k, res.size(), res)) res.pop_back();
		}

		// cout << "kplex_size = " << res.size() << endl;
		if(res.size() > kplex.size()) {
			kplex = res;
		}
	}
	cout << "\t hec_kplex_size = " << kplex.size() << endl;
	// for(auto u : P) {
	// 	printf("%d ", u);
	// }
	// printf("\n");
}

//1
bool Graph::move_u_to_S_with_prune(ui S_end, ui &R_end, ui level) {
    assert(S_end > 0);
    ui u = SR[S_end-1];

    ui neighbors_n = 0, nonneighbors_n = 0;
    for(ui i = 0;i < R_end;i ++) if(i != S_end-1) {
        if(s_matrix[u][SR[i]]) neighbors[neighbors_n++] = SR[i];
        else nonneighbors[nonneighbors_n++] = SR[i];
    }

    for(ui i = 0;i < neighbors_n;i ++) ++s_degree_in_S[neighbors[i]];

    bool is_balanced = true;
    //
    if(is_balanced == false) return true;

    assert(Qv.empty());
    if(s_degree_in_S[u] + K == S_end) { // only neighbors of u in R can be candidates --- RR2
        ui i = 0;
        while(i < nonneighbors_n&&SR_rid[nonneighbors[i]] < S_end) ++ i;
        for(;i < nonneighbors_n;i ++) { // remove non-neighbors from R
            assert(level_id[nonneighbors[i]] > level);
            level_id[nonneighbors[i]] = level;
            Qv.push(nonneighbors[i]);
        }
    }
    else { // only non-neighbors of u may change their allowance --- RR1
        ui i = 0;
        while(i < nonneighbors_n&&SR_rid[nonneighbors[i]] < S_end) ++ i;
        for(;i < nonneighbors_n;i ++) if(s_degree_in_S[nonneighbors[i]] + K <= S_end) {
            assert(level_id[nonneighbors[i]] > level);
            level_id[nonneighbors[i]] = level;
            Qv.push(nonneighbors[i]);
        }
    }

    // RR2
    for(ui i = 0;i < nonneighbors_n&&SR_rid[nonneighbors[i]] < S_end;i ++) if(s_degree_in_S[nonneighbors[i]] + K == S_end) {
        ui v = nonneighbors[i];
        for(ui j = S_end;j < R_end;j ++) if(level_id[SR[j]] > level&&!s_matrix[v][SR[j]]) {
            level_id[SR[j]] = level;
            Qv.push(SR[j]);
        }
    }
    return remove_vertices_and_edges_with_prune(S_end, R_end, level);
}

//ok
bool Graph::remove_vertices_and_edges_with_prune(ui S_end, ui &R_end, ui level) {
    while(!Qv.empty()) {
        while(!Qv.empty()) {
            ui u = Qv.front(); Qv.pop(); // remove u
            assert(SR[SR_rid[u]] == u);
            assert(SR_rid[u] >= S_end&&SR_rid[u] < R_end);
            -- R_end;
            swap_pos(SR_rid[u], R_end);

            bool terminate = false;
            ui neighbors_n = 0;
            for(ui i = 0;i < R_end;i ++) if(s_matrix[u][SR[i]]) {
                ui w = SR[i];
                neighbors[neighbors_n++] = w;
                -- s_degree[w];
                if(s_degree[w] + K <= s_solution_size) {
                    if(i < S_end) terminate = true; // UB1
                    else if(level_id[w] > level) { // RR3
                        level_id[w] = level;
                        Qv.push(w);
                    }
                }
            }
            // UB1
            if(terminate) {
                for(ui i = 0;i < neighbors_n;i ++) ++ s_degree[neighbors[i]];
                level_id[u] = s_n;
                ++ R_end;
                return true;
            }
        }
    }
    return false;
}

//ok
void Graph::restore_SR_and_edges(ui S_end, ui &R_end, ui old_R_end, ui level) {
    while(!Qv.empty()) {
        ui u = Qv.front(); Qv.pop();
        assert(level_id[u] == level&&SR_rid[u] < R_end);
        level_id[u] = s_n;
    }
    while(R_end < old_R_end) { // insert u back into R
        ui u = SR[R_end];
        assert(level_id[u] == level&&SR_rid[u] == R_end);
        level_id[u] = s_n;

        ui neighbors_n = 0;
        for(ui i = 0;i < R_end;i ++) if(s_matrix[u][SR[i]]) {
            ui w = SR[i];
            neighbors[neighbors_n ++] = w;
            ++ s_degree[w];
        }
        ++ R_end;
    }
}

//ok
void Graph::move_u_to_R_wo_prune(ui &S_end, ui &R_end, ui level) {
    assert(S_end);
    ui u = SR[--S_end];
    ui neighbors_n = 0;
    for(ui i = 0;i < R_end;i++) if(s_matrix[u][SR[i]]) neighbors[neighbors_n++] = SR[i];
    for(ui i = 0;i < neighbors_n;i++) --s_degree_in_S[neighbors[i]];
}

//ok
bool Graph::remove_u_from_S_with_prune(ui &S_end, ui &R_end, ui level) {
    assert(S_end);
    ui u = SR[S_end-1];
    -- S_end; -- R_end;
    swap_pos(S_end, R_end);
    level_id[u] = level;

    bool ret = false;
    ui neighbors_n = 0;
    for(ui i = 0;i < R_end;i ++) if(s_matrix[u][SR[i]]) neighbors[neighbors_n ++] = SR[i];
    for(ui i = 0;i < neighbors_n;i ++) {
        -- s_degree_in_S[neighbors[i]];
        -- s_degree[neighbors[i]];
        if(s_degree[neighbors[i]] + K <= s_solution_size) {
            if(SR_rid[neighbors[i]] < S_end) ret =  true;
            else {
                assert(level_id[neighbors[i]] > level);
                level_id[neighbors[i]] = level;
                Qv.push(neighbors[i]);
            }
        }
    }
    if(ret) return true;
    return false;
}

//ok
bool Graph::collect_removable_vertices_and_edges(ui S_end, ui R_end, ui level) {
    for(ui i = 0;i < S_end;i ++) if(s_degree[SR[i]] + K <= s_solution_size) return true;

    for(ui i = S_end;i < R_end;i ++) if(level_id[SR[i]] > level){
        if(S_end - s_degree_in_S[SR[i]] >= K||s_degree[SR[i]] + K <= s_solution_size) {
            assert(level_id[SR[i]] > level);
            level_id[SR[i]] = level;
            Qv.push(SR[i]);
            continue;
        }
        for(ui j = 0;j < S_end;j ++) {
            if(S_end - s_degree_in_S[SR[j]] == K&&!s_matrix[SR[i]][SR[j]])
            {
                assert(level_id[SR[i]] > level);
                level_id[SR[i]] = level;
                Qv.push(SR[i]);
                break;
            }
        }
    }

    return false;
}

//ok
ui Graph::choose_branch_vertex_based_on_non_neighbors(ui S_end, ui R_end) {
    assert(SR[S_end] != s_n);
    return SR[S_end];
    // ui u = s_n, min_degree_in_S = s_n;
    // for(ui i = 0;i < R_end;i ++) {
    //     ui v = SR[i];
    //     if(s_degree[v] + K >= R_end) continue;
    //     if(s_degree_in_S[v] < min_degree_in_S) {
    //         u = v;
    //         min_degree_in_S = s_degree_in_S[v];
    //     }
    //     else if(s_degree_in_S[v] == min_degree_in_S) {
    //         if((i < S_end&&s_degree[v] < s_degree[u])||(i >= S_end&&s_degree[v] > s_degree[u])) {
    //             u = v;
    //             min_degree_in_S = s_degree_in_S[v];
    //         }
    //     }
    // }
    // assert(u != s_n);
    // if(SR_rid[u] < S_end) {
    //     u = s_n;
    //     ui max_degree = 0;
    //     for(ui i = S_end;i < R_end;i ++) if(!s_matrix[u][SR[i]]&&s_degree[SR[i]] > max_degree) {
    //         max_degree = s_degree[SR[i]];
    //         u = SR[i];
    //     }
    // }
    // assert(u != s_n);
    // return u;
}

//ok
void Graph::kplex_search(ui S_end, ui R_end, ui level, bool choose_zero)
{
    if(S_end > s_solution_size) { // find a larger solution
        s_solution_size = S_end;
        for(ui i = 0; i < s_solution_size; i++) s_solution[i] = SR[i];
    }
    if(R_end <= s_solution_size) return;

    ui u = s_n; // u is the branching vertex
    bool must_include = false;
    if(choose_zero) {
        assert(S_end == 0);
        if(SR_rid[0] >= R_end) return;
        u = 0;
        must_include = true;
    }
    else {
        if(u == s_n) {
            u = choose_branch_vertex_based_on_non_neighbors(S_end, R_end);
        }
    }
    assert(s_degree[u] + K > s_solution_size&&s_degree[u] + K > S_end);

    // the first branch includes u into S
    assert(SR[SR_rid[u]] == u&&SR_rid[u] >= S_end&&SR_rid[u] < R_end);
    swap_pos(S_end, SR_rid[u]);
    ++S_end;

    ui pre_s_solution_size = s_solution_size, old_R_end = R_end;
    if(!move_u_to_S_with_prune(S_end, R_end, level)) kplex_search(S_end, R_end, level+1, false);
    restore_SR_and_edges(S_end, R_end, old_R_end, level);

    if(must_include) {
        move_u_to_R_wo_prune(S_end, R_end, level);
        return;
    }

    // the second branch exclude u from S
    assert(Qv.empty());
    bool pruned = remove_u_from_S_with_prune(S_end, R_end, level);
    if(!pruned && s_solution_size > pre_s_solution_size) pruned = collect_removable_vertices_and_edges(S_end, R_end, level);
    if(!pruned) {
        if(!remove_vertices_and_edges_with_prune(S_end, R_end, level)) {
            kplex_search(S_end, R_end, level+1, false);
        }
        restore_SR_and_edges(S_end, R_end, old_R_end, level);
    }
}

void Graph::s_init(ui &R_end) {
    s_degree = new ui[n];
    s_degree_in_S = new ui[n];
    s_solution = new ui[n];
    SR = new ui[n];
    SR_rid = new ui[n];
    neighbors = new ui[n];
    nonneighbors = new ui[n];
    level_id = new ui[n];

    memset(s_degree_in_S, 0, sizeof(ui)*s_n);
    R_end = 0;
    for(ui i = 0;i < s_n;i ++) SR_rid[i] = s_n;
    for(ui i = 0;i < s_n;i ++) {
        SR[R_end] = i; SR_rid[i] = R_end;
        ++ R_end;
    }

    if(SR_rid[0] == s_n) {
        R_end = 0;
        return;
    }

    for(ui i = 0;i < R_end;i ++) {
        ui u = SR[i];
        s_degree[u] = 0;
        for(ui j = 0;j < R_end;j ++) if(s_matrix[u][SR[j]]) ++ s_degree[u];
    }

    memset(level_id, 0, sizeof(ui)*s_n);
    for(ui i = 0;i < R_end;i ++) level_id[SR[i]] = s_n;

    assert(Qv.empty());

    if(remove_vertices_and_edges_with_prune(0, R_end, 0)) R_end = 0;
}

void Graph::matrix_kplex() {
    s_solution_size = kplex.size();
    ui R_end;
    s_init(R_end);
    if(R_end) kplex_search(0, R_end, 1, true);
    if(s_solution_size > kplex.size()) {
        kplex.clear();
        for(int i = 0;i < s_solution_size;i ++) kplex.push_back(s_solution[i]);
    }
}

void Graph::swap_pos(ui i, ui j) {
    swap(SR[i], SR[j]);
    SR_rid[SR[i]] = i;
    SR_rid[SR[j]] = j;
}

void Graph::subgraph_init() {
	s_matrix = new int*[n];
	for(ui i = 0; i < n; i++) s_matrix[i] = new int[n];

    s_solution = NULL;
    s_degree = s_degree_in_S = NULL;

    SR = SR_rid = NULL;
    level_id = NULL;

    neighbors = nonneighbors = NULL;
}

void Graph::find_signed_kplex() {
	Timer t;
	t.restart();

    lb = 0, ub = n;

    int ans = 0;
	// // find heuristic signed k-plex
	// heu_signed_kplex(1, k);
	// ans = P.size();


    if((int)kplex.size() < ub) {
        lb = max((int)kplex.size(), 2*K-2);
        get_k_core(lb+1-K);
		get_degree();
		get_tricnt();
		CTCP(-1, 1, lb+1-K, max(0, lb+1-2*K));

        subgraph_init();

        while(n > 0) {
            //get u
            ui u = 0;
            for(ui i = 1; i < n; i++) {
				if(degree[u] > degree[i]) {
					u = i;
				}
            }
			printf("u = %d\n", u);
            
            //get g
			get_g(u);

            //kplex
            matrix_kplex();
			bool lb_changed = 0;
			if(kplex.size() >= lb) {
				ans = kplex.size();
				lb = kplex.size();
				if(kplex.size() > lb) lb_changed = 1;
			}
            printf("kplex = %d \n", kplex.size());

            // CTCP
            CTCP(u, lb_changed, lb+1-K, lb+1-2*K);            
        }

    }

    cout<<"\t find_signed_kplex, T : "<<integer_to_string(t.elapsed())<<",\t kplex_size="<<ans<<endl;
}