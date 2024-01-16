// degeneracy-based k-plex
// return an upper bound of the maximum k-plex size
// return dOrder
void kplex_degen(ListLinearHeap *heap, int k, int *dOrder)
{
	Timer t;
	int *peel_sequence = new int[G.n];
	int *vis = new int[G.n];
	
	for(int i = 0; i < G.n; i++) peel_sequence[i] = i;
	for(int i = 0; i < G.n; i++) vis[i] = 0;

	heap->init(G.n, G.n-1, peel_sequence, G.degree);
	for(int i = 0; i < G.n; i++) {
		int u, deg;
		heap->pop_min(u, deg);
		dOrder[i] = u;

		// if(deg+k >= G.n-i+1 && G.n-i+1 > lb) lb = G.n-i+1;
		ub = max(ub, min(deg+k, G.n-i+1));

		for(int j = G.pstart[u]; j < G.pend[u]; j++) {
			int v = G.edges[j];
			if(vis[v] == 0) heap->decrement(v, 1);
		}
		vis[u] = 1;
	}

	// for(int i = 0; i < G.n; i++) {
	// 	printf("%d ", dOrder[i]);
	// }
	// printf("\n");

	delete [] peel_sequence;
	delete [] vis;

}

//get degree for all nodes of G
void get_G_deg()
{
	// for (int i = 0; i < G.n; i++) {
    //     G.p_degree[i] = G.p_pstart[i+1] - G.p_pstart[i];
    //     G.n_degree[i] = G.n_pstart[i+1] - G.n_pstart[i];
    //     G.degree[i] = G.p_degree[i] + G.n_degree[i];
	// }
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

int Graph::cal(int k, int siz, vector<ui> s)
{
	int * deg = new int[siz];
	memset(deg, 0, sizeof(int)*siz);

    for(ui i = 0; i < siz; i++) {
        for(ui j = i + 1; j < siz; j++) {
            ui u = s[i], v = s[j];
            if(s_matrix[u][v]) {
                deg[i]++;
                deg[j]++;
            }
        }
    }
    for(ui i = 0; i < siz; i++) {
        if(deg[i] < siz - k) return 0;
    }
	delete [] deg;

	for(ui i = 0; i < siz; i++) {
		ui u = s[i];
		for(ui j = i + 1; j < siz; j++) {
			ui v = s[j];
			if(s_matrix[u][v]) {
				for(ui k = j + 1; k < siz; k++) {
					ui w = s[k];
					if(s_matrix[v][w] && s_matrix[u][w]) {
						int tmp = s_matrix[u][v] + s_matrix[v][w] + s_matrix[u][w];
						if(tmp == 1 || tmp == -3) {
							return 0;
						}
					}
				}
			}
		}
	}
    return 1;
}