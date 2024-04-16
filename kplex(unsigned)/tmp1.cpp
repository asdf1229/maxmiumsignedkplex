//get k-core
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
