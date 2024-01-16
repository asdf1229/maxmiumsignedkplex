// void Graph::heu_signed_kplex(int rounds, int k)
// {
// 	if(rounds < 1) return;
// 	priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> kset;
// 	for(int i = 0; i < rounds; i++) kset.push(make_pair(degree[i], i));

//     for(int i = rounds; i < n; i++){
//         if(degree[i] > kset.top().first){
//             kset.pop();
//             kset.push(make_pair(degree[i], i));
//         }
//     }
//     vector<pair<int, int>> ordV(rounds);
//     for(int i = 0; i < rounds; i++){
//         ordV[i] = kset.top();
//         kset.pop();
//     }
//     assert(kset.empty());
// 	sort(ordV.begin(), ordV.end(), greater<pair<int, int>>()); //decreasing order

//     int * label = new int[n];
//     int * vs_deg = new int[n];
// 	int * inR = new int[n];
	
// 	for(int round = 0; round < rounds && round < n; round++) {
//         int u = ordV[round].second;

// 		get_g(u);
		
//         memset(label, 0, sizeof(int)*n);
// 		memset(inR, 0, sizeof(int)*n);
//         vector<ui> res;
//         res.push_back(u);
// 		inR[u] = 1;
//         vector<int> vsP, vsN;
//         for(int i = p_pstart[u]; i < p_pend[u]; i++){
//             int v = p_edges[i];
//             vsP.push_back(v);
//             label[v] = 1;
//         }
//         for(int i = n_pstart[u]; i < n_pend[u]; i++){
//             int v = n_edges[i];
//             vsN.push_back(v);
//             label[v] = 2;
//         }
//         for(auto e : vsP) vs_deg[e] = 0;
//         for(auto e : vsN) vs_deg[e] = 0;
//         for(auto e : vsP){
//             for(int i = p_pstart[e]; i < p_pend[e]; i++){
//                 int v = p_edges[i];
//                 if(label[v] == 1) ++ vs_deg[e];
//             }
//             for(int i = n_pstart[e]; i < n_pend[e]; i++){
//                 int v = n_edges[i];
//                 if(label[v] == 2) ++ vs_deg[e];
//             }
//         }
//         for(auto e : vsN){
//             for(int i = p_pstart[e]; i < p_pend[e]; i++){
//                 int v = p_edges[i];
//                 if(label[v] == 2) ++ vs_deg[e];
//             }
//             for(int i = n_pstart[e]; i < n_pend[e]; i++){
//                 int v = n_edges[i];
//                 if(label[v] == 1) ++ vs_deg[e];
//             }
//         }

// 		while (!vsP.empty() || !vsN.empty()) {
// 			if(!vsP.empty()) {
//                 int tmp_deg = 0;
//                 int next_v;
//                 for(int i = 0; i < vsP.size(); i++){
//                     if(vs_deg[vsP[i]] >= tmp_deg){
//                         tmp_deg = vs_deg[vsP[i]];
//                         next_v = vsP[i];
//                     }
//                 }
//                 res.push_back(next_v);
// 				inR[next_v] = 1;
//                 vector<int> new_vsP, new_vsN;
//                 assert(label[next_v] == 1);
//                 for(int i = p_pstart[next_v]; i < p_pend[next_v]; i++){
//                     int v = p_edges[i];
//                     if(label[v] == 1) new_vsP.push_back(v);
//                 }
//                 for(int i = n_pstart[next_v]; i < n_pend[next_v]; i++){
//                     int v = n_edges[i];
//                     if(label[v] == 2) new_vsN.push_back(v);
//                 }
//                 for(auto e : vsP) label[e] = 0;
//                 for(auto e : vsN) label[e] = 0;
//                 vsP = new_vsP;
//                 vsN = new_vsN;
//                 for(auto e : vsP) label[e] = 1;
//                 for(auto e : vsN) label[e] = 2;
//                 for(auto e : vsP) vs_deg[e] = 0;
//                 for(auto e : vsN) vs_deg[e] = 0;
//                 for(auto e : vsP){
//                     for(int i = p_pstart[e]; i < p_pend[e]; i++){
//                         int v = p_edges[i];
//                         if(label[v] == 1) ++ vs_deg[e];
//                     }
//                     for(int i = n_pstart[e]; i < n_pend[e]; i++){
//                         int v = n_edges[i];
//                         if(label[v] == 2) ++ vs_deg[e];
//                     }
//                 }
//                 for(auto e : vsN){
//                     for(int i = p_pstart[e]; i < p_pend[e]; i++){
//                         int v = p_edges[i];
//                         if(label[v] == 2) ++ vs_deg[e];
//                     }
//                     for(int i = n_pstart[e]; i < n_pend[e]; i++){
//                         int v = n_edges[i];
//                         if(label[v] == 1) ++ vs_deg[e];
//                     }
//                 }
// 			}
// 			else if(!vsN.empty()) {
//                 int tmp_deg = 0;
//                 int next_v;
//                 for(int i = 0; i < vsN.size(); i++){
//                     if(vs_deg[vsN[i]] >= tmp_deg){
//                         tmp_deg = vs_deg[vsN[i]];
//                         next_v = vsN[i];
//                     }
//                 }
//                 res.push_back(next_v);
// 				inR[next_v] = 1;
//                 vector<int> new_vsP, new_vsN;
//                 assert(label[next_v] == 2);
//                 for(int i = p_pstart[next_v]; i < p_pend[next_v]; i++){
//                     int v = p_edges[i];
//                     if(label[v] == 2) new_vsN.push_back(v);
//                 }
//                 for(int i = n_pstart[next_v]; i < n_pend[next_v]; i++){
//                     int v = n_edges[i];
//                     if(label[v] == 1) new_vsP.push_back(v);
//                 }
//                 for(auto e : vsP) label[e] = 0;
//                 for(auto e : vsN) label[e] = 0;
//                 vsP = new_vsP;
//                 vsN = new_vsN;
//                 for(auto e : vsP) label[e] = 1;
//                 for(auto e : vsN) label[e] = 2;
//                 for(auto e : vsP) vs_deg[e] = 0;
//                 for(auto e : vsN) vs_deg[e] = 0;
//                 for(auto e : vsP){
//                     for(int i = p_pstart[e]; i < p_pend[e]; i++){
//                         int v = p_edges[i];
//                         if(label[v] == 1) ++ vs_deg[e];
//                     }
//                     for(int i = n_pstart[e]; i < n_pend[e]; i++){
//                         int v = n_edges[i];
//                         if(label[v] == 2) ++ vs_deg[e];
//                     }
//                 }
//                 for(auto e : vsN){
//                     for(int i = p_pstart[e]; i < p_pend[e]; i++){
//                         int v = p_edges[i];
//                         if(label[v] == 2) ++ vs_deg[e];
//                     }
//                     for(int i = n_pstart[e]; i < n_pend[e]; i++){
//                         int v = n_edges[i];
//                         if(label[v] == 1) ++ vs_deg[e];
//                     }
//                 }
// 			}
// 		}

// 		// cout << "clique_size = " << res.size() << endl;
// 		if(res.size() > kplex.size()) {
// 			kplex = res;
// 		}

// 		ordV.clear();
// 		for(int i = 0; i < n; i++) {
// 			if(inR[i]) continue;
// 			ordV.push_back(make_pair(degree[i], i));
// 			sort(ordV.begin(), ordV.end(), greater<pair<int, int>>()); //decreasing order
// 		}
// 		for(auto pa : ordV) {
// 			int v = pa.second;
// 			// printf("%d ", v);
// 			res.push_back(v);
// 			if(!cal(k, res.size(), res)) res.pop_back();
// 		}

// 		// cout << "kplex_size = " << res.size() << endl;
// 		if(res.size() > kplex.size()) {
// 			kplex = res;
// 		}
// 	}
// 	cout << "\t hec_kplex_size = " << kplex.size() << endl;
// 	// for(auto u : P) {
// 	// 	printf("%d ", u);
// 	// }
// 	// printf("\n");
// }

// //1
// bool Graph::move_u_to_S_with_prune(ui S_end, ui &R_end, ui level) {
//     assert(S_end > 0);
//     ui u = SR[S_end-1];

//     ui neighbors_n = 0, nonneighbors_n = 0;
//     for(ui i = 0;i < R_end;i ++) if(i != S_end-1) {
//         if(s_matrix[u][SR[i]]) neighbors[neighbors_n++] = SR[i];
//         else nonneighbors[nonneighbors_n++] = SR[i];
//     }

//     for(ui i = 0;i < neighbors_n;i ++) ++s_degree_in_S[neighbors[i]];

//     bool is_balanced = true;
//     //
//     if(is_balanced == false) return true;

//     assert(Qv.empty());
//     if(s_degree_in_S[u] + K == S_end) { // only neighbors of u in R can be candidates --- RR2
//         ui i = 0;
//         while(i < nonneighbors_n&&SR_rid[nonneighbors[i]] < S_end) ++ i;
//         for(;i < nonneighbors_n;i ++) { // remove non-neighbors from R
//             assert(level_id[nonneighbors[i]] > level);
//             level_id[nonneighbors[i]] = level;
//             Qv.push(nonneighbors[i]);
//         }
//     }
//     else { // only non-neighbors of u may change their allowance --- RR1
//         ui i = 0;
//         while(i < nonneighbors_n&&SR_rid[nonneighbors[i]] < S_end) ++ i;
//         for(;i < nonneighbors_n;i ++) if(s_degree_in_S[nonneighbors[i]] + K <= S_end) {
//             assert(level_id[nonneighbors[i]] > level);
//             level_id[nonneighbors[i]] = level;
//             Qv.push(nonneighbors[i]);
//         }
//     }

//     // RR2
//     for(ui i = 0;i < nonneighbors_n&&SR_rid[nonneighbors[i]] < S_end;i ++) if(s_degree_in_S[nonneighbors[i]] + K == S_end) {
//         ui v = nonneighbors[i];
//         for(ui j = S_end;j < R_end;j ++) if(level_id[SR[j]] > level&&!s_matrix[v][SR[j]]) {
//             level_id[SR[j]] = level;
//             Qv.push(SR[j]);
//         }
//     }
//     return remove_vertices_and_edges_with_prune(S_end, R_end, level);
// }

// //ok
// bool Graph::remove_vertices_and_edges_with_prune(ui S_end, ui &R_end, ui level) {
//     while(!Qv.empty()) {
//         while(!Qv.empty()) {
//             ui u = Qv.front(); Qv.pop(); // remove u
//             assert(SR[SR_rid[u]] == u);
//             assert(SR_rid[u] >= S_end&&SR_rid[u] < R_end);
//             -- R_end;
//             swap_pos(SR_rid[u], R_end);

//             bool terminate = false;
//             ui neighbors_n = 0;
//             for(ui i = 0;i < R_end;i ++) if(s_matrix[u][SR[i]]) {
//                 ui w = SR[i];
//                 neighbors[neighbors_n++] = w;
//                 -- s_degree[w];
//                 if(s_degree[w] + K <= s_solution_size) {
//                     if(i < S_end) terminate = true; // UB1
//                     else if(level_id[w] > level) { // RR3
//                         level_id[w] = level;
//                         Qv.push(w);
//                     }
//                 }
//             }
//             // UB1
//             if(terminate) {
//                 for(ui i = 0;i < neighbors_n;i ++) ++ s_degree[neighbors[i]];
//                 level_id[u] = s_n;
//                 ++ R_end;
//                 return true;
//             }
//         }
//     }
//     return false;
// }