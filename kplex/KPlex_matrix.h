/*
This file contains code from the Maximum-kPlex project, which is licensed under the MIT License.
The original code and license can be found at: https://github.com/LijunChang/Maximum-kPlex
*/

#ifndef _KPLEX_MATRIX_
#define _KPLEX_MATRIX_

#include "Utility.h"
#include "Timer.h"
#include "LinearHeap.h"

#define _SECOND_ORDER_PRUNING_

class KPLEX_MATRIX {
private:
    long long n;

    int *matrix;
    long long matrix_size;

#ifdef _SECOND_ORDER_PRUNING_
    ui *cn;
    std::queue<std::pair<ui,ui> > Qe;
    std::vector<std::pair<ui, ui> > removed_edges;
    long long removed_edges_n;
#endif

    ui *degree;
    ui *degree_in_S;

    ui K;
    ui *best_solution;
    ui best_solution_size;

    ui *neighbors;
    ui *nonneighbors;

    ui *SR; // union of S and R, where S is at the front
    ui *SR_rid; // reverse ID for SR
    std::queue<ui> Qv;
    ui *level_id;

public:
    KPLEX_MATRIX()
    {
    	n = 0;
        matrix = NULL;
        matrix_size = 0;

        degree = degree_in_S = NULL;
        
        best_solution = NULL;
        K = best_solution_size = 0;

        neighbors = nonneighbors = NULL;

        SR = SR_rid = NULL;
        level_id = NULL;
    }

    ~KPLEX_MATRIX()
    {
        if(matrix != NULL){
            delete[] matrix;
            matrix = NULL;
        }
#ifdef _SECOND_ORDER_PRUNING_
        if(cn != NULL) {
        	delete[] cn;
        	cn = NULL;
        }
#endif
        if(degree != NULL){
            delete[] degree;
            degree = NULL;
        }
        if(degree_in_S != NULL){
            delete[] degree_in_S;
            degree_in_S = NULL;
        }
        if(best_solution != NULL){
            delete[] best_solution;
            best_solution = NULL;
        }
        if(SR != NULL){
            delete[] SR;
            SR = NULL;
        }
        if(SR_rid != NULL){
            delete[] SR_rid;
            SR_rid = NULL;
        }
        if(neighbors != NULL){
            delete[] neighbors;
            neighbors = NULL;
        }
        if(nonneighbors != NULL){
            delete[] nonneighbors;
            nonneighbors = NULL;
        }
        if(level_id != NULL){
        	delete[] level_id;
        	level_id = NULL;
        }
    }

    void allocateMemory(ui n, ui m)
    {
        matrix = new int[m*2];
        matrix_size = m*2;
#ifdef _SECOND_ORDER_PRUNING_
        cn = new ui[m*2];
#endif
        degree = new ui[n];
        degree_in_S = new ui[n];
        best_solution = new ui[n];

        neighbors = new ui[n];
        nonneighbors = new ui[n];

        SR = new ui[n];
        SR_rid = new ui[n];
        level_id = new ui[n];
    }

    // initialize matrix, degree
    void load_graph(ui _n, const std::vector<std::pair<int,int> > &vp, const std::vector<int> sgn) 
    {
        n = _n;
        if(((long long)n)*n > matrix_size) {
        	do {
        		matrix_size *= 2;
        	} while(((long long)n)*n > matrix_size);
        	delete[] matrix; matrix = new int[matrix_size];
#ifdef _SECOND_ORDER_PRUNING_
        	delete[] cn; cn = new ui[matrix_size];
#endif
        }

#ifdef _SECOND_ORDER_PRUNING_
        memset(cn, 0, sizeof(ui)*((long long)n)*n);
#endif
        memset(matrix, 0, sizeof(int)*matrix_size);
        memset(degree, 0, sizeof(ui)*n);

        assert(vp.size() == sgn.size());
        for(ui i = 0; i < vp.size(); i++) {
            assert(vp[i].first >= 0&&vp[i].first < n&&vp[i].second >= 0&&vp[i].second < n);
        	ui a = vp[i].first, b = vp[i].second;
            degree[a]++; degree[b]++;
            matrix[a*n + b] = matrix[b*n + a] = sgn[i];
        }

#ifndef NDEBUG
        printf("load graph of size n=%lld, m=%lu\n", n, vp.size());
        //for(ui i = 0;i < vp.size();i ++) printf("%d %d %d\n", vp[i].first, vp[i].second, sgn[i]);
#endif
    }

    void kPlex(ui K_, std::vector<ui> &kplex, bool must_include_0)
    {
        K = K_;
        best_solution_size = kplex.size();
        ui R_end = 0;
        initialization(R_end, must_include_0);
        if(R_end) kplex_search(0, R_end, 1, must_include_0, n);
        if(best_solution_size > kplex.size()) {
            kplex.clear();
            for(ui i = 0; i < best_solution_size; i++) kplex.push_back(best_solution[i]);
        }
    }

    void heu_kPlex(ui K_, std::vector<ui> &kplex)
    {
        K = K_;
        best_solution_size = 0;
        ui *label = level_id;
        // ui *dOrder = SR;
        memset(label, 0, sizeof(ui)*n);
        memset(degree_in_S, 0, sizeof(ui)*n);

        // ListLinearHeap *heap = new ListLinearHeap(n, n-1);
        // degen(heap, dOrder);

        {
            ui u = n, max_degree = 0;
            for(int i = 0; i < n; i++) {
                assert(label[i] == 0);
                if(degree[i] > max_degree) {
                    max_degree = degree[i];
                    u = i;
                }
            }
            printf("u = %d, maxdegree = %d\n", u, max_degree);
            assert(u != n);
            label[u] = 3;
            best_solution[best_solution_size++] = u;
            int *t_matrix = matrix + u*n;

            ui neighbors_n = 0, nonneighbors_n = 0;
            for(ui i = 0; i < n; i++) if(i != u) {
                if(t_matrix[i]) neighbors[neighbors_n++] = i;
                else nonneighbors[nonneighbors_n++] = i;
            }

            for(ui i = 0; i < neighbors_n; i++) degree_in_S[neighbors[i]]++;

            for(ui i = 0; i < neighbors_n; i++) {
                ui v = neighbors[i];
                if(t_matrix[v] == 1) {
                    label[v] = 1;
                }
                else if(t_matrix[v] == -1) {
                    label[v] = 2;
                }
            }
            int c[10];
            for(int i = 0; i < 6; i++) c[i] = 0;
            for(ui i = 0; i < n; i++) {
                assert(label[i] >= 0 && label[i] < 6);
                c[label[i]]++;
            }
            for(int i = 0; i < 6; i++) {
                printf("%d ", c[i]);
            }
            printf("\n");
        }

        while(1) {
            // u is the next vertex
            ui u = n, max_degree = 0, inP = 0;
            
            for(ui i = 0; i < n; i++) {
                assert(label[i] >= 0 && label[i] < 6);
                // if(label[i] > 3) {
                //     u = i;
                //     inP = label[i];
                //     break;
                // }



                if(label[i] == 1 && degree[i] > max_degree) {
                    max_degree = degree[i];
                    inP = 1;
                    u = i;
                }
                else if(label[i] == 2 && degree[i] > max_degree) {
                    max_degree = degree[i];
                    inP = 2;
                    u = i;
                }
            }
            printf("u = %d, maxdegree = %d\n", u, max_degree);
            if(u == n) break;
            best_solution[best_solution_size++] = u;
            int *t_matrix = matrix + u*n;
            assert(label[u] == 1 || label[u] == 2);
            assert((inP == 1 || inP == 2));

            ui neighbors_n = 0, nonneighbors_n = 0;
            for(ui i = 0; i < n; i++) if(i != u) {
                if(t_matrix[i]) neighbors[neighbors_n++] = i;
                else nonneighbors[nonneighbors_n++] = i;
            }

            for(ui i = 0; i < neighbors_n; i++) degree_in_S[neighbors[i]]++;

            // after adding u, it's necessary to check the non neighbors of u and u

            // check if the nonneighbors of u in R can be candidates
            if(degree_in_S[u] + K == best_solution_size) {
                for(ui i = 0; i < nonneighbors_n; i++) {
                    int v = nonneighbors[i];
                    if(label[v] < 3) label[v] = 5;
                }
            }
            else {
                for(ui i = 0; i < nonneighbors_n; i++) {
                    int v = nonneighbors[i];
                    if(label[v] < 3 && degree_in_S[v] + K <= best_solution_size) label[v] = 5;
                }
            }

            // check the neighbors of nodes in S
            for(ui i = 0; i < nonneighbors_n; i++) {
                ui v = nonneighbors[i];
                if(!(label[v] == 3 || label[v] == 4)) continue;
                if(degree_in_S[v] + K == best_solution_size) {
                    int *tt_matrix = matrix + v*n;
                    for(ui j = 0; j < n; j++) if(v != j) {
                        if(label[j] < 3 && !tt_matrix[j]) label[j] = 5;
                    }
                }
            }

            {
                int c[10];
                for(int i = 0; i < 6; i++) c[i] = 0;
                for(ui i = 0; i < n; i++) {
                    assert(label[i] >= 0 && label[i] < 6);
                    c[label[i]]++;
                }
                for(int i = 0; i < 6; i++) {
                    printf("%d ", c[i]);
                }
                printf("#\n");
            }


            if(inP == 1) {
                label[u] = 3;
                for(ui i = 0; i < n; i++) if(u != i) {
                    assert(label[i] >= 0 && label[i] < 6);
                    if(t_matrix[i] == 1) {
                        if(label[i] == 0) label[i] = 1;
                        if(label[i] == 2) label[i] = 5;
                    }
                    else if(t_matrix[i] == -1) {
                        if(label[i] == 0) label[i] = 2;
                        if(label[i] == 1) label[i] = 5;
                    }
                }
            }
            else if(inP == 2) {
                label[u] = 4;
                for(ui i = 0; i < n; i++) if(u != i) {
                    assert(label[i] >= 0 && label[i] < 6);
                    if(t_matrix[i] == 1) {
                        if(label[i] == 0) label[i] = 2;
                        if(label[i] == 1) label[i] = 5;
                    }
                    else if(t_matrix[i] == -1) {
                        if(label[i] == 0) label[i] = 1;
                        if(label[i] == 2) label[i] = 5;
                    }
                }
            }

            // {
            //     int c[10];
            //     for(int i = 0; i < 6; i++) c[i] = 0;
            //     for(ui i = 0; i < n; i++) {
            //         assert(label[i] >= 0 && label[i] < 6);
            //         c[label[i]]++;
            //     }
            //     for(int i = 0; i < 6; i++) {
            //         printf("%d ", c[i]);
            //     }
            //     printf("\n");
            // }
        }
        if(best_solution_size > kplex.size()) {
            kplex.clear();
            for(ui i = 0; i < best_solution_size; i++) kplex.push_back(best_solution[i]);
        }
    }

private:
    // initialize degree_in_S, R_end, level_id
    void initialization(ui &R_end, bool must_include_0) {
        memset(degree_in_S, 0, sizeof(ui)*n);
        R_end = 0;
        for(ui i = 0; i < n; i++) SR_rid[i] = n;
        for(ui i = 0; i < n; i++) if(degree[i] + K > best_solution_size){
            SR[R_end] = i; SR_rid[i] = R_end;
            R_end++;
        }

        for(ui i = 0; i < n; i++) level_id[i] = n;

        assert(Qv.empty());

#ifdef _SECOND_ORDER_PRUNING_
        for(ui i = 0;i < R_end;i ++) {
        	ui neighbors_n = 0;
        	int *t_matrix = matrix + SR[i]*n;
        	for(ui j = 0;j < R_end;j ++) if(t_matrix[SR[j]]) neighbors[neighbors_n ++] = SR[j];
        	for(ui j = 0;j < neighbors_n;j ++) for(ui k = j+1;k < neighbors_n;k ++) {
        		++ cn[neighbors[j]*n + neighbors[k]];
        		++ cn[neighbors[k]*n + neighbors[j]];
        	}
        }

        while(!Qe.empty()) Qe.pop();
        for(ui i = 0;i < R_end;i ++) for(ui j = i+1;j < R_end;j ++) {
        	if(matrix[SR[i]*n + SR[j]]&&upper_bound_based_prune(0, SR[i], SR[j])) {
        		Qe.push(std::make_pair(SR[i], SR[j]));
        	}
        }
        removed_edges_n = 0;
#endif
    }

    void kplex_search(ui S_end, ui R_end, ui level, bool choose_zero, ui last_choice) {
        if(S_end > best_solution_size) { // find a larger solution
            best_solution_size = S_end;
            for(ui i = 0; i < best_solution_size; i++) best_solution[i] = SR[i];
        }
        if(R_end <= best_solution_size) return;

        if(!calc_upper_bound_partition(S_end, R_end)) return;

        // choose branching vertex
        bool must_include = false;
        ui u = n; // u is the branching vertex
        if(choose_zero) {
        	assert(S_end == 0);
        	if(SR_rid[0] >= R_end) return;
        	u = 0;
        	must_include = true;
        }
        else {
        	u = choose_branch_vertex(S_end, R_end);
            assert(check_balance(S_end, u));
        }
        if(u == n) return;
        // assert(u != n);
        assert(degree[u] + K > best_solution_size&&degree[u] + K > S_end);

        // the first branch includes u into S
        bool pruned = true;
        ui pre_best_solution_size = best_solution_size, old_R_end = R_end;
        ui old_removed_edges_n = 0;
#ifdef  _SECOND_ORDER_PRUNING_
        old_removed_edges_n = removed_edges_n;
#endif

        assert(SR[SR_rid[u]] == u&&SR_rid[u] >= S_end&&SR_rid[u] < R_end);
        swap_pos(S_end, SR_rid[u]);
        S_end++;
        pruned = move_u_to_S(S_end, R_end, level);
        if(!pruned) kplex_search(S_end, R_end, level+1, false, u);
        restore_SR(S_end, R_end, old_R_end, level, old_removed_edges_n);
#ifdef  _SECOND_ORDER_PRUNING_
        assert(removed_edges_n == old_removed_edges_n);
#endif

        if(must_include) {
        	move_u_to_R(S_end, R_end, level);
        	return;
        }

        // the second branch exclude u from S
        assert(Qv.empty());
#ifdef _SECOND_ORDER_PRUNING_
        while(!Qe.empty()) Qe.pop();
#endif
        pruned = remove_u_from_SR(S_end, R_end, level);
        if(!pruned && best_solution_size > pre_best_solution_size) pruned = collect_removable_vertices(S_end, R_end, level);
        if(!pruned) {
            if(!remove_vertices_from_R(S_end, R_end, level)) {
            	kplex_search(S_end, R_end, level+1, false, last_choice);
            }
        }
        restore_SR(S_end, R_end, old_R_end, level, old_removed_edges_n);
    }

    bool move_u_to_S(ui S_end, ui &R_end, ui level)
    {
    	assert(S_end > 0);
        ui u = SR[S_end-1];
        int *t_matrix = matrix + u*n;

        ui neighbors_n = 0, nonneighbors_n = 0;
        for(ui i = 0; i < R_end; i++) if(i != S_end-1) {
            if(t_matrix[SR[i]]) neighbors[neighbors_n++] = SR[i];
            else nonneighbors[nonneighbors_n++] = SR[i];
        }

        for(ui i = 0; i < neighbors_n; i++) degree_in_S[neighbors[i]]++;

        // after adding u, it's necessary to check the non-neighbors of u and u
        assert(Qv.empty());
        // check if the nonneighbors of u in R can be candidates
        if(degree_in_S[u] + K == S_end) {
        	ui i = 0;
        	while(i < nonneighbors_n && SR_rid[nonneighbors[i]] < S_end) i++;
            for(; i < nonneighbors_n; i++) {
            	assert(level_id[nonneighbors[i]] > level);
            	level_id[nonneighbors[i]] = level;
                Qv.push(nonneighbors[i]);
            }
        }
        else {
        	ui i = 0;
        	while(i < nonneighbors_n && SR_rid[nonneighbors[i]] < S_end) i++;
            for(; i < nonneighbors_n; i++) if(degree_in_S[nonneighbors[i]] + K <= S_end) {
            	assert(level_id[nonneighbors[i]] > level);
            	level_id[nonneighbors[i]] = level;
                Qv.push(nonneighbors[i]);
            }
        }

        // check the neighbors of nodes in S
        for(ui i = 0; i < nonneighbors_n && SR_rid[nonneighbors[i]] < S_end; i++) if(degree_in_S[nonneighbors[i]] + K == S_end) {
            int *tt_matrix = matrix + nonneighbors[i]*n;
            for(ui j = S_end; j < R_end; j++) if(level_id[SR[j]] > level && !tt_matrix[SR[j]]) {
            	level_id[SR[j]] = level;
                Qv.push(SR[j]);
            }
        }

#ifdef _SECOND_ORDER_PRUNING_
        // RR4
        for(ui i = 0;i < nonneighbors_n;i ++) {
            int v = nonneighbors[i];
            assert(!t_matrix[v]);
            if(SR_rid[v] < S_end||level_id[v] == level||t_matrix[v]) continue;
            if(upper_bound_based_prune(S_end, u, v)) {
                level_id[v] = level;
                Qv.push(v);
            }
        }

        // update cn(.,.)
        for(ui i = 0;i < neighbors_n;i ++) { // process common neighbors of u
            for(ui j = i+1;j < neighbors_n;j ++) {
// #ifndef NDEBUG
//                 if(!cn[neighbors[i]*n + neighbors[j]]) {
//                     printf("cn[neighbors[i]*n + neighbors[j]]: %u %u\n", cn[neighbors[i]*n + neighbors[j]], cn[neighbors[j]*n + neighbors[i]]);
//                 }
// #endif
                assert(cn[neighbors[i]*n + neighbors[j]]);
                -- cn[neighbors[i]*n + neighbors[j]];
                -- cn[neighbors[j]*n + neighbors[i]];
            }
        }

        while(!Qe.empty()) Qe.pop();
        int new_n = 0;
        for(ui i = 0;i < nonneighbors_n;i ++) if(level_id[nonneighbors[i]] > level) nonneighbors[new_n ++] = nonneighbors[i];
        nonneighbors_n = new_n;
        for(ui i = 1;i < nonneighbors_n;i ++) { // process common non-neighbors of u
            ui w = nonneighbors[i];
            for(ui j = 0;j < i;j ++) {
                ui v = nonneighbors[j];
                if(!upper_bound_based_prune(S_end, v, w)) continue;
                if(SR_rid[w] < S_end) return true; // v, w \in S --- UB2
                else if(SR_rid[v] >= S_end) { // v, w, \in R --- RR5
                    if(matrix[v*n + w]) Qe.push(std::make_pair(v,w));
                }
                else { // RR4
                    assert(level_id[w] > level);
                    level_id[w] = level;
                    Qv.push(w);
                    break;
                }
            }
        }
#endif

        // // check balance
        // {
        //     ui i = 0;
        //     while(i < neighbors_n && SR_rid[neighbors[i]] < S_end) i++;
        //     for(ui j = 0; j < i; j++) {
        //         ui v = neighbors[j];
        //         for(ui k = i; k < neighbors_n; k++) {
        //             ui w = neighbors[k];
        //             if(!matrix[v*n + w]) continue;
        //             if(level_id[w] > level) {
        //                 ui tri = matrix[v*n + w] + t_matrix[v] + t_matrix[w];
        //                 assert(tri == 3 || tri == 1 || tri == -1 || tri == -3);

        //                 if(tri == 1 || tri == -3) {
        //                     level_id[w] = level;
        //                     Qv.push(w);
        //                 }
        //             }
        //         }
        //     }
        // }
        // check balance
        for(ui i = S_end; i < R_end; i++) {
            if(level_id[SR[i]] > level && !check_balance(S_end, SR[i])) {
            	level_id[SR[i]] = level;
                Qv.push(SR[i]);
            }
        }
        return remove_vertices_from_R(S_end, R_end, level);
    }

    bool remove_vertices_from_R(ui S_end, ui &R_end, ui level)
    {
#ifdef _SECOND_ORDER_PRUNING_
        while(!Qv.empty()||!Qe.empty()) {
#else
        while(!Qv.empty()) {
#endif
            while(!Qv.empty()) {
                ui u = Qv.front(); Qv.pop(); // remove u
                assert(SR[SR_rid[u]] == u);
                assert(SR_rid[u] >= S_end&&SR_rid[u] < R_end);
                R_end--;
                swap_pos(SR_rid[u], R_end);

                bool terminate = false;
                ui neighbors_n = 0;
                int *t_matrix = matrix + u*n;
                for(ui i = 0; i < R_end; i++) if(t_matrix[SR[i]]) {
                    ui w = SR[i];
                    neighbors[neighbors_n++] = w;
                    degree[w]--;
                    if(degree[w] + K <= best_solution_size) {
                        if(i < S_end) terminate = true; // UB1
                        else if(level_id[w] > level) { // RR3
                            level_id[w] = level;
                            Qv.push(w);
                        }
                    }
                }

                if(terminate) {
                    for(ui i = 0; i < neighbors_n; i++) degree[neighbors[i]]++;
                    level_id[u] = n;
                    R_end++;
                    return true;
                }

#ifdef _SECOND_ORDER_PRUNING_
                for(ui i = 1;i < neighbors_n;i ++) {
                    ui w = neighbors[i];
                    for(ui j = 0;j < i;j ++) {
                        ui v = neighbors[j];
                        assert(cn[v*n+w]);
#ifndef NDEBUG
                        ui common_neighbors = 0;
                        for(ui k = S_end;k <= R_end;k ++) if(matrix[SR[k]*n + v]&&matrix[SR[k]*n + w]) ++ common_neighbors;
                        assert(cn[v*n + w] == common_neighbors);
                        assert(cn[w*n + v] == common_neighbors);
#endif
                        -- cn[v*n + w];
                        -- cn[w*n + v];
#ifndef NDEBUG
                        common_neighbors = 0;
                        for(ui k = S_end;k < R_end;k ++) if(matrix[SR[k]*n + v]&&matrix[SR[k]*n + w]) ++ common_neighbors;
                        assert(cn[v*n + w] == common_neighbors);
                        assert(cn[w*n + v] == common_neighbors);
#endif

                        if(!upper_bound_based_prune(S_end, v, w)) continue;

                        if(SR_rid[w] < S_end) terminate = true; // v, w \in S --- UB2
                        else if(SR_rid[v] >= S_end) { // v, w, \in R --- RR5
                            if(matrix[v*n + w]) Qe.push(std::make_pair(v,w));
                        }
                        else if(level_id[w] > level) { // RR4
                            level_id[w] = level;
                            Qv.push(w);
                        }
                    }
                }
                if(terminate) return true;
#endif
            }

#ifdef _SECOND_ORDER_PRUNING_
        	if(Qe.empty()) break;

        	ui v = Qe.front().first, w =  Qe.front().second; Qe.pop();
        	if(level_id[v] <= level||level_id[w] <= level||!matrix[v*n + w]) continue;
        	assert(SR_rid[v] >= S_end&&SR_rid[v] < R_end&&SR_rid[w] >= S_end&&SR_rid[w] < R_end);

        	if(degree[v] + K <= best_solution_size + 1) {
        		level_id[v] = level;
        		Qv.push(v);
        	}
        	if(degree[w] + K <= best_solution_size + 1) {
        		level_id[w] = level;
        		Qv.push(w);
        	}
        	if(!Qv.empty()) continue;

#ifndef NDEBUG
        	//printf("remove edge between %u and %u\n", v, w);
#endif

        	assert(matrix[v*n + w]);
        	matrix[v*n + w] = matrix[w*n + v] = 0;
        	-- degree[v]; -- degree[w];

        	if(removed_edges.size() == removed_edges_n) {
        		removed_edges.push_back(std::make_pair(v,w));
        		++ removed_edges_n;
        	}
        	else removed_edges[removed_edges_n ++] = std::make_pair(v,w);

        	int *t_matrix = matrix + v*n;
        	for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) {
        		-- cn[w*n + SR[i]];
        		-- cn[SR[i]*n + w];
        		if(!upper_bound_based_prune(S_end, w, SR[i])) continue;
        		if(i < S_end) {
        			if(level_id[w] > level) {
        				level_id[w] = level;
        				Qv.push(w);
        			}
        		}
        		else if(matrix[w*n + SR[i]]) Qe.push(std::make_pair(w, SR[i]));
        	}
        	t_matrix = matrix + w*n;
        	for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) {
        		-- cn[v*n + SR[i]];
        		-- cn[SR[i]*n + v];
        		if(!upper_bound_based_prune(S_end, v, SR[i])) continue;
        		if(i < S_end) {
        			if(level_id[v] > level) {
        				level_id[v] = level;
        				Qv.push(v);
        			}
        		}
        		else if(matrix[v*n + SR[i]]) Qe.push(std::make_pair(v, SR[i]));
        	}
#endif
        }

        return false;
    }

    void restore_SR(ui S_end, ui &R_end, ui old_R_end, ui level, ui old_removed_edges_n)
    {
        while(!Qv.empty()) {
            ui u = Qv.front(); Qv.pop();
            assert(level_id[u] == level&&SR_rid[u] < R_end);
            level_id[u] = n;
        }
        while(R_end < old_R_end) { // insert u back into R
            ui u = SR[R_end];
            assert(level_id[u] == level&&SR_rid[u] == R_end);
            level_id[u] = n;

            int *t_matrix = matrix + u*n;
            for(ui i = 0; i < R_end; i++) if(t_matrix[SR[i]]) {
                degree[SR[i]]++;
            }

#ifdef _SECOND_ORDER_PRUNING_
            ui neighbors_n = 0;
            for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) {
                neighbors[neighbors_n++] = SR[i];
            }
            for(ui i = 0;i < neighbors_n;i ++) {
                ui v = neighbors[i];
                for(ui j = i + 1;j < neighbors_n;j ++) {
                    ui w = neighbors[j];
                    ++ cn[v*n + w];
                    ++ cn[w*n + v];
                }
            }
            ui *t_cn = cn + u*n;
            for(ui i = 0;i < R_end;i ++) t_cn[SR[i]] = 0;
            for(ui i = 0;i < neighbors_n;i ++) if(SR_rid[neighbors[i]] >= S_end) {
            	ui v = neighbors[i];
            	int *t_matrix = matrix + v*n;
            	for(ui j = 0;j < R_end;j ++) if(t_matrix[SR[j]]) ++ t_cn[SR[j]];
            }
            for(ui i = 0;i < R_end;i ++) {
            	cn[SR[i]*n + u] = t_cn[SR[i]];
#ifndef NDEBUG
            	ui common_neighbors = 0, v = SR[i], w = u;
            	for(ui k = S_end;k < R_end;k ++) if(matrix[SR[k]*n + v]&&matrix[SR[k]*n + w]) ++ common_neighbors;
            	if(t_cn[SR[i]] != common_neighbors) printf("t_cn[SR[i]] = %u, comon_neighbors = %u\n", t_cn[SR[i]], common_neighbors);
            	assert(t_cn[SR[i]] == common_neighbors);
#endif
            }
#endif

            R_end++;
        }
    }
    
    void move_u_to_R(ui &S_end, ui &R_end, ui level)
    {
    	assert(S_end);
        ui u = SR[--S_end];
        ui neighbors_n = 0;
        int *t_matrix = matrix + u*n;
        for(ui i = 0; i < R_end; i++) if(t_matrix[SR[i]]) neighbors[neighbors_n++] = SR[i];
        for(ui i = 0; i < neighbors_n; i++) degree_in_S[neighbors[i]]--;

#ifdef _SECOND_ORDER_PRUNING_
        for(ui i = 0;i < neighbors_n;i ++) {
        	ui v = neighbors[i];
        	for(ui j = i+1;j < neighbors_n;j ++) {
        		ui w = neighbors[j];
        		++ cn[v*n + w];
        		++ cn[w*n + v];
        	}
        }
#endif
    }

    bool remove_u_from_SR(ui &S_end, ui &R_end, ui level)
    {
    	assert(S_end);
		ui u = SR[S_end-1];
		S_end--; R_end--;
		swap_pos(S_end, R_end);
		level_id[u] = level;

		bool terminate = false;
        ui neighbors_n = 0;
        int *t_matrix = matrix + u*n;
        for(ui i = 0; i < R_end; i++) if(t_matrix[SR[i]]) {
            ui v = SR[i];
        	degree_in_S[v]--;
        	degree[v]--;
        	if(degree[v] + K <= best_solution_size) {
        		if(i < S_end) terminate = true;
        		else {
        			assert(level_id[v] > level);
        			level_id[v] = level;
        			Qv.push(v);
        		}
        	}
        }
        if(terminate) return true;

#ifdef _SECOND_ORDER_PRUNING_
        for(ui i = 1;i < neighbors_n;i ++) if(level_id[neighbors[i]] > level) {
        	ui w = neighbors[i];
        	for(ui j = 0;j < i;j ++) {
        		ui v = neighbors[j];
        		if(!upper_bound_based_prune(S_end, v, w)) continue;

        		if(SR_rid[w] < S_end) return true; // v, w \in S
				else if(SR_rid[v] >= S_end) { // v, w, \in R
					if(matrix[v*n + w]) Qe.push(std::make_pair(v,w));
				}
				else {
					assert(level_id[w] > level);
					level_id[w] = level;
					Qv.push(w);
					break;
				}
        	}
        }
#endif
		return false;
	}

    bool collect_removable_vertices(ui S_end, ui R_end, ui level)
    {
    	for(ui i = 0;i < S_end;i ++) if(degree[SR[i]] + K <= best_solution_size) return true;

#ifdef _SECOND_ORDER_PRUNING_
    	for(ui i = 0;i < S_end;i ++) for(ui j = i+1;j < S_end;j ++) {
    		if(upper_bound_based_prune(S_end, SR[i], SR[j])) return true;
    	}
#endif

    	for(ui i = S_end; i < R_end; i++) if(level_id[SR[i]] > level) {
            ui v = SR[i];
    		if(degree_in_S[v] + K <= S_end || degree[v] + K <= best_solution_size) { //RR1, RR3
    			level_id[v] = level;
    			Qv.push(v);
    			continue;
    		}
    		int *t_matrix = matrix + v*n;
    		for(ui j = 0; j < S_end; j++) {
                ui w = SR[j];
#ifdef _SECOND_ORDER_PRUNING_
    			if((S_end - degree_in_S[SR[j]] == K&&!t_matrix[SR[j]])||upper_bound_based_prune(S_end, SR[i], SR[j])) {
#else
    			if(degree_in_S[w] + K == S_end && !t_matrix[w]) { //RR2
#endif
    				level_id[v] = level;
    				Qv.push(v);
    				break;
    			}
    		}
    	}

#ifdef _SECOND_ORDER_PRUNING_
    	for(ui i = S_end;i < R_end;i ++) if(level_id[SR[i]] > level) {
    		for(ui j = i+1;j < R_end;j ++) if(level_id[SR[j]] > level&&matrix[SR[i]*n + SR[j]]) {
    			if(upper_bound_based_prune(S_end, SR[i], SR[j])) Qe.push(std::make_pair(SR[i], SR[j]));
    		}
    	}
#endif

        return false;
    }

#ifdef _SECOND_ORDER_PRUNING_
    bool upper_bound_based_prune(ui S_end, ui u, ui v) {
    	// ui ub = S_end + 2*K - (S_end - degree_in_S[u]) - (S_end - degree_in_S[v]) + cn[u*n + v];
    	ui ub = 2*K + degree_in_S[u] - S_end + degree_in_S[v] + cn[u*n + v];
    	if(SR_rid[u] >= S_end) {
    		-- ub; // S_end ++
    		if(matrix[u*n+v]) ++ ub; // degree_in_S[v] ++
    	}
    	if(SR_rid[v] >= S_end) {
    		-- ub;
    		if(matrix[v*n+u]) ++ ub;
    	}
    	return ub <= best_solution_size;
    }
#endif

    void swap_pos(ui i, ui j)
    {
        std::swap(SR[i], SR[j]);
        SR_rid[SR[i]] = i;
        SR_rid[SR[j]] = j;
    }

    ui choose_branch_vertex(ui S_end, ui R_end)
    {
        ui u = n, min_degree_in_S = n;
        for(ui i = S_end;i < R_end;i ++) {
            ui v = SR[i];
            // if(check_balance(S_end, v)) {
                u = v;
                break;
            // }
        }
        // assert(u != n);
        return u;
    }
    // ui choose_branch_vertex(ui S_end, ui R_end)
    // {
    // 	ui u = n, min_degree_in_S = n;
    // 	for(ui i = 0;i < R_end;i ++) {
    // 		ui v = SR[i];
    // 		if(degree_in_S[v] < min_degree_in_S) {
    // 			u = v;
    // 			min_degree_in_S = degree_in_S[v];
    // 		}
    // 		else if(degree_in_S[v] == min_degree_in_S) {
    // 			if((i < S_end&&degree[v] < degree[u])||(i >= S_end&&degree[v] > degree[u])) {
    // 				u = v;
    // 				min_degree_in_S = degree_in_S[v];
    // 			}
    // 		}
    // 	}
    // 	assert(u != n);
    // 	if(SR_rid[u] < S_end) {
    // 		int *t_matrix = matrix+u*n;
    // 		u = n;
    // 		ui max_degree = 0;
    // 		for(ui i = S_end;i < R_end;i ++) if(!t_matrix[SR[i]]&&degree[SR[i]] > max_degree) {
    // 			max_degree = degree[SR[i]];
    // 			u = SR[i];
    // 		}
    // 	}
    //     if(u == n) u = SR[S_end];
    // 	assert(u != n);
    // 	return u;
    // }

    bool check_balance(ui S_end, ui u)
    {
        ui neighbors_n = 0;
        int *t_matrix = matrix + u*n;
        for(ui i = 0; i < S_end; i++) if(t_matrix[SR[i]]) neighbors[neighbors_n++] = SR[i];
        for(ui i = 0; i < neighbors_n; i++) {
            for(ui j = i+1; j < neighbors_n; j++) {
                ui v = neighbors[i], w = neighbors[j];
                if(!matrix[v*n + w]) continue;
                ui tri = matrix[v*n + w] + t_matrix[v] + t_matrix[w];
                assert(tri == 3 || tri == 1 || tri == -1 || tri == -3);
                if(tri == 1 || tri == -3) return false;
            }
        }
        return true;
    }

    // degeneracy-based k-plex
    // return an upper bound of the maximum k-plex size
    // return dOrder
    // void degen(ListLinearHeap *heap, ui *dOrder)
    // {
    // 	int *peel_sequence = new ui[n];
    // 	int *vis = new ui[n];
        // 
    // 	for(ui i = 0; i < n; i++) peel_sequence[i] = i;
    // 	for(ui i = 0; i < n; i++) vis[i] = 0;
// 
    // 	heap->init(n, n-1, peel_sequence, degree);
    // 	for(ui i = 0; i < n; i++) {
    // 		int u, deg;
    // 		heap->pop_min(u, deg);
    // 		dOrder[n - i] = u;
// 
    //         int *t_matrix = matrix + u*n;
    //         for(ui j = 0; j < n; j++) if(j != u && t_matrix[j]) {
    // 			if(vis[j] == 0) heap->decrement(j, 1);
    // 		}
    // 		vis[u] = 1;
    // 	}
// 
    // 	// for(ui i = 0; i < n; i++) {
    // 	// 	printf("%d ", dOrder[i]);
    // 	// }
    // 	// printf("\n");
// 
    // 	delete [] peel_sequence;
    // 	delete [] vis;
    // }

    bool calc_upper_bound_partition(ui S_end, ui R_end) {
        ui ub = 0, pi0 = R_end - S_end;
        ui *label = neighbors;
        ui *missing_edges = nonneighbors;
        memset(label, 0, sizeof(ui)*n);

        //calculate missing_edges
        for(ui i = 0; i < S_end; i++) {
            missing_edges[i] = S_end - degree_in_S[SR[i]] - 1;
        }

        for(ui i = 0; i < S_end; i++) {
            ui u = SR[i];
            ui cn = 0;
            int *t_matrix = matrix + u*n;
            for(ui j = S_end; j < R_end; j++) if(!label[SR[j]] && !t_matrix[SR[j]]) {
                label[SR[j]] = 1;
                cn++;
                pi0--;
            }
            ub = ub + min(K - 1 - missing_edges[i], cn);
        }

        ub = S_end + pi0 + ub;

        // printf("%d %d\n", ub, best_solution_size);

        if(ub <= best_solution_size) {
            printf("!\n");
            return false;
        }
        return true;
    }
};

#endif