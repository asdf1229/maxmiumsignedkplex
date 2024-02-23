/*
This file contains code from the Maximum-kPlex project, which is licensed under the MIT License.
The original code and license can be found at: https://github.com/LijunChang/Maximum-kPlex
*/

#ifndef _KPLEX_MATRIX_
#define _KPLEX_MATRIX_

#include "Utility.h"
#include "Timer.h"
#include "LinearHeap.h"

class KPLEX_MATRIX {
private:
    long long n;

    int *matrix;
    long long matrix_size;

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
        }

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
        if(R_end) _kplex_search(0, R_end, 1);
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
    }

    void _kplex_search(ui S_end, ui R_end, ui level) {
        //check
        

        if(S_end > best_solution_size) { // find a larger solution
            best_solution_size = S_end;
            for(ui i = 0; i < best_solution_size; i++) best_solution[i] = SR[i];
        }
        if(R_end <= best_solution_size) return;

        // choose branching vertex
        ui u = SR[S_end];

        // the first branch includes u into S
        bool pruned = true;
        ui pre_best_solution_size = best_solution_size, old_R_end = R_end;

        assert(SR[SR_rid[u]] == u);
        assert(SR[SR_rid[u]] == u&&SR_rid[u] >= S_end&&SR_rid[u] < R_end);
        
        swap_pos(S_end, SR_rid[u]); S_end++;
        move_u_to_S(S_end, R_end, level);
        if(check_balance(S_end-1, u) && check_kplex(S_end)) {
            _kplex_search(S_end, R_end, level+1);
        }
        restore_SR(S_end, R_end, old_R_end, level);

        // the second branch exclude u from S
        remove_u_from_SR(S_end, R_end, level);
        _kplex_search(S_end, R_end, level+1);
        restore_SR(S_end, R_end, old_R_end, level);
    }

    void kplex_search(ui S_end, ui R_end, ui level, bool choose_zero) {
        if(S_end > best_solution_size) { // find a larger solution
            best_solution_size = S_end;
            for(ui i = 0; i < best_solution_size; i++) best_solution[i] = SR[i];
        }
        if(R_end <= best_solution_size) return;

        //choose branching vertex
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
            // assert(check_balance(S_end, u));
        }
        if(u == n) return;
        // assert(u != n);
        assert(degree[u] + K > best_solution_size&&degree[u] + K > S_end);

        // the first branch includes u into S
        bool pruned = true;
        ui pre_best_solution_size = best_solution_size, old_R_end = R_end;

        assert(SR[SR_rid[u]] == u&&SR_rid[u] >= S_end&&SR_rid[u] < R_end);
        swap_pos(S_end, SR_rid[u]);
        S_end++;
        pruned = move_u_to_S(S_end, R_end, level);
        if(!pruned) kplex_search(S_end, R_end, level+1, false);
        restore_SR(S_end, R_end, old_R_end, level);

        if(must_include) {
        	move_u_to_R(S_end, R_end, level);
        	return;
        }

        // the second branch exclude u from S
        assert(Qv.empty());
        pruned = remove_u_from_SR(S_end, R_end, level);
        if(!pruned && best_solution_size > pre_best_solution_size) pruned = collect_removable_vertices(S_end, R_end, level);
        if(!pruned) {
            if(!remove_vertices_from_R(S_end, R_end, level)) {
            	kplex_search(S_end, R_end, level+1, false);
            }
        }
        restore_SR(S_end, R_end, old_R_end, level);
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

        return true;
    }

    bool remove_vertices_from_R(ui S_end, ui &R_end, ui level)
    {
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
		}
        return false;
    }

    void restore_SR(ui S_end, ui &R_end, ui old_R_end, ui level)
    {
        while(R_end < old_R_end) { // insert u back into R
            ui u = SR[R_end];
            assert(level_id[u] == level&&SR_rid[u] == R_end);
            level_id[u] = n;
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
    }

    bool remove_u_from_SR(ui &S_end, ui &R_end, ui level)
    {
    	assert(S_end);
		ui u = SR[S_end-1];
		S_end--; R_end--;
		swap_pos(S_end, R_end);
		level_id[u] = level;

        int *t_matrix = matrix + u*n;
        for(ui i = 0; i < R_end; i++) if(t_matrix[SR[i]]) {
            ui v = SR[i];
        	degree_in_S[v]--;
        }
		return false;
	}

    bool collect_removable_vertices(ui S_end, ui R_end, ui level)
    {
    	for(ui i = 0;i < S_end;i ++) if(degree[SR[i]] + K <= best_solution_size) return true;

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
    			if(degree_in_S[w] + K == S_end && !t_matrix[w]) { //RR2
    				level_id[v] = level;
    				Qv.push(v);
    				break;
    			}
    		}
    	}
        return false;
    }

    void swap_pos(ui i, ui j)
    {
        std::swap(SR[i], SR[j]);
        SR_rid[SR[i]] = i;
        SR_rid[SR[j]] = j;
    }

    ui choose_branch_vertex(ui S_end, ui R_end)
    {
        ui u = n;
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

    bool check_kplex(ui S_end)
    {
        for(ui i = 0; i < S_end; i++) {
            if(degree_in_S[SR[i]] + K < S_end) {
                return false;
            }
        }
        return true;
    }
};

#endif