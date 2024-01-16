#ifndef _KPLEX_BB_MATRIX_
#define _KPLEX_BB_MATRIX_

#include "Utility.h"
#include "Timer.h"

class KPLEX_BB_MATRIX {
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
    KPLEX_BB_MATRIX()
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

    ~KPLEX_BB_MATRIX()
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

    //?> initialize matrix, degree
    void load_graph(ui _n, const std::vector<std::pair<int,int> > &vp) 
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
        for(ui i = 0; i < vp.size(); i++) {
            assert(vp[i].first >= 0&&vp[i].first < n&&vp[i].second >= 0&&vp[i].second < n);
        	ui a = vp[i].first, b = vp[i].second;
            degree[a]++; degree[b]++;
            matrix[a*n + b] = matrix[b*n + a] = 1;
        }

#ifndef NDEBUG
        printf("load graph of size n=%lld, m=%lu\n", n, vp.size());
        //for(ui i = 0;i < vp.size();i ++) printf("%d %d\n", vp[i].first, vp[i].second);
#endif
    }

    void kPlex(ui K_, std::vector<ui> &kplex, bool must_include_0)
    {
        K = K_;
        best_solution_size = kplex.size();
        ui R_end = 0;
        initialization(R_end, must_include_0);
        if(R_end) kplex_search(0, R_end, 1, must_include_0);
        if(best_solution_size > kplex.size()) {
            kplex.clear();
            for(int i = 0; i < best_solution_size; i++) kplex.push_back(best_solution[i]);
        }
    }

private:
    // initialize degree_in_S, level_id
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

    void kplex_search(ui S_end, ui R_end, ui level, bool choose_zero) {
        if(S_end > best_solution_size) { // find a larger solution
            best_solution_size = S_end;
            for(ui i = 0; i < best_solution_size; i++) best_solution[i] = SR[i];
        }
        if(R_end <= best_solution_size) return ;

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
        }
        assert(u != n);
        assert(degree[u] + K > best_solution_size&&degree[u] + K > S_end);

        // the first branch includes u into S
        assert(SR[SR_rid[u]] == u&&SR_rid[u] >= S_end&&SR_rid[u] < R_end);
        swap_pos(S_end, SR_rid[u]);
        ++ S_end;

        ui pre_best_solution_size = best_solution_size, old_R_end = R_end;
        if(!move_u_to_S_with_prune(S_end, R_end, level)) kplex_search(S_end, R_end, level+1, false);
        restore_SR_and_edges(S_end, R_end, old_R_end, level);

        if(must_include) {
        	move_u_to_R_wo_prune(S_end, R_end, level);
        	return ;
        }

        // the second branch exclude u from S
        assert(Qv.empty());
        bool pruned = remove_u_from_S_with_prune(S_end, R_end, level);
        if(!pruned&&best_solution_size > pre_best_solution_size) pruned = collect_removable_vertices_and_edges(S_end, R_end, level);
        if(!pruned) {
            if(!remove_vertices_and_edges_with_prune(S_end, R_end, level)) {
            	kplex_search(S_end, R_end, level+1, false);
            }
        }
        restore_SR_and_edges(S_end, R_end, old_R_end, level);
    }

    bool move_u_to_S_with_prune(ui S_end, ui &R_end, ui level)
    {
    	assert(S_end > 0);
        ui u = SR[S_end-1];
        int *t_matrix = matrix + u*n;

        ui neighbors_n = 0, nonneighbors_n = 0;
        for(ui i = 0;i < R_end;i ++) if(i != S_end-1) {
            if(t_matrix[SR[i]]) neighbors[neighbors_n++] = SR[i];
            else nonneighbors[nonneighbors_n++] = SR[i];
        }

        for(ui i = 0;i < neighbors_n;i ++) ++ degree_in_S[neighbors[i]];

        assert(Qv.empty());
        if(degree_in_S[u] + K == S_end) { // only neighbors of u in R can be candidates --- RR2
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
            for(;i < nonneighbors_n;i ++) if(S_end - degree_in_S[nonneighbors[i]] >= K) {
            	assert(level_id[nonneighbors[i]] > level);
            	level_id[nonneighbors[i]] = level;
                Qv.push(nonneighbors[i]);
            }
        }

        // RR2
        for(ui i = 0;i < nonneighbors_n&&SR_rid[nonneighbors[i]] < S_end;i ++) if(degree_in_S[nonneighbors[i]] + K == S_end) {
            int *tt_matrix = matrix + nonneighbors[i]*n;
            for(ui j = S_end;j < R_end;j ++) if(level_id[SR[j]] > level&&!tt_matrix[SR[j]]) {
            	level_id[SR[j]] = level;
                Qv.push(SR[j]);
            }
        }
        return remove_vertices_and_edges_with_prune(S_end, R_end, level);
    }

    bool remove_vertices_and_edges_with_prune(ui S_end, ui &R_end, ui level)
    {
		while(!Qv.empty()) {
			ui u = Qv.front(); Qv.pop(); // remove u
			assert(SR[SR_rid[u]] == u);
			assert(SR_rid[u] >= S_end&&SR_rid[u] < R_end);
			-- R_end;
			swap_pos(SR_rid[u], R_end);

			bool terminate = false;
			ui neighbors_n = 0;
			int *t_matrix = matrix + u*n;
			for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) {
				ui w = SR[i];
				neighbors[neighbors_n++] = w;
				-- degree[w];
				if(degree[w] + K <= best_solution_size) {
					if(i < S_end) terminate = true; // UB1
					else if(level_id[w] > level) { // RR3
						level_id[w] = level;
						Qv.push(w);
					}
				}
			}
			// UB1
			if(terminate) {
				for(ui i = 0;i < neighbors_n;i ++) ++ degree[neighbors[i]];
				level_id[u] = n;
				++ R_end;
				return true;
			}
		}

        return false;
    }

    void restore_SR_and_edges(ui S_end, ui &R_end, ui old_R_end, ui level)
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

            ui neighbors_n = 0;
            int *t_matrix = matrix + u*n;
            for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) {
                ui w = SR[i];
                neighbors[neighbors_n ++] = w;
                ++ degree[w];
            }

            ++ R_end;
        }
    }
    
    void move_u_to_R_wo_prune(ui &S_end, ui &R_end, ui level)
    {
    	assert(S_end);
        ui u = SR[-- S_end];
        ui neighbors_n = 0;
        int *t_matrix = matrix + u*n;
        for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) neighbors[neighbors_n ++] = SR[i];
        for(ui i = 0;i < neighbors_n;i ++) -- degree_in_S[neighbors[i]];
    }

    bool remove_u_from_S_with_prune(ui &S_end, ui &R_end, ui level)
    {
    	assert(S_end);
		ui u = SR[S_end-1];
		-- S_end; -- R_end;
		swap_pos(S_end, R_end);
		level_id[u] = level;

		bool ret = false;
        ui neighbors_n = 0;
        int *t_matrix = matrix + u*n;
        for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) neighbors[neighbors_n ++] = SR[i];
        for(ui i = 0;i < neighbors_n;i ++) {
        	-- degree_in_S[neighbors[i]];
        	-- degree[neighbors[i]];
        	if(degree[neighbors[i]] + K <= best_solution_size) {
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

    bool collect_removable_vertices_and_edges(ui S_end, ui R_end, ui level) {
    	for(ui i = 0;i < S_end;i ++) if(degree[SR[i]] + K <= best_solution_size) return true;

    	for(ui i = S_end;i < R_end;i ++) if(level_id[SR[i]] > level){
    		if(S_end - degree_in_S[SR[i]] >= K||degree[SR[i]] + K <= best_solution_size) {
    			assert(level_id[SR[i]] > level);
    			level_id[SR[i]] = level;
    			Qv.push(SR[i]);
    			continue;
    		}
    		int *t_matrix = matrix + SR[i]*n;
    		for(ui j = 0;j < S_end;j ++) {
    			if(S_end - degree_in_S[SR[j]] == K&&!t_matrix[SR[j]])
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

    void swap_pos(ui i, ui j) {
        std::swap(SR[i], SR[j]);
        SR_rid[SR[i]] = i;
        SR_rid[SR[j]] = j;
    }

    ui choose_branch_vertex(ui S_end, ui R_end)
    {
        ui u = n, min_degree_in_S = n;
        for(ui i = S_end;i < R_end;i ++) {
            ui v = SR[i];
            if(degree[v] + K > best_solution_size && degree[v] + K > S_end) {
                u = v;
                break;
            }
        }
        assert(u != n);
        return u;
    }
};

#endif