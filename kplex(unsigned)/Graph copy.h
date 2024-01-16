#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "Utility.h"
#include "Timer.h"
#include "LinearHeap.h"

class Graph {
private:
    vector<ui> kplex;

    int K;

    ui n; //number of vertices
    ept m; //number of edges
    ept pm; //number of positive edges
    ept nm; //number of negative edges

    /*Store edges in linear arrays*/
    ept *pstart; //start of edge number of a point
    ept *pend; //end of edge number of a point
    ui *edges; //edges

	ept *p_pstart; //start of positive edge number of a point
	ept *p_pend; //end of positive edge number of a point
	ui *p_edges; //positive edges

	ept *n_pstart; //start of negative edge number of a point
	ept *n_pend; //end of negative edge number of a point
	ui *n_edges; //negative edges

	ui *degree; //degree of a point
	ui *p_degree; //positive degree of a point
	ui *n_degree; //negative degree of a point
    ept *tri_cnt;
    int lb, ub;


    // int max_core; //length of s_matrix
    int s_n;
    int **s_matrix; //adjacent matrix of subgraph
    ui s_solution_size;
    ui *s_solution;
    queue<ui> Qv;
    ui *s_degree_in_S;
    ui *s_degree;
    ui *SR;
    ui *SR_rid;
    ui *neighbors;
    ui *nonneighbors;
    ui *level_id;

public:
	Graph(const int _k);
	~Graph();

    void load_graph(string input_graph);
    void swap_pos(ui i, ui j);
    void subgraph_init();
    void find_signed_kplex();

private:
    void get_k_core(int k);
    void get_degree();
    void get_tricnt();
    void get_g(ui u);
    void rebuild_graph(bool *v_del);
    void CTCP(int del_v, bool lb_changed, ui tv, ui te);
    void heu_signed_kplex(int rounds, int k);
    int cal(int k, int siz, vector<ui> s);
    void matrix_kplex();
    bool move_u_to_S_with_prune(ui S_end, ui &R_end, ui level);
    bool remove_vertices_and_edges_with_prune(ui S_end, ui &R_end, ui level);
    void restore_SR_and_edges(ui S_end, ui &R_end, ui old_R_end, ui level);
    void move_u_to_R_wo_prune(ui &S_end, ui &R_end, ui level);
    bool remove_u_from_S_with_prune(ui &S_end, ui &R_end, ui level);
    bool collect_removable_vertices_and_edges(ui S_end, ui R_end, ui level);
    ui choose_branch_vertex_based_on_non_neighbors(ui S_end, ui R_end);
    void kplex_search(ui S_end, ui R_end, ui dep, bool choose_zero);
    void s_init(ui &R_end);
};

#endif