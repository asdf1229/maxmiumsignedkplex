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

    int s_n;

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
    void get_g(ui u, vector<pair<int,int> > &vp, vector<int> &sgn);
    void rebuild_graph(bool *v_del);
    void CTCP(int del_v, bool lb_changed, int tv, int te);
    void heu_signed_kplex(int rounds, int k);
};

#endif