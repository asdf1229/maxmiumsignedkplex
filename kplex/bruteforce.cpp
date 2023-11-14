#include <bits/stdc++.h>
#include "Timer.h"
#include "Utility.h"
#include "LinearHeap.h"
using namespace std;

vector<int> P;
// vector<int> newP;
int ans, newans;
int lb, ub;

struct Edge {
	int idx, a, b, l;
};

struct Graph {
	int n; //number of vertices
	int m; //number of edges
	int pm; //number of positive edges
	int nm; //number of negative edges

	/*Store edges in linear arrays*/
	int * pstart; //start of edge number of a point
	int * pend; //end of edge number of a point
	int * edges; //edges

	int * p_pstart; //start of positive edge number of a point
	int * p_pend; //end of positive edge number of a point
	int * p_edges; //positive edges

	int * n_pstart; //start of negative edge number of a point
	int * n_pend; //end of negative edge number of a point
	int * n_edges; //negative edges

	int * degree; //degree of a point
	int * p_degree; //positive degree of a point
	int * n_degree; //negative degree of a point

	int ** Matrix;

	Edge * edges_pair;

	//Sinttable for situations with small data scales
	int * tri_cnt; 
}G, g;

int * bit_del; //mark whether the node is deleted
int * bit_sel; //mark whether the node is selected
int * edge_del;
int * deleted;

void load_graph(string input_graph)
{
	Timer t;
	t.restart();
    string buffer;
    ifstream input_file(input_graph, ios::in);

    if (!input_file.is_open()) {
        cout<<"cannot open file : "<<input_graph<<endl;exit(1);
    }
	else {
		input_file >> G.n >> G.m;
		map<int, int> * s_G = new map<int, int>[G.n];
		int u, v;
		int flag;
		while (input_file >> u >> v >> flag) {
			if (u == v) continue;
            assert(u >= 0 && u < n);
            assert(v >= 0 && v < n);
            assert(flag == 1 || flag == -1);
            s_G[u].insert(make_pair(v, flag));
            s_G[v].insert(make_pair(u, flag));
		}
		G.m = 0; G.pm = 0; G.nm = 0;
		for (int i = 0; i < G.n; i++) {
			const map<int, int> & nei = s_G[i];
			for (auto e : nei) {
				if (e.second == 1)
					++G.pm;
				else
					++G.nm;
			}
			G.m += nei.size();
		}
		assert(m%2 == 0);assert(pm%2 == 0);assert(nm%2 == 0);
		G.m /= 2; G.pm /= 2; G.nm /= 2;

		input_file.close();

		G.pstart = new int[G.n+1];
		G.pend = new int[G.n];
		G.edges = new int[2*G.m];
		G.edges_pair = new Edge[2*G.m];

		G.p_pstart = new int[G.n+1];
        G.p_pend = new int[G.n];
		G.p_edges = new int[2*G.pm];
        G.n_pstart = new int[G.n+1];
        G.n_pend = new int[G.n];
        G.n_edges = new int[2*G.nm];

        G.degree = new int[G.n];
        G.p_degree = new int[G.n];
        G.n_degree = new int[G.n];

		{
			g.pstart = new int[G.n+1];
			g.pend = new int[G.n];
			g.edges = new int[2*G.m];
			g.p_pstart = new int[G.n+1];
			g.p_pend = new int[G.n];
			g.p_edges = new int[2*G.pm];
			g.n_pstart = new int[G.n+1];
			g.n_pend = new int[G.n];
			g.n_edges = new int[2*G.nm];

			g.degree = new int[G.n];
			g.p_degree = new int[G.n];
			g.n_degree = new int[G.n];
		}

		//construct positive edges
		G.p_pstart[0] = 0;
		for (int i = 0; i < G.n; i++) {
			const map<int, int> & nei = s_G[i];
			int start_idx = G.p_pstart[i];
			int d = 0;
			for (auto e : nei) {
				if (e.second == 1) {
					G.p_edges[start_idx++] = e.first;
					d++;
				}
			}
			G.p_pstart[i+1] = start_idx;
			G.p_degree[i] = d;
		}
		assert(G.p_pstart[G.n] == 2*G.pm);

		//construct negative edges
		G.n_pstart[0] = 0;
		for (int i = 0; i < G.n; i++) {
			const map<int, int> & nei = s_G[i];
			int start_idx = G.n_pstart[i];
			int d = 0;
			for (auto e : nei) {
				if (e.second == -1) {
					G.n_edges[start_idx++] = e.first;
					d++;
				}
			}
			G.n_pstart[i+1] = start_idx;
			G.n_degree[i] = d;
		}
		assert(G.n_pstart[G.n] == 2*G.nm);

		for (int i = 0; i < G.n; i++) {
            G.p_pend[i] = G.p_pstart[i+1];
            G.n_pend[i] = G.n_pstart[i+1];
		}

		//construct edges
		G.pstart[0] = 0;
		for (int u = 0; u < G.n; u++) {
			int start_idx = G.pstart[u];
			for (int j = G.p_pstart[u]; j < G.p_pend[u]; j++) {
				int v = G.p_edges[j];
				G.edges_pair[start_idx] = Edge{start_idx, u, v, 1};
				G.edges[start_idx++] = v;
			}
			for (int j = G.n_pstart[u]; j < G.n_pend[u]; j++) {
				int v = G.n_edges[j];
				G.edges_pair[start_idx] = Edge{start_idx, u, v, -1};
				G.edges[start_idx++] = v;
			}
			G.pend[u] = start_idx;
			G.pstart[u + 1] = start_idx;
		}
		assert(G.pstart[G.n] == 2*G.m);

		for (int i = 0; i < G.n; i++) {
			G.degree[i] = G.p_degree[i] + G.n_degree[i];
		}
		delete [] s_G;
	}

	cout<<"\t load_graph: time cost = "<<integer_to_string(t.elapsed())<<endl;
	cout<<"\t G : n = "<<G.n<<", m = "<<G.m<<", pm = "<<G.pm<<", nm = "<<G.nm<<endl;
}

// degeneracy-based k-plex
// return an upper bound of the maximum k-plex size
// return dOrder
void kplex_degen(ListLinearHeap *heap, int k, int *dOrder) {
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

void kplex_hec(int k, int *dOrder) {
	// P.clear();
	// for(int i = G.n - 1; i >= 0; i--) {
	// 	int u = dOrder[i];
	// 	int nowsize = P.size();
	// 	for(int j = 0; j < nowsize; j++) {
	// 		for(int k = j + 1; k < nowsize; k++) {

	// 		}
	// 	}
	// }
}

//get G's k-core
void get_G_core(int k)
{
	printf("\t k = %d\n", k);
	Timer t;
	if(k < 2) { //threshold should be at least 2
		cout<<"\t get_G_core, T : "<<integer_to_string(t.elapsed())<<",\t n="<<G.n<<", m="<<G.m<<endl;
		return;
	}
	t.restart();
	int threshold = k;
	int del_count = 0;
	bit_del = new int[G.n];
	memset(bit_del, 0, sizeof(int)*G.n);
	queue<int> q;

	for (int i = 0; i < G.n; i++) if(G.degree[i] < threshold) q.push(i);
	while (!q.empty()) {
		int u = q.front();
		q.pop();
		del_count++;
		bit_del[u] = 1;
        for (int i = G.p_pstart[u]; i < G.p_pend[u]; i++) {
			int v = G.p_edges[i];
			G.p_degree[v]--;
			if (G.p_degree[v] + G.n_degree[v] == threshold - 1) {
				q.push(v);
			}
        }
        for (int i = G.n_pstart[u]; i < G.n_pend[u]; i++) {
			int v = G.n_edges[i];
			G.n_degree[v]--;
			if (G.p_degree[v] + G.n_degree[v] == threshold - 1) {
				q.push(v);
			}
        }
	}


	//rebuild
	int * mapping = new int[G.n];
	int idx = 0;
	for (int i = 0; i < G.n; i++) {
		if (!bit_del[i]) {
			mapping[i] = idx++;
		}
	}

    int * t_p_pstart = new int[G.n+1];
    int * t_p_edges = new int[2*G.pm];
    int * t_n_pstart = new int[G.n+1];
    int * t_n_edges = new int[2*G.nm];
	
	int new_i = 0;
	t_p_pstart[0] = 0;
	for (int i = 0; i < G.n; i++) if (!bit_del[i]) {
		int start_idx = t_p_pstart[new_i];
		for (int j = G.p_pstart[i]; j < G.p_pend[i]; j++) {
			int v = G.p_edges[j];
			if (!bit_del[v]) {
				t_p_edges[start_idx++] = mapping[v];
			}
		}
		t_p_pstart[++new_i] = start_idx;
	}
	assert(new_i == idx);
	new_i = 0;
	t_n_pstart[0] = 0;
	for (int i = 0; i < G.n; i++) if (!bit_del[i]) {
		int start_idx = t_n_pstart[new_i];
		for (int j = G.n_pstart[i]; j < G.n_pend[i]; j++) {
			int v = G.n_edges[j];
			if (!bit_del[v]) {
				t_n_edges[start_idx++] = mapping[v];
			}
		}
		t_n_pstart[++new_i] = start_idx;
	}
	assert(new_i == idx);
	G.n = idx;

    delete [] G.p_pstart;
    delete [] G.p_edges;
    delete [] G.n_pstart;
    delete [] G.n_edges;
    delete [] mapping;
	delete [] bit_del;

    G.p_pstart = t_p_pstart;
    G.p_edges = t_p_edges;
    G.n_pstart = t_n_pstart;
    G.n_edges = t_n_edges;

    for (int i = 0; i < G.n; i++){
        G.p_pend[i] = G.p_pstart[i+1];
        G.n_pend[i] = G.n_pstart[i+1];
    }

	//construct edges
	G.pstart[0] = 0;
	for (int u = 0; u < G.n; u++) {
		int start_idx = G.pstart[u];
		for (int j = G.p_pstart[u]; j < G.p_pend[u]; j++) {
			int v = G.p_edges[j];
			G.edges_pair[start_idx] = Edge{start_idx, u, v, 1};
			G.edges[start_idx++] = v;
		}
		for (int j = G.n_pstart[u]; j < G.n_pend[u]; j++) {
			int v = G.n_edges[j];
			G.edges_pair[start_idx] = Edge{start_idx, u, v, -1};
			G.edges[start_idx++] = v;
		}
		G.pend[u] = start_idx;
		G.pstart[u + 1] = start_idx;
	}
	assert(G.pstart[G.n] == 2*G.m);

	//rebuild degree
	for (int i = 0; i < G.n; i++) {
        G.p_degree[i] = G.p_pstart[i+1] - G.p_pstart[i];
        G.n_degree[i] = G.n_pstart[i+1] - G.n_pstart[i];
        G.degree[i] = G.p_degree[i] + G.n_degree[i];
	}

    int now_m = 0, now_pm = 0, now_nm = 0;
    for(int i = 0; i < G.n; i++) {now_m += G.degree[i]; now_pm += G.p_degree[i]; now_nm += G.n_degree[i];}
    assert(now_m%2 == 0); assert(now_pm%2 == 0); assert(now_nm%2 == 0);
    now_m /= 2; now_pm /= 2; now_nm /= 2;
	G.m = now_m; G.pm = now_pm; G.nm = now_nm;
	
	cout<<"\t get_G_core, T : "<<integer_to_string(t.elapsed())<<",\t n="<<G.n<<", m="<<G.m<<endl;
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

//get triangle count for all edges of G
void get_G_tricnt()
{
	long long cnt = 0;
	int *adj = new int[G.n];
	G.tri_cnt = new int[G.m*2];

	int *rid = adj;
	for(int i = 0; i < G.n; i++) rid[i] = i;
	for(int i = 0; i < G.n; i++) {
		int &end = G.pend[i] = G.pstart[i];
		for(int j = G.pstart[i]; j < G.p_pstart[i + 1]; j++) {
			if(rid[G.edges[j]] > rid[i]) {
				G.edges[end++] = G.edges[j];
			}
		}
	}

	memset(adj, 0, sizeof(int)*G.n);
	memset(G.tri_cnt, 0, sizeof(int)*G.m);
	for(int u = 0; u < G.n; u++) {
		for(int j = G.pstart[u]; j < G.pend[u]; j++) adj[G.edges[j]] = j+1;

		for(int j = G.pstart[u]; j < G.pend[u]; j++) {
			int v = G.edges[j];
			for(int k = G.pstart[v]; k < G.pend[v]; k++) if(adj[G.edges[k]]) {
				G.tri_cnt[j]++;
				G.tri_cnt[k]++;
				G.tri_cnt[adj[G.edges[k]]-1]++;
				cnt++;
			}
		}
	}

	delete [] adj;
}

void get_g(int u)
{
	Timer t;
	bit_sel = new int[G.n];
	memset(bit_sel, 0, sizeof(int)*G.n);
    bit_sel[u] = 1;

	for(int i = G.p_pstart[u]; i < G.p_pend[u]; i++) {
		int v1 = G.p_edges[i];
		if(deleted[v1]) continue;
		bit_sel[v1] = 1;
		for(int j = G.p_pstart[v1]; j < G.p_pend[v1]; j++) {
			int v2 = G.p_edges[j];
			if(deleted[v2]) continue;
			bit_sel[v2] = 1;
		}
		for(int j = G.n_pstart[v1]; j < G.n_pend[v1]; j++) {
			int v2 = G.n_edges[j];
			if(deleted[v2]) continue;
			bit_sel[v2] = 1;
		}
	}
	for(int i = G.n_pstart[u]; i < G.n_pend[u]; i++) {
		int v1 = G.n_edges[i];
		if(deleted[v1]) continue;
		bit_sel[v1] = 1;
		for(int j = G.p_pstart[v1]; j < G.p_pend[v1]; j++) {
			int v2 = G.p_edges[j];
			if(deleted[v2]) continue;
			bit_sel[v2] = 1;
		}
		for(int j = G.n_pstart[v1]; j < G.n_pend[v1]; j++) {
			int v2 = G.n_edges[j];
			if(deleted[v2]) continue;
			bit_sel[v2] = 1;
		}
	}

	//build g
	int * mapping = new int[G.n];
	int idx = 0;
	for (int i = 0; i < G.n; i++) {
		if(deleted[i]) continue;
		if (bit_sel[i]) {
			mapping[i] = idx++;
		}
	}
	
	int new_i = 0;
	g.p_pstart[0] = 0;
	for (int i = 0; i < G.n; i++) if (bit_sel[i]) {
		int start_idx = g.p_pstart[new_i];
		for (int j = G.p_pstart[i]; j < G.p_pend[i]; j++) {
			int v = G.p_edges[j];
			if (bit_sel[v]) {
				g.p_edges[start_idx++] = mapping[v];
			}
		}
		g.p_pstart[++new_i] = start_idx;
	}
	assert(new_i == idx);
	new_i = 0;
	g.n_pstart[0] = 0;
	for (int i = 0; i < G.n; i++) if (bit_sel[i]) {
		int start_idx = g.n_pstart[new_i];
		for (int j = G.n_pstart[i]; j < G.n_pend[i]; j++) {
			int v = G.n_edges[j];
			if (bit_sel[v]) {
				g.n_edges[start_idx++] = mapping[v];
			}
		}
		g.n_pstart[++new_i] = start_idx;
	}
	assert(new_i == idx);
	g.n = idx;

	delete [] mapping;
	delete [] bit_sel;

    for (int i = 0; i < g.n; i++){
        g.p_pend[i] = g.p_pstart[i+1];
        g.n_pend[i] = g.n_pstart[i+1];
    }

	//build degree
	for (int i = 0; i < g.n; i++) {
        g.p_degree[i] = g.p_pstart[i+1] - g.p_pstart[i];
        g.n_degree[i] = g.n_pstart[i+1] - g.n_pstart[i];
        g.degree[i] = g.p_degree[i] + g.n_degree[i];
	}

    int now_m = 0, now_pm = 0, now_nm = 0;
    for(int i = 0; i < g.n; i++) {now_m += g.degree[i]; now_pm += g.p_degree[i]; now_nm += g.n_degree[i];}
    assert(now_m%2 == 0); assert(now_pm%2 == 0); assert(now_nm%2 == 0);
    now_m /= 2; now_pm /= 2; now_nm /= 2;
	g.m = now_m; g.pm = now_pm; g.nm = now_nm;

	//build adjacent matrix
	/// ?
	// for(int i = 0; i < g.n; i++) {
	// 	if(g.Matrix[i] != NULL) {
	// 		delete [] g.Matrix[i];
	// 		g.Matrix[i] = NULL;
	// 	}
	// }
	// if(g.Matrix != NULL) {
	// 	delete [] g.Matrix;
	// 	g.Matrix = NULL;
	// }

	g.Matrix = new int*[g.n];
	for(int i = 0; i < g.n; i++) g.Matrix[i] = new int[g.n];
	for(int i = 0; i < g.n; i++)
		for(int j = 0; j < g.n; j++)
			g.Matrix[i][j] = 0;

	for(int i = 0; i < g.n; i++) {
		int u = i;
		for(int j = g.p_pstart[i]; j < g.p_pend[i]; j++) {
			int v = g.p_edges[j];
			g.Matrix[u][v] = 1;
		}
		for(int j = g.n_pstart[i]; j < g.n_pend[i]; j++) {
			int v = g.n_edges[j];
			g.Matrix[u][v] = -1;
		}
	}
	// for(int i = 0; i < g.n; i++) {
	// 	for(int j = 0; j < g.n; j++) {
	// 		printf("%d ", g.Matrix[i][j]);
	// 	}
	// 	printf("\n");
	// }

	cout<<"\t get_g, T : "<<integer_to_string(t.elapsed())<<",\t n="<<g.n<<", m="<<g.m<<endl;
}

void truss_peeling(queue<int> *qv, queue<Edge> *qe, int tv, int te) {
	while(!qv->empty()) {
		int u = qv->front(); qv->pop();
		for(int i = G.pstart[u]; i < G.pend[u]; i++) {
			int v = G.edges[i];
			if(deleted[v]) continue;
			G.degree[v]--;
			if(G.degree[v] == tv - 1) qv->push(v);
		}
		for(int i = G.p_pstart[u]; i < G.p_pend[u]; i++) {
			int v = G.p_edges[i];
			if(deleted[v]) continue;
			G.p_degree[v]--;
		}
		for(int i = G.n_pstart[u]; i < G.n_pend[u]; i++) {
			int v = G.n_edges[i];
			if(deleted[v]) continue;
			G.n_degree[v]--;
		}
		deleted[u] = 1;
	}
	// while(!qv.empty() && !qe.empty()) {
		// 	while(!qe.empty()) {

		// 	}
	// }
}

//core-truss co-pruning
void CTCP(int del, int lb_changed, int tv, int te) {
	queue<int> qv;
	queue<Edge> qe;
	if(del != -1) qv.push(del);
	// if(lb_changed) {
	// 	for(int i = 0; i < G.m; i++) {
	// 		if(G.tri_cnt[i] < te && !del[i]) {
	// 			qe.push(G.edges_pair[i]);
	// 		}
	// 	}
	// }
	truss_peeling(&qv, &qe, tv, te);
	while(!qv.empty()) {
		truss_peeling(&qv, &qe, tv, te);
	}
}

int cal(int k, int siz, vector<int> s) {
	int * deg = new int[siz];
	memset(deg, 0, sizeof(int)*siz);

    for(int i = 0; i < siz; i++) {
        for(int j = i + 1; j < siz; j++) {
            int u = s[i], v = s[j];
            if(g.Matrix[u][v]) {
                deg[i]++;
                deg[j]++;
            }
        }
    }
    for(int i = 0; i < siz; i++) {
        if(deg[i] < siz - k) return 0;
    }
	delete [] deg;

	for(int i = 0; i < siz; i++) {
		int u = s[i];
		for(int j = i + 1; j < siz; j++) {
			int v = s[j];
			if(g.Matrix[u][v]) {
				for(int k = j + 1; k < siz; k++) {
					int w = s[k];
					if(g.Matrix[v][w] && g.Matrix[u][w]) {
						int tmp = g.Matrix[u][v] + g.Matrix[v][w] + g.Matrix[u][w];
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

void kplex_search(int dep, int siz, int k, vector<int> s)
{
    if(dep == g.n) {
		if(s.size() <= newans || s.size() <= lb) return;
        if(cal(k, siz, s)) {
			for(auto u : s) {
				printf("%d ", u);
			}
			printf("\n");
            newans = max(newans, siz);
        }
        return;
    }
    for(int i = dep; i < g.n; i++) {
        int u = i;
        //reduce
        kplex_search(i + 1, siz, k, s);

        s.push_back(u);
        kplex_search(i + 1, siz + 1, k, s);
        s.pop_back();
    }

}

int main(int argc, const char * argv[])
{
    if(argc < 2) {
        cout<<"\t Usage: [0]exe [1]input_graph [2]k (optional)\t"<<endl; exit(1);
    }

	load_graph(argv[1]);

	int k; //size constraint
    k = 3;
    if(argc > 2) k = atoi(argv[2]);	

	cout<<"\t Graph: "<<argv[1]<<",\t k: "<<k<<endl;

	//(lb, ub) <- kplex_degen(g, k)
	ListLinearHeap *heap = new ListLinearHeap(G.n, G.n-1);
	int *dOrder = new int[G.n];
	kplex_degen(heap, k, dOrder);
	delete heap;
	printf("\t degen: lb = %d, ub = %d\n", lb, ub);
	kplex_hec(k, dOrder);

	deleted = new int[G.n];
	memset(deleted, 0, sizeof(int)*G.n);

	if((int)P.size() < ub) {
		lb = max((int)P.size(), 2*k-2);
		get_G_core(lb+1-k);
		get_G_deg();
		get_G_tricnt();
		CTCP(-1, 1, lb+1-k, lb+1-2*k);

		while(G.n > 0) {
		// for(int i = 0; i < G.n; i++) {
			//get u
			// int u = i;
			int u = -1;
			for(int i = 0; i < G.n; i++) {
				if(deleted[i]) continue;
				if(u == -1 || G.degree[u] > G.degree[i]) {
					u = i;
				}
			}
			printf("u = %d\n", u);
			if(u == -1) break;
			//get g
			get_g(u);
			//kplex
			vector<int> s; s.clear();
			kplex_search(0, 0, k, s);
			int lb_changed = 0;
			if(newans >= lb) {
				ans = newans;
				lb = newans;
				if(newans > lb) lb_changed = 1;
			}
			// CTCP
			CTCP(u, lb_changed, lb+1-k, lb+1-2*k);
		}
	}

	printf("ans = %d\n", ans);
	// delete_memo();
    return 0;
}
