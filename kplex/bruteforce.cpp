#include <bits/stdc++.h>
#include "Timer.h"
#include "Utility.h"
#include "LinearHeap.h"
using namespace std;

vector<int> P;
// vector<int> newP;
int ans, newans;
int lb, ub = 10000;

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
	// int * tri_cnt;
	map<int, int> * tri_cnt;

	//orient graph
	int * o_pstart;
	int * o_pend;
	int * o_edges;

	int ** Matrix;
}G, g;

int * v_del; //mark whether the node is deleted
int * v_sel; //mark whether the node is selected
// int * e_del;
int * mark;
map<int, int> * e_del;
// map<int, int> * e_idx;

void init_hash()
{
	// e_idx = new map<int, int>[G.n];
	e_del = new map<int, int>[G.n];
	G.tri_cnt = new map<int, int>[G.n];
    for(int u = 0; u < G.n; u++){
        for(int i = G.pstart[u]; i < G.pend[u]; i++){
			int v = G.edges[i];
			if(u < v) {
				// e_idx[u][v] = i; e_idx[v][u] = i;
				e_del[u][v] = 0; e_del[v][u] = 0;
				G.tri_cnt[u][v] = 0; G.tri_cnt[v][u] = 0;
			}
        }
    }
}

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
			G.p_pend[i] = G.p_pstart[i+1] = start_idx;
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
			G.n_pend[i] = G.n_pstart[i+1] = start_idx;
			G.n_degree[i] = d;
		}
		assert(G.n_pstart[G.n] == 2*G.nm);

		//construct edges
		G.pstart[0] = 0;
		for (int u = 0; u < G.n; u++) {
			int start_idx = G.pstart[u];
			for (int j = G.p_pstart[u]; j < G.p_pend[u]; j++) {
				int v = G.p_edges[j];
				G.edges[start_idx++] = v;
			}
			for (int j = G.n_pstart[u]; j < G.n_pend[u]; j++) {
				int v = G.n_edges[j];
				G.edges[start_idx++] = v;
			}
			G.pend[u] = start_idx;
			G.pstart[u + 1] = start_idx;
			G.degree[u] = G.p_degree[u] + G.n_degree[u];
		}
		assert(G.pstart[G.n] == 2*G.m);

		delete [] s_G;
	}

	cout<<"\t load_graph: time cost = "<<integer_to_string(t.elapsed())<<endl;
	cout<<"\t G : n = "<<G.n<<", m = "<<G.m<<", pm = "<<G.pm<<", nm = "<<G.nm<<endl;
}

//get G's k-core
void get_G_core(int k)
{
	// printf("\t k = %d\n", k);
	Timer t;
	if(k < 2) { //threshold should be at least 2
		cout<<"\t get_G_core, T : "<<integer_to_string(t.elapsed())<<",\t n="<<G.n<<", m="<<G.m<<endl;
		return;
	}
	t.restart();
	int threshold = k;
	int del_count = 0;
	memset(v_del, 0, sizeof(int)*G.n);
	queue<int> q;

	for (int i = 0; i < G.n; i++) if(G.degree[i] < threshold) q.push(i);
	while (!q.empty()) {
		int u = q.front();
		q.pop();
		del_count++;
		v_del[u] = 1;
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
		if (!v_del[i]) {
			mapping[i] = idx++;
		}
	}

    int * t_p_pstart = new int[G.n+1];
    int * t_p_edges = new int[2*G.pm];
    int * t_n_pstart = new int[G.n+1];
    int * t_n_edges = new int[2*G.nm];
	
	int new_i = 0;
	t_p_pstart[0] = 0;
	for (int i = 0; i < G.n; i++) if (!v_del[i]) {
		int start_idx = t_p_pstart[new_i];
		for (int j = G.p_pstart[i]; j < G.p_pend[i]; j++) {
			int v = G.p_edges[j];
			if (!v_del[v]) {
				t_p_edges[start_idx++] = mapping[v];
			}
		}
		t_p_pstart[++new_i] = start_idx;
	}
	assert(new_i == idx);
	new_i = 0;
	t_n_pstart[0] = 0;
	for (int i = 0; i < G.n; i++) if (!v_del[i]) {
		int start_idx = t_n_pstart[new_i];
		for (int j = G.n_pstart[i]; j < G.n_pend[i]; j++) {
			int v = G.n_edges[j];
			if (!v_del[v]) {
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
	memset(v_del, 0, sizeof(int)*G.n);

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
			G.edges[start_idx++] = v;
		}
		for (int j = G.n_pstart[u]; j < G.n_pend[u]; j++) {
			int v = G.n_edges[j];
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

//get triangle count for all edges of G
void get_G_tricnt()
{
	mark = new int[G.n];

	// G.o_pstart = new int[G.n+1];
	// G.o_pend = new int[G.n];
	// G.o_edges = new int[G.m*2];

	// G.o_pstart[0] = 0;
	// int *rid = adj;
	// for(int u = 0; u < G.n; u++) rid[u] = u;
	// for(int u = 0; u < G.n; u++) {
	// 	int &end = G.o_pend[u] = G.o_pstart[u];
	// 	for(int i = G.pstart[u]; i < G.pend[u]; i++) {
	// 		int v = G.edges[i];
	// 		if(u < v) {
	// 			G.o_edges[end++] = v;
	// 		}
	// 	}
	// 	G.o_pstart[u+1] = G.o_pend[u];
	// }
	
	memset(mark, 0, sizeof(int)*G.n);
	for(int u = 0; u < G.n; u++) {
		for(int j = G.pstart[u]; j < G.pend[u]; j++) mark[G.edges[j]] = 1;

		for(int j = G.pstart[u]; j < G.pend[u]; j++) {
			int v = G.edges[j];
			for(int k = G.pstart[v]; k < G.pend[v]; k++) {
				int w = G.edges[k];
				if(mark[w] == 1 && v < w) {
					G.tri_cnt[u][v]++; G.tri_cnt[v][u]++;
					G.tri_cnt[u][w]++; G.tri_cnt[w][u]++;
					G.tri_cnt[v][w]++; G.tri_cnt[w][v]++;
				}
			}
		}

		for(int j = G.pstart[u]; j < G.pend[u]; j++) mark[G.edges[j]] = 0;
	}

	delete [] mark;
}

int get_g(int u)
{
	Timer t;
	v_sel = new int[G.n];
	memset(v_sel, 0, sizeof(int)*G.n);
    v_sel[u] = 1;

	for(int i = G.p_pstart[u]; i < G.p_pend[u]; i++) {
		int v1 = G.p_edges[i];
		if(v_del[v1]) continue;
		v_sel[v1] = 1;
		for(int j = G.p_pstart[v1]; j < G.p_pend[v1]; j++) {
			int v2 = G.p_edges[j];
			if(v_del[v2]) continue;
			v_sel[v2] = 1;
		}
		for(int j = G.n_pstart[v1]; j < G.n_pend[v1]; j++) {
			int v2 = G.n_edges[j];
			if(v_del[v2]) continue;
			v_sel[v2] = 1;
		}
	}
	for(int i = G.n_pstart[u]; i < G.n_pend[u]; i++) {
		int v1 = G.n_edges[i];
		if(v_del[v1]) continue;
		v_sel[v1] = 1;
		for(int j = G.p_pstart[v1]; j < G.p_pend[v1]; j++) {
			int v2 = G.p_edges[j];
			if(v_del[v2]) continue;
			v_sel[v2] = 1;
		}
		for(int j = G.n_pstart[v1]; j < G.n_pend[v1]; j++) {
			int v2 = G.n_edges[j];
			if(v_del[v2]) continue;
			v_sel[v2] = 1;
		}
	}

	//build g
	int nowu = -1;
	int * mapping = new int[G.n];
	int idx = 0;
	for (int i = 0; i < G.n; i++) {
		if(v_del[i]) continue;
		if(v_sel[i]) {
			if(u == i) nowu = i;
			mapping[i] = idx++;
		}
	}
	
	int new_i = 0;
	g.p_pstart[0] = 0;
	for (int i = 0; i < G.n; i++) if (v_sel[i]) {
		int start_idx = g.p_pstart[new_i];
		for (int j = G.p_pstart[i]; j < G.p_pend[i]; j++) {
			int v = G.p_edges[j];
			if (v_sel[v]) {
				g.p_edges[start_idx++] = mapping[v];
			}
		}
		g.p_pstart[++new_i] = start_idx;
	}
	assert(new_i == idx);
	new_i = 0;
	g.n_pstart[0] = 0;
	for (int i = 0; i < G.n; i++) if (v_sel[i]) {
		int start_idx = g.n_pstart[new_i];
		for (int j = G.n_pstart[i]; j < G.n_pend[i]; j++) {
			int v = G.n_edges[j];
			if (v_sel[v]) {
				g.n_edges[start_idx++] = mapping[v];
			}
		}
		g.n_pstart[++new_i] = start_idx;
	}
	assert(new_i == idx);
	g.n = idx;

	delete [] mapping;
	delete [] v_sel;

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
	return nowu;
}

void truss_peeling(queue<int> *qv, queue<pair<int, int>> *qe, int tv, int te)
{
	memset(mark, 0, sizeof(int)*G.n);
	while(!qv->empty()) {
		int u = qv->front(); qv->pop();
		for(int i = G.pstart[u]; i < G.pend[u]; i++) {
			int v = G.edges[i];
			if(!e_del[u][v] && !v_del[v]) {
				if((G.degree[v]--) == tv) qv->push(v);
			}
			e_del[u][v] = 1; e_del[v][u] = 1;
		}

		for(int i = G.pstart[u]; i < G.pend[u]; i++) {
			int v = G.edges[i];
			if(!e_del[u][v] && !v_del[v]) {
				mark[v] = 1;
			}
		}
		for(int i = G.pstart[u]; i < G.pend[u]; i++) {
			int v = G.edges[i];
			if(e_del[u][v] || v_del[v]) continue;
			for(int j = G.pstart[v]; j < G.pend[v]; j++) {
				int w = G.edges[j];
				if(!e_del[v][w] && !v_del[w]) {
					if(mark[w] == 1 && w > v) {
						if((G.tri_cnt[u][v]--) == te) qe->push(make_pair(u, v));
						if((G.tri_cnt[u][w]--) == te) qe->push(make_pair(v, w));
						if((G.tri_cnt[v][w]--) == te) qe->push(make_pair(u, w));
						G.tri_cnt[v][u]--;
						G.tri_cnt[w][u]--;
						G.tri_cnt[w][v]--;
					}
				}
			}
		}
		for(int i = G.pstart[u]; i < G.pend[u]; i++) {
			int v = G.edges[i];
			if(!e_del[u][v] && !v_del[v]) {
				mark[v] = 0;
			}
		}
		v_del[u] = 1;
	}
	while(!qe->empty()) {
		pair<int, int> ue = qe->front(); qe->pop();
		int u = ue.first;
		int v = ue.second;
		if(v_del[u] || v_del[v] || e_del[u][v]) continue;
		for(int i = G.pstart[u]; i < G.pend[u]; i++){
			int w = G.edges[i];
			if(!v_del[w] && !e_del[u][v]) {
				if(e_del[v].find(w) != e_del[v].end()) {
					if((G.tri_cnt[u][w]--) == te) qe->push(make_pair(u, w));
					if((G.tri_cnt[v][w]--) == te) qe->push(make_pair(v, w));
					G.tri_cnt[w][u]--;
					G.tri_cnt[w][v]--;
				}
			}
		}
		e_del[u][v] = 1; e_del[v][u] = 1;
	}
}

//core-truss co-pruning
void CTCP(int del, int lb_changed, int tv, int te)
{
	printf("\t CTCP: tv = %d, te = %d\n", tv, te);
	queue<int> qv;
	queue<pair<int, int>> qe;
	mark = new int[G.n];
	if(del != -1) qv.push(del);
	if(lb_changed) {
		for(int u = 0; u < G.n; u++) if(!v_del[u]) {
			for(int i = G.pstart[u]; i < G.pend[u]; i++) {
				int v = G.edges[i];
				if(u < v && !e_del[u][v] && !v_del[v]) {
					if(G.tri_cnt[u][v] < te) {
						qe.push(make_pair(u, v));
					}
				}
			}
		}
	}
	truss_peeling(&qv, &qe, tv, te);
	while(!qv.empty()) {
		truss_peeling(&qv, &qe, tv, te);
	}
	delete [] mark;
}

int cal(int k, int siz, vector<int> s)
{
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

void heu_signed_kplex(int rounds, int k)
{
	if(rounds < 1) return;
	priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> kset;
	for(int i = 0; i < rounds; i++) kset.push(make_pair(G.degree[i], i));

    for(int i = rounds; i < G.n; i++){
        if(G.degree[i] > kset.top().first){
            kset.pop();
            kset.push(make_pair(G.degree[i], i));
        }
    }
    vector<pair<int, int>> ordV(rounds);
    for(int i = 0; i < rounds; i++){
        ordV[i] = kset.top();
        kset.pop();
    }
    assert(kset.empty());
	sort(ordV.begin(), ordV.end(), greater<pair<int, int>>()); //decreasing order

    int * label = new int[G.n];
    int * vs_deg = new int[G.n];
	int * inR = new int[G.n];
	
	for(int round = 0; round < rounds && round < G.n; round++) {
        int u = ordV[round].second;

		u = get_g(u);
		
        memset(label, 0, sizeof(int)*g.n);
		memset(inR, 0, sizeof(int)*g.n);
        vector<int> res;
        res.push_back(u);
		inR[u] = 1;
        vector<int> vsP, vsN;
        for(int i = g.p_pstart[u]; i < g.p_pend[u]; i++){
            int v = g.p_edges[i];
            vsP.push_back(v);
            label[v] = 1;
        }
        for(int i = g.n_pstart[u]; i < g.n_pend[u]; i++){
            int v = g.n_edges[i];
            vsN.push_back(v);
            label[v] = 2;
        }
        for(auto e : vsP) vs_deg[e] = 0;
        for(auto e : vsN) vs_deg[e] = 0;
        for(auto e : vsP){
            for(int i = g.p_pstart[e]; i < g.p_pend[e]; i++){
                int v = g.p_edges[i];
                if(label[v] == 1) ++ vs_deg[e];
            }
            for(int i = g.n_pstart[e]; i < g.n_pend[e]; i++){
                int v = g.n_edges[i];
                if(label[v] == 2) ++ vs_deg[e];
            }
        }
        for(auto e : vsN){
            for(int i = g.p_pstart[e]; i < g.p_pend[e]; i++){
                int v = g.p_edges[i];
                if(label[v] == 2) ++ vs_deg[e];
            }
            for(int i = g.n_pstart[e]; i < g.n_pend[e]; i++){
                int v = g.n_edges[i];
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
                for(int i = g.p_pstart[next_v]; i < g.p_pend[next_v]; i++){
                    int v = g.p_edges[i];
                    if(label[v] == 1) new_vsP.push_back(v);
                }
                for(int i = g.n_pstart[next_v]; i < g.n_pend[next_v]; i++){
                    int v = g.n_edges[i];
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
                    for(int i = g.p_pstart[e]; i < g.p_pend[e]; i++){
                        int v = g.p_edges[i];
                        if(label[v] == 1) ++ vs_deg[e];
                    }
                    for(int i = g.n_pstart[e]; i < g.n_pend[e]; i++){
                        int v = g.n_edges[i];
                        if(label[v] == 2) ++ vs_deg[e];
                    }
                }
                for(auto e : vsN){
                    for(int i = g.p_pstart[e]; i < g.p_pend[e]; i++){
                        int v = g.p_edges[i];
                        if(label[v] == 2) ++ vs_deg[e];
                    }
                    for(int i = g.n_pstart[e]; i < g.n_pend[e]; i++){
                        int v = g.n_edges[i];
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
                for(int i = g.p_pstart[next_v]; i < g.p_pend[next_v]; i++){
                    int v = g.p_edges[i];
                    if(label[v] == 2) new_vsN.push_back(v);
                }
                for(int i = g.n_pstart[next_v]; i < g.n_pend[next_v]; i++){
                    int v = g.n_edges[i];
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
                    for(int i = g.p_pstart[e]; i < g.p_pend[e]; i++){
                        int v = g.p_edges[i];
                        if(label[v] == 1) ++ vs_deg[e];
                    }
                    for(int i = g.n_pstart[e]; i < g.n_pend[e]; i++){
                        int v = g.n_edges[i];
                        if(label[v] == 2) ++ vs_deg[e];
                    }
                }
                for(auto e : vsN){
                    for(int i = g.p_pstart[e]; i < g.p_pend[e]; i++){
                        int v = g.p_edges[i];
                        if(label[v] == 2) ++ vs_deg[e];
                    }
                    for(int i = g.n_pstart[e]; i < g.n_pend[e]; i++){
                        int v = g.n_edges[i];
                        if(label[v] == 1) ++ vs_deg[e];
                    }
                }
			}
		}

		// cout << "clique_size = " << res.size() << endl;
		if(res.size() > P.size()) {
			P = res;
		}

		ordV.clear();
		for(int i = 0; i < g.n; i++) {
			if(inR[i]) continue;
			ordV.push_back(make_pair(g.degree[i], i));
			sort(ordV.begin(), ordV.end(), greater<pair<int, int>>()); //decreasing order
		}
		for(auto pa : ordV) {
			int v = pa.second;
			// printf("%d ", v);
			res.push_back(v);
			if(!cal(k, res.size(), res)) res.pop_back();
		}

		// cout << "kplex_size = " << res.size() << endl;
		if(res.size() > P.size()) {
			P = res;
		}
	}
	cout << "\t hec_kplex_size = " << P.size() << endl;
	// for(auto u : P) {
	// 	printf("%d ", u);
	// }
	// printf("\n");
}

int main(int argc, const char * argv[])
{
    if(argc < 2) {
        cout<<"\t Usage: [0]exe [1]input_graph [2]k (optional)\t"<<endl; exit(1);
    }

	// load graph
	load_graph(argv[1]);
	int k = 3; //size constraint
    if(argc > 2) k = atoi(argv[2]);	
	cout<<"\t Graph: "<<argv[1]<<",\t k: "<<k<<endl;

	init_hash();
	v_del = new int[G.n];
	memset(v_del, 0, sizeof(int)*G.n);

	// find heuristic signed k-plex
	heu_signed_kplex(1, k);
	ans = P.size();

	if((int)P.size() < ub) {
		lb = max((int)P.size(), 2*k-2);
		get_G_core(lb+1-k);
		get_G_tricnt();
		CTCP(-1, 1, lb+1-k, lb+1-2*k);

		int nn = 0, nm = 0;
		for(int u = 0; u < G.n; u++) {
			if(!v_del[u]) {
				// printf("u = %d\n", u);
				nn++;
				for(int i = G.pstart[u]; i < G.pend[u]; i++) {
					int v = G.edges[i];
					
					if(!v_del[v] && !e_del[u][v] && u < v) {
						// printf("%d %d\n", u, v);
						nm++;
					} 
				}
			}
		}
		printf("%d %d\n", nn, nm);
		// return 0;

		while(G.n > 0) {
			int u = -1;
			for(int i = 0; i < G.n; i++) if(!v_del[i]) {
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
