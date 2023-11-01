#include <bits/stdc++.h>
#include "Timer.h"
#include "Utility.h"
using namespace std;

vector<ui> P;
vector<ui> newP;
int lb, ub;

struct Graph {
	int n; //number of vertices
	int m; //number of edges
	int pm; //number of positive edges
	int nm; //number of negative edges

	/*Store edges in linear arrays*/
	ui * p_pstart; //start of positive edge number of a point
	ui * p_pend; //end of positive edge number of a point
	ui * p_edges; //positive edges

	ui * n_pstart; //start of negative edge number of a point
	ui * n_pend; //end of negative edge number of a point
	ui * n_edges; //negative edges

	ui * degree; //degree of a point
	ui * p_degree; //positive degree of a point
	ui * n_degree; //negative degree of a point

	int ** Matrix;

	//Suitable for situations with small data scales
	int ** Tricnt; 
}G, g;

ui * bit_del; //mark whether the node is deleted
ui * bit_sel; //mark whether the node is selected

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
		map<ui, int> * s_G = new map<ui, int>[G.n];
		ui u, v;
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
		for (ui i = 0; i < G.n; i++) {
			const map<ui, int> & nei = s_G[i];
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

		G.p_pstart = new ui[G.n+1];
        G.p_edges = new ui[2*G.pm];
        G.p_pend = new ui[G.n];
        G.n_pstart = new ui[G.n+1];
        G.n_edges = new ui[2*G.nm];
        G.n_pend = new ui[G.n];
        
        G.degree = new ui[G.n];
        G.p_degree = new ui[G.n];
        G.n_degree = new ui[G.n];

		//construct positive edges
		G.p_pstart[0] = 0;
		for (ui i = 0; i < G.n; i++) {
			const map<ui, int> & nei = s_G[i];
			ui start_idx = G.p_pstart[i];
			ui d = 0;
			for (auto e : nei) {
				if (e.second == 1) {
					G.p_edges[start_idx++] = e.first;
					d++;
				}
			}
			G.p_pstart[i+1] = start_idx;
			G.p_degree[i] = d;
		}
		assert(p_pstart[n] == 2*pm);

		//construct negative edges
		G.n_pstart[0] = 0;
		for (ui i = 0; i < G.n; i++) {
			const map<ui, int> & nei = s_G[i];
			ui start_idx = G.n_pstart[i];
			ui d = 0;
			for (auto e : nei) {
				if (e.second == -1) {
					G.n_edges[start_idx++] = e.first;
					d++;
				}
			}
			G.n_pstart[i+1] = start_idx;
			G.n_degree[i] = d;
		}
		assert(n_pstart[n] == 2*nm);

		for (ui i = 0; i < G.n; i++) {
            G.p_pend[i] = G.p_pstart[i+1];
            G.n_pend[i] = G.n_pstart[i+1];
		}

		for (ui i = 0; i < G.n; i++) {
			G.degree[i] = G.p_degree[i] + G.n_degree[i];
		}
		delete [] s_G;
	}
	cout<<"\t load_graph: time cost = "<<integer_to_string(t.elapsed())<<endl;
	cout<<"\t G : n = "<<G.n<<", m = "<<G.m<<", pm = "<<G.pm<<", nm = "<<G.nm<<endl;
}

//get G's k-core
void get_core(int k) {
	if(k < 2) return; //threshold should be at least 2
	Timer t;
	t.restart();	
	ui threshold = k - 1;
	ui del_count = 0;
	bit_del = new ui[G.n];
	memset(bit_del, 0, sizeof(ui)*G.n);

	queue<ui> q;

	for (ui i = 0; i < G.n; i++) if(G.degree[i] < threshold) q.push(i);
	while (!q.empty()) {
		ui u = q.front();
		q.pop();
		del_count++;
		bit_del[u] = 1;
        for (ui i = G.p_pstart[u]; i < G.p_pend[u]; i++) {
			ui v = G.p_edges[i];
			G.p_degree[v]--;
			if (G.p_degree[v] + G.n_degree[v] == threshold) {
				q.push(v);
			}
        }
        for (ui i = G.n_pstart[u]; i < G.n_pend[u]; i++) {
			ui v = G.n_edges[i];
			G.n_degree[v]--;
			if (G.p_degree[v] + G.n_degree[v] == threshold) {
				q.push(v);
			}
        }
	}

	//rebuild
	ui * mapping = new ui[G.n];
	ui idx = 0;
	for (ui i = 0; i < G.n; i++) {
		if (!bit_del[i]) {
			mapping[i] = idx++;
		}
	}

    ui * t_p_pstart = new ui[G.n+1];
    ui * t_p_edges = new ui[2*G.pm];
    ui * t_n_pstart = new ui[G.n+1];
    ui * t_n_edges = new ui[2*G.nm];
	
	ui new_i = 0;
	t_p_pstart[0] = 0;
	for (ui i = 0; i < G.n; i++) if (!bit_del[i]) {
		ui start_idx = t_p_pstart[new_i];
		for (ui j = G.p_pstart[i]; j < G.p_pend[i]; j++) {
			ui v = G.p_edges[j];
			if (!bit_del[v]) {
				t_p_edges[start_idx++] = mapping[v];
			}
		}
		t_p_pstart[++new_i] = start_idx;
	}
	assert(new_i == idx);
	new_i = 0;
	t_n_pstart[0] = 0;
	for (ui i = 0; i < G.n; i++) if (!bit_del[i]) {
		ui start_idx = t_n_pstart[new_i];
		for (ui j = G.n_pstart[i]; j < G.n_pend[i]; j++) {
			ui v = G.n_edges[j];
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
	delete bit_del;

    G.p_pstart = t_p_pstart;
    G.p_edges = t_p_edges;
    G.n_pstart = t_n_pstart;
    G.n_edges = t_n_edges;

    for (ui i = 0; i < G.n; i++){
        G.p_pend[i] = G.p_pstart[i+1];
        G.n_pend[i] = G.n_pstart[i+1];
    }

	// //rebuild degree
	// for (ui i = 0; i < G.n; i++) {
    //     G.p_degree[i] = G.p_pstart[i+1] - G.p_pstart[i];
    //     G.n_degree[i] = G.n_pstart[i+1] - G.n_pstart[i];
    //     G.degree[i] = G.p_degree[i] + G.n_degree[i];
	// }

    // ui now_m = 0, now_pm = 0, now_nm = 0;
    // for(ui i = 0; i < G.n; i++) {now_m += G.degree[i];now_pm += G.p_degree[i]; now_nm += G.n_degree[i];}
    // assert(now_m%2 == 0); assert(now_pm%2 == 0); assert(now_nm%2 == 0);
    // now_m /= 2; now_pm /= 2; now_nm /= 2;
	// cout<<"\t vertex reduce, T : "<<integer_to_string(t.elapsed())<<",\t n="<<G.n<<", m="<<now_m<<endl;
}

//get degree for all nodes of G
void get_deg() {
	for (ui i = 0; i < G.n; i++) {
        G.p_degree[i] = G.p_pstart[i+1] - G.p_pstart[i];
        G.n_degree[i] = G.n_pstart[i+1] - G.n_pstart[i];
        G.degree[i] = G.p_degree[i] + G.n_degree[i];
	}
}

//get triangle count for all edges of G
void get_tricnt() {

}

void get_g(ui u) {
	bit_sel = new ui[G.n];
	memset(bit_sel, 0, sizeof(ui)*G.n);
    bit_sel[u] = 1;

	for(ui i = G.p_pstart[u]; i < G.p_pend[u]; i++) {
		ui v1 = G.p_edges[i];
		bit_sel[v1] = 1;
		for(ui j = G.p_pstart[v1]; j < G.p_pend[v1]; j++) {
			ui v2 = G.p_edges[j];
			bit_sel[v2] = 1;
		}
		for(ui j = G.n_pstart[v1]; j < G.n_pend[v1]; j++) {
			ui v2 = G.n_edges[j];
			bit_sel[v2] = 1;
		}
	}
	for(ui i = G.n_pstart[u]; i < G.n_pend[u]; i++) {
		ui v1 = G.n_edges[i];
		bit_sel[v1] = 1;
		for(ui j = G.p_pstart[v1]; j < G.p_pend[v1]; j++) {
			ui v2 = G.p_edges[j];
			bit_sel[v2] = 1;
		}
		for(ui j = G.n_pstart[v1]; j < G.n_pend[v1]; j++) {
			ui v2 = G.n_edges[j];
			bit_sel[v2] = 1;
		}
	}

	//build g
	ui * mapping = new ui[G.n];
	ui idx = 0;
	for (ui i = 0; i < G.n; i++) {
		if (bit_sel[i]) {
			mapping[i] = idx++;
		}
	}

    ui * t_p_pstart = new ui[G.n+1];
    ui * t_p_edges = new ui[2*G.pm];
    ui * t_n_pstart = new ui[G.n+1];
    ui * t_n_edges = new ui[2*G.nm];
	
	ui new_i = 0;
	t_p_pstart[0] = 0;
	for (ui i = 0; i < G.n; i++) if (bit_sel[i]) {
		ui start_idx = t_p_pstart[new_i];
		for (ui j = G.p_pstart[i]; j < G.p_pend[i]; j++) {
			ui v = G.p_edges[j];
			if (bit_sel[v]) {
				t_p_edges[start_idx++] = mapping[v];
			}
		}
		t_p_pstart[++new_i] = start_idx;
	}
	assert(new_i == idx);
	new_i = 0;
	t_n_pstart[0] = 0;
	for (ui i = 0; i < G.n; i++) if (bit_sel[i]) {
		ui start_idx = t_n_pstart[new_i];
		for (ui j = G.n_pstart[i]; j < G.n_pend[i]; j++) {
			ui v = G.n_edges[j];
			if (bit_sel[v]) {
				t_n_edges[start_idx++] = mapping[v];
			}
		}
		t_n_pstart[++new_i] = start_idx;
	}
	assert(new_i == idx);
	g.n = idx;

    delete [] g.p_pstart;
    delete [] g.p_edges;
    delete [] g.n_pstart;
    delete [] g.n_edges;
    delete [] mapping;
	delete bit_sel;

    g.p_pstart = t_p_pstart;
    g.p_edges = t_p_edges;
    g.n_pstart = t_n_pstart;
    g.n_edges = t_n_edges;

    for (ui i = 0; i < g.n; i++){
        g.p_pend[i] = g.p_pstart[i+1];
        g.n_pend[i] = g.n_pstart[i+1];
    }

	//build degree
	for (ui i = 0; i < g.n; i++) {
        g.p_degree[i] = g.p_pstart[i+1] - g.p_pstart[i];
        g.n_degree[i] = g.n_pstart[i+1] - g.n_pstart[i];
        g.degree[i] = g.p_degree[i] + g.n_degree[i];
	}
}

//core-truss co-pruning
void CTCP(int lb_changed, int tv, int te) {

}

bool dfs() {
    for(int i = 0; i < g.n; i++) {
        for(int j = i + 1; j < g.n; j++) {
            for(int k = j + 1; k < g.n; k++) {
                if(Matrix[i][j] && Matrix[j][k] && Matrix[i][k]) {
                    int cn = Matrix[i][j] + Matrix[j][k] + Matrix[i][k];
                    if(cn == 1 || cn == -3) return false;
                }
            }
        }
    }
    return true;
}

int cal(ui k, ui siz, vector<ui> s) {
    vector<int> deg(siz);
    for(int i = 0; i < siz; i++) {
        deg[i] = 0;
    }
    for(int i = 0; i < siz; i++) {
        for(int j = i + 1; j < siz; j++) {
            ui u = s[i], v = s[j];
            if(Matrix[u][v]) {
                deg[i]++;
                deg[j]++;
            }
        }
    }
    for(int i = 0; i < siz; i++) {
        if(deg[i] < siz - k) return 0;
    }

    
    if(dfs()) return 1;
    return 0;
}

void kplex(ui dep, ui siz, vector<ui> s) {
    if(dep == g.n) {
        if(cal(k, siz, s)) {
            newans = max(newans, siz);
        }
        return;
    }
    for(ui i = dep; i < g.n; i++) {
        ui u = i;
        //reduce

        kplex(i + 1, siz, s);

        s.push_back(u);
        kplex(i + 1, siz + 1, s);
        s.pop_back();
    }

}

int main(int argc, const char * argv[]) {
	
    if(argc < 2) {
        cout<<"\t Usage: [0]exe [1]input_graph [2]k (optional)\t"<<endl; exit(1);
    }

	load_graph(argv[1]);

	int k; //size constraint
    k = 2;
    if(argc > 2) k = atoi(argv[2]);	

	cout<<"\t Graph: "<<argv[1]<<",\t k: "<<k<<endl;

	//(P*, ub) <- kplex_degen(g, k)
	//P*, ub

	if((int)P.size() < ub) {
		lb = max((int)P.size(), 2*k-2);
		get_core(lb+1-k);
		get_deg();
		get_tricnt();
		CTCP(1, lb+1-k, lb+1-2*k);

		while(G.n > 0) {
			//get u
			ui u = 0;
			for(int i = 0; i < G.n; i++) {
				if(G.degree[u] > G.degree[i]) {
					u = i;
				}
			}
			//kplex
			newP.clear();

			vector<ui> s;
			kplex(0, 0, s);
			int lb_changed = 0;
			if(newP.size() > lb) {
				P = newP;
				lb = newP.size();
				lb_changed = 1;
			}
			//CTCP
			CTCP(lb_changed, lb+1-k, lb+1-2*k);
			}
	}

	printf("lb = %d\n", lb);
	// delete_memo();
    return 0;
}
