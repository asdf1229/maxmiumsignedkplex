#include <bits/stdc++.h>
#include "Timer.h"
#include "Utility.h"

int k; //size constraint
int tot, ans;
int Empty = 1;

ui n; //number of vertices
ui m; //number of edges
ui pm; //number of positive edges
ui nm; //number of negative edges

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

ui * bit_del; //mark whether the node is deleted

int ** Matrix;

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
		input_file >> n >> m;
		map<ui, int> * s_G = new map<ui, int>[n];
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
		m = 0; pm = 0; nm = 0;
		for (ui i = 0; i < n; i++) {
			const map<ui, int> & nei = s_G[i];
			for (auto e : nei) {
				if (e.second == 1)
					++pm;
				else
					++nm;
			}
			m += nei.size();
		}
		assert(m%2 == 0);assert(pm%2 == 0);assert(nm%2 == 0);
		m /= 2; pm /= 2; nm /= 2;

		input_file.close();

		p_pstart = new ui[n+1];
        p_edges = new ui[2*pm];
        p_pend = new ui[n];
        n_pstart = new ui[n+1];
        n_edges = new ui[2*nm];
        n_pend = new ui[n];
        
        degree = new ui[n];
        p_degree = new ui[n];
        n_degree = new ui[n];

		//construct positive edges
		p_pstart[0] = 0;
		for (ui i = 0; i < n; i++) {
			const map<ui, int> & nei = s_G[i];
			ui start_idx = p_pstart[i];
			ui d = 0;
			for (auto e : nei) {
				if (e.second == 1) {
					p_edges[start_idx++] = e.first;
					d++;
				}
			}
			p_pstart[i+1] = start_idx;
			p_degree[i] = d;
		}
		assert(p_pstart[n] == 2*pm);

		//construct negative edges
		n_pstart[0] = 0;
		for (ui i = 0; i < n; i++) {
			const map<ui, int> & nei = s_G[i];
			ui start_idx = n_pstart[i];
			ui d = 0;
			for (auto e : nei) {
				if (e.second == -1) {
					n_edges[start_idx++] = e.first;
					d++;
				}
			}
			n_pstart[i+1] = start_idx;
			n_degree[i] = d;
		}
		assert(n_pstart[n] == 2*nm);

		for (ui i = 0; i < n; i++) {
            p_pend[i] = p_pstart[i+1];
            n_pend[i] = n_pstart[i+1];
		}

		for (ui i = 0; i < n; i++) {
			degree[i] = p_degree[i] + n_degree[i];
		}
		delete [] s_G;
	}
	cout<<"\t load_graph: time cost = "<<integer_to_string(t.elapsed())<<endl;
	cout<<"\t G : n = "<<n<<", m = "<<m<<", pm = "<<pm<<", nm = "<<nm<<endl;
}

void ComputeCore(ui thr) {
	if(thr < 2) return; //threshold should be at least 2
	Timer t;
	t.restart();	
	ui threshold = thr - 1;
	ui del_count = 0;
	bit_del = new ui[n];
	memset(bit_del, 0, sizeof(ui)*n);
	queue<ui> q;

	for (ui i = 0; i < n; i++) if(degree[i] < threshold) q.push(i);
	while (!q.empty()) {
		ui u = q.front();
		q.pop();
		del_count++;
		bit_del[u] = 1;
        for (ui i = p_pstart[u]; i < p_pend[u]; i++) {
			ui v = p_edges[i];
			p_degree[v]--;
			if (p_degree[v] + n_degree[v] == threshold) {
				q.push(v);
			}
        }
        for (ui i = n_pstart[u]; i < n_pend[u]; i++) {
			ui v = n_edges[i];
			n_degree[v]--;
			if (p_degree[v] + n_degree[v] == threshold) {
				q.push(v);
			}
        }
	}

	//rebuild
	ui * mapping = new ui[n];
	ui idx = 0;
	for (ui i = 0; i < n; i++) {
		if (!bit_del[i]) {
			mapping[i] = idx++;
		}
	}

    ui * t_p_pstart = new ui[n+1];
    ui * t_p_edges = new ui[2*pm];
    ui * t_n_pstart = new ui[n+1];
    ui * t_n_edges = new ui[2*nm];
	
	ui new_i = 0;
	t_p_pstart[0] = 0;
	for (ui i = 0; i < n; i++) if (!bit_del[i]) {
		ui start_idx = t_p_pstart[new_i];
		for (ui j = p_pstart[i]; j < p_pend[i]; j++) {
			ui v = p_edges[j];
			if (!bit_del[v]) {
				t_p_edges[start_idx++] = mapping[v];
			}
		}
		t_p_pstart[++new_i] = start_idx;
	}
	assert(new_i == idx);
	new_i = 0;
	t_n_pstart[0] = 0;
	for (ui i = 0; i < n; i++) if (!bit_del[i]) {
		ui start_idx = t_n_pstart[new_i];
		for (ui j = n_pstart[i]; j < n_pend[i]; j++) {
			ui v = n_edges[j];
			if (!bit_del[v]) {
				t_n_edges[start_idx++] = mapping[v];
			}
		}
		t_n_pstart[++new_i] = start_idx;
	}
	assert(new_i == idx);
	n = idx;

    delete [] p_pstart;
    delete [] p_edges;
    delete [] n_pstart;
    delete [] n_edges;
    delete [] mapping;

    p_pstart = t_p_pstart;
    p_edges = t_p_edges;
    n_pstart = t_n_pstart;
    n_edges = t_n_edges;

    for (ui i = 0; i < n; i++){
        p_pend[i] = p_pstart[i+1];
        n_pend[i] = n_pstart[i+1];
    }

	//rebuild degree
	for (ui i = 0; i < n; i++) {
        p_degree[i] = p_pstart[i+1] - p_pstart[i];
        n_degree[i] = n_pstart[i+1] - n_pstart[i];
        degree[i] = p_degree[i] + n_degree[i];
	}

    ui now_m = 0, now_pm = 0, now_nm = 0;
    for(ui i = 0; i < n; i++) {now_m += degree[i];now_pm += p_degree[i]; now_nm += n_degree[i];}
    assert(now_m%2 == 0); assert(now_pm%2 == 0); assert(now_nm%2 == 0);
    now_m /= 2; now_pm /= 2; now_nm /= 2;
	cout<<"\t vertex reduce, T : "<<integer_to_string(t.elapsed())<<",\t n="<<n<<", m="<<now_m<<endl;
}

void NepRunning() {

}

void Construct_Matrix() {
	Matrix = new int*[n];
	for(ui i = 0; i < n; i++) Matrix[i] = new int[n];
	for(ui i = 0; i < n; i++)
		for(ui j = 0; j < n; j++)
			Matrix[i][j] = 0;

	for(ui i = 0; i < n; i++) {
		ui u = i;
		for(ui j = p_pstart[i]; j < p_pend[i]; j++) {
			ui v = p_edges[j];
			Matrix[u][v] = 1;
		}
		for(ui j = n_pstart[i]; j < n_pend[i]; j++) {
			ui v = n_edges[j];
			Matrix[u][v] = -1;
		}
	}
}

void EnumSKCE(vector<ui> st, vector<ui> sd, vector<ui> c, vector<ui> ct, vector<ui> cd, vector<ui> x, vector<ui> xt, vector<ui> xd) {
	ui st_size = st.size();
	ui sd_size = sd.size();
	ui c_size = c.size();
	ui ct_size = ct.size();
	ui cd_size = cd.size();
	printf("%d %d %d\n", st_size + sd_size, c_size, (int)x.size());

	// if(st_size + sd_size + c_size < ans) return;
	if(st_size + sd_size + c_size < k) return;
	if(c_size == 0) {
		if(x.size() == 0) {
			tot++;
			// vector<ui> tmp;
			// for(ui i = 0; i < st_size; i++) {
			// 	tmp.push_back(st[i]);
			// 	// printf("%d ", st[i]);
			// }
			// for(ui i = 0; i < sd_size; i++) {
			// 	tmp.push_back(sd[i]);
			// 	// printf("%d ", sd[i]);
			// }
			// sort(tmp.begin(), tmp.end());
			// for(ui i = 0; i < tmp.size(); i++) {
			// 	printf("%d ", tmp[i]);
			// }
			// printf("\n");
			ans = max(ans, (int)(st_size + sd_size));
		}
		return;
	}
	if(Empty == 1) {
		Empty = 0;
		for(ui i = 0; i < c_size; i++) {
			ui u = c[i];
			vector<ui> next_st, next_sd;
			next_st.clear(); next_sd.clear();
			next_st.push_back(u);

			vector<ui> next_c, next_ct, next_cd;
			next_c.clear(); next_ct.clear(); next_cd.clear();
			for(ui j = 0; j < c_size; j++) {
				ui v = c[j];
				if(Matrix[u][v] == 1) {
					next_ct.push_back(v);
					next_c.push_back(v);
				}
				if(Matrix[u][v] == -1) {
					next_cd.push_back(v);
					next_c.push_back(v);
				}
			}

			vector<ui> next_x, next_xt, next_xd;
			for(ui j = 0; j < x.size(); j++) {
				ui v = x[j];
				if(Matrix[u][v] == 1) {
					next_xt.push_back(v);
					next_x.push_back(v);
				}
				if(Matrix[u][v] == -1) {
					next_xd.push_back(v);
					next_x.push_back(v);
				}
			}
			EnumSKCE(next_st, next_sd, next_c, next_ct, next_cd, next_x, next_xt, next_xd);
			x.push_back(u);
		}
	}
	else {
        vector<int> pivot_indicator(c_size);
        vector<int> exp_indicator(c_size, 1);

		for(int i = 0; i < c_size; i++) {
			pivot_indicator[i] = 1;
		}

		for(ui i = 0; i < ct_size; i++) {
			ui u = ct[i];

			if(pivot_indicator[i] == 0) continue;

			vector<ui> next_st(st), next_sd(sd);
			next_st.push_back(u);

			for(ui j = 0; j < ct_size; j++) {
				ui v = ct[j];
				if(Matrix[u][v] != 0) pivot_indicator[j] = 0;
			}
			for(ui j = 0; j < cd_size; j++) {
				ui v = cd[j];
				if(Matrix[u][v] != 0) pivot_indicator[j + ct_size] = 0;
			}

			vector<ui> next_c, next_ct, next_cd;
			for(ui j = 0; j < ct_size; j++) {
				ui v = ct[j];
				if(Matrix[u][v] == 1) {
					next_ct.push_back(v);
					next_c.push_back(v);
				}
			}
			for(ui j = 0; j < cd_size; j++) {
				ui v = cd[j];
				if(Matrix[u][v] == -1) {
					next_cd.push_back(v);
					next_c.push_back(v);
				}
			}

			vector<ui> next_x, next_xt, next_xd;
			for(ui j = 0; j < xt.size(); j++) {
				ui v = xt[j];
				if(Matrix[u][v] == 1) {
					next_xt.push_back(v);
					next_x.push_back(v);
				}
			}
			for(ui j = 0; j < xd.size(); j++) {
				ui v = xd[j];
				if(Matrix[u][v] == -1) {
					next_xd.push_back(v);
					next_x.push_back(v);
				}
			}

			EnumSKCE(next_st, next_sd, next_c, next_ct, next_cd, next_x, next_xt, next_xd);
			xt.push_back(u);
		}
		
		for(ui i = 0; i < cd_size; i++) {
			ui u = cd[i];

			if(pivot_indicator[i] == 0) continue;

			vector<ui> next_st(st), next_sd(sd);
			next_sd.push_back(u);

			for(ui j = 0; j < ct_size; j++) {
				ui v = ct[j];
				if(Matrix[u][v] != 0) pivot_indicator[j] = 0;
			}
			for(ui j = 0; j < cd_size; j++) {
				ui v = cd[j];
				if(Matrix[u][v] != 0) pivot_indicator[j + ct_size] = 0;
			}

			vector<ui> next_c, next_ct, next_cd;
			for(ui j = 0; j < ct_size; j++) {
				ui v = ct[j];
				if(Matrix[u][v] == -1) {
					next_ct.push_back(v);
					next_c.push_back(v);
				}
			}
			for(ui j = 0; j < cd_size; j++) {
				ui v = cd[j];
				if(Matrix[u][v] == 1) {
					next_cd.push_back(v);
					next_c.push_back(v);
				}
			}

			vector<ui> next_x, next_xt, next_xd;
			for(ui j = 0; j < xt.size(); j++) {
				ui v = xt[j];
				if(Matrix[u][v] == -1) {
					next_xt.push_back(v);
					next_x.push_back(v);
				}
			}
			for(ui j = 0; j < xd.size(); j++) {
				ui v = xd[j];
				if(Matrix[u][v] == 1) {
					next_xd.push_back(v);
					next_x.push_back(v);
				}
			}

			EnumSKCE(next_st, next_sd, next_c, next_ct, next_cd, next_x, next_xt, next_xd);
			xd.push_back(u);
		}	
	}
}

void delete_memo()
{
    if(p_pstart != nullptr){
        delete [] p_pstart;
        p_pstart = nullptr;
    }
    if(p_edges != nullptr){
        delete [] p_edges;
        p_edges = nullptr;
    }
    if(p_pend != nullptr){
        delete [] p_pend;
        p_pend = nullptr;
    }
    if(n_pstart != nullptr)
    {
        delete [] n_pstart;
        n_pstart = nullptr;
    }
    if(n_edges != nullptr){
        delete [] n_edges;
        n_edges = nullptr;
    }
    if(n_pend != nullptr){
        delete [] n_pend;
        n_pend = nullptr;
    }
    if(degree != nullptr){
        delete [] degree;
        degree = nullptr;
    }
    if(p_degree != nullptr){
        delete [] p_degree;
        p_degree = nullptr;
    }
    if(n_degree != nullptr){
        delete [] n_degree;
        n_degree = nullptr;
    }
    if(bit_del != nullptr){
        delete [] bit_del;
        bit_del = nullptr;
    }
}

int main(int argc, const char * argv[]) {

    if(argc < 2) {
        cout<<"\t Usage: [0]exe [1]input_graph [2]k (optional)\t"<<endl; exit(1);
    }

	load_graph(argv[1]);

    k = 3;
    if(argc > 2) k = atoi(argv[2]);	

	cout<<"\t Graph: "<<argv[1]<<",\t k: "<<k<<endl;

	ComputeCore(k);

	// NepRunning();

	Construct_Matrix();

	Empty = 1;

	vector<ui> st;
	vector<ui> sd;
	vector<ui> c;
	vector<ui> ct;
	vector<ui> cd;
	vector<ui> x;
	vector<ui> xt;
	vector<ui> xd;

	st.clear();
	sd.clear();
	c.clear();
	ct.clear();
	cd.clear();
	x.clear();
	xt.clear();
	xd.clear();
	
	for(int i = 0; i < n; i++) {
		c.push_back(i);
	}
	Timer t;
	EnumSKCE(st, sd, c, ct, cd, x, xt, xd);
	cout<<"\t EnumSKCE: time cost  = "<<integer_to_string(t.elapsed())<<endl;
	printf("tot = %d\n", tot);
	printf("ans = %d\n", ans);

	delete_memo();

	return 0;
}
