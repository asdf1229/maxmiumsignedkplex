#include <bits/stdc++.h>
#include "Timer.h"
#include "Utility.h"
using namespace std;

int k; //size constraint
int tot, ans;
int Empty;

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

void ComputeCore() {
	if(k < 1) {cout<<"\t k should be at least 1"<<endl; exit(1); }
	Timer t;
	t.restart();	
	ui threshold = k - 1;
	ui del_count = 0;
	bit_del = new ui[n];
	memset(bit_del, 0, sizeof(ui)*n);

    // int pt = k - 1;
    // int nt = k;
	// queue<ui> Q;

    // for(ui i = 0; i < n; i++) if(p_degree[i] < pt || n_degree[i] < nt) Q.push(i);
    // while (!Q.empty()) {
    //     ui u = Q.front();
    //     Q.pop();
    //     ++ del_count;
    //     bit_del[u] = 1;
    //     for(ui i = p_pstart[u]; i < p_pend[u]; i++){
    //         if((p_degree[p_edges[i]]--) == pt && n_degree[p_edges[i]] >= nt){
    //             Q.push(p_edges[i]);
    //         }
    //     }
    //     for(ui i = n_pstart[u]; i < n_pend[u]; i++){
    //         if((n_degree[n_edges[i]]--) == nt && p_degree[n_edges[i]] >= pt){
    //             Q.push(n_edges[i]);
    //         }
    //     }
    // }

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
	cout<<"\t vertex reduce: time cost = "<<integer_to_string(t.elapsed())<<",\t n="<<n<<", m="<<now_m<<endl;
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

void EnumBase(vector<ui> s, vector<ui> ct, vector<ui> cd) {
	ui c_size = ct.size() + cd.size();
	// if(s.size() + c_size < ans) return;
	if(s.size() + c_size < k) return;
	if(c_size == 0) {
		tot++;
		ans = max(ans, (int)(s.size()));
		return;
	}

	ui ct_size = ct.size();
	ui cd_size = cd.size();

	if(Empty == 1) {
		Empty = 0;
		for(ui i = 0; i < ct_size; i++) {
			ui u = ct[i];
			vector<ui> next_s;
			next_s.push_back(u);

			vector<ui> next_ct, next_cd;

			for(ui j = 0; j < ct_size; j++) {
				ui v = ct[j];
				if(Matrix[u][v] == 1) next_ct.push_back(v);
				if(Matrix[u][v] == -1) next_cd.push_back(v);
			}

			EnumBase(next_s, next_ct, next_cd);
		}
	}
	else {
		vector<ui> pivot_indicator(c_size);
		vector<ui> exp_indicator(c_size);

		for(ui i = 0; i < c_size; i++) {
			pivot_indicator[i] = 1;
			exp_indicator[i] = 1;
		}

		for(ui i = 0; i < ct_size; i++) {
			ui u = ct[i];
			if(pivot_indicator[i] == 0) continue;

			vector<ui> next_s(s);
			next_s.push_back(u);

			for(ui j = 0; j < ct_size; j++) {
				ui v = ct[j];
				if(Matrix[u][v] != 0) pivot_indicator[j] = 0;
			}
			for(ui j = 0; j < cd_size; j++) {
				ui v = cd[j];
				if(Matrix[u][v] != 0) pivot_indicator[j + ct_size] = 0;
			}

			vector<ui> next_ct, next_cd;

			for(ui j = 0; j < ct_size; j++) {
				ui v = ct[j];
				if(exp_indicator[j] == 1) {
					if(Matrix[u][v] == 1) next_ct.push_back(v);
				}
			}
			for(ui j = 0; j < cd_size; j++) {
				ui v = cd[j];
				if(exp_indicator[j + ct_size] == 1) {
					if(Matrix[u][v] == -1) next_cd.push_back(v);
				}
			}

			EnumBase(next_s, next_ct, next_cd);
			exp_indicator[i] = 0;
		}

		for(ui i = 0; i < cd_size; i++) {
			ui u = cd[i];
			if(pivot_indicator[i + ct_size] == 0) continue;

			vector<ui> next_s(s);
			next_s.push_back(u);

			for(ui j = 0; j < ct_size; j++) {
				ui v = ct[j];
				if(Matrix[u][v] != 0) pivot_indicator[j] = 0;
			}
			for(ui j = 0; j < cd_size; j++) {
				ui v = cd[j];
				if(Matrix[u][v] != 0) pivot_indicator[j + ct_size] = 0;
			}

			vector<ui> next_ct, next_cd;

			for(ui j = 0; j < ct_size; j++) {
				ui v = ct[j];
				if(exp_indicator[j] == 1) {
					if(Matrix[u][v] == -1) next_ct.push_back(v);
				}
			}
			for(ui j = 0; j < cd_size; j++) {
				ui v = cd[j];
				if(exp_indicator[j + ct_size] == 1) {
					if(Matrix[u][v] == 1) next_cd.push_back(v);
				}
			}

			EnumBase(next_s, next_ct, next_cd);
			exp_indicator[i] = 0;
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

	ComputeCore();

	Construct_Matrix();

	printf("\t n = %d\n", n);

	Empty = 1;

	vector<ui> s;
	vector<ui> ct;
	vector<ui> cd;

	s.clear();
	ct.clear();
	cd.clear();

	for(int i = 0; i < n; i++) {
		ct.push_back(i);
	}

	Timer t;
	EnumBase(s, ct, cd);
	cout<<"\t EnumBase: time cost = "<<integer_to_string(t.elapsed())<<endl;
	printf("\t tot = %d\n", tot);
	printf("\t ans = %d\n", ans);

	delete_memo();
    return 0;
}
