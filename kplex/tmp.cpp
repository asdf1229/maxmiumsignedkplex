// degeneracy-based k-plex
// return an upper bound of the maximum k-plex size
// return dOrder
void kplex_degen(ListLinearHeap *heap, int k, int *dOrder)
{
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

//get degree for all nodes of G
void get_G_deg()
{
	// for (int i = 0; i < G.n; i++) {
    //     G.p_degree[i] = G.p_pstart[i+1] - G.p_pstart[i];
    //     G.n_degree[i] = G.n_pstart[i+1] - G.n_pstart[i];
    //     G.degree[i] = G.p_degree[i] + G.n_degree[i];
	// }
}