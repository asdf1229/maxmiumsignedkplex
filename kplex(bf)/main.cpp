#include "Graph.h"
// #include "Utility.h"

using namespace std;

int main(int argc, const char *argv[])
{
    if(argc < 2) {
        cout<<"\t Usage: [0]exe [1]input_graph [2]k (optional)\t"<<endl; exit(1);
    }
    // load graph
    int k = 3; //size constraint
    if(argc > 2) k = atoi(argv[2]);

    Graph *graph = new Graph(k);
    graph->load_graph(argv[1]);

	cout<<"\t Graph: "<<argv[1]<<",\t k: "<<k<<endl;

    graph->find_signed_kplex();

    // graph->output_one_kplex();

    delete graph;
    return 0;
}