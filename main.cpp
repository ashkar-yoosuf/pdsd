#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

#define INITIAL_CAPACITY 2

struct Graph
{
	int V;
	int** arrayOfArrays;
	int* degrees;
	int* capacities;
};

void push(int **arr, int index, int value, int *size, int *capacity)
{
	if((*size + 1) > *capacity){
		*arr = (int*)realloc(*arr, ((*size + 1) * sizeof(int)));
		*capacity = (*size + 1);
	}
	(*arr)[index] = value;
	*size = *size + 1;
}

void addEdge(struct Graph* graph, int src, int dest)
{
	push(&(graph->arrayOfArrays[src]), graph->degrees[src], dest, &(graph->degrees[src]), &(graph->capacities[src]));

	push(&(graph->arrayOfArrays[dest]), graph->degrees[dest], src, &(graph->degrees[dest]), &(graph->capacities[dest]));

}

struct Graph* createGraph(int V)
{
	struct Graph* graph =
			(struct Graph*) malloc(sizeof(struct Graph));

	graph->V = V;

	graph->arrayOfArrays =
			(int **) malloc(V * sizeof(*(graph->arrayOfArrays)));

	graph->degrees =
			(int *) malloc(V * sizeof(*(graph->degrees)));

	graph->capacities =
			(int *) malloc(V * sizeof(*(graph->capacities)));

	int i;
	for (i = 0; i < V; ++i) {
		graph->arrayOfArrays[i] = (int *)malloc(INITIAL_CAPACITY * sizeof(int));
		graph->degrees[i] = 0;
		graph->capacities[i] = INITIAL_CAPACITY;
	}
	return graph;
}

void printGraph(struct Graph* graph)
{
	int v;
	int i;
	struct Graph* p_graph = graph;
	for (v = 0; v < p_graph->V; ++v)
	{
		if (p_graph->degrees[v] != 0)
		{
			cout << "[" << v << "] --> [";
			for (i = 0; i < p_graph->degrees[v]; ++i)
			{
				if (i != p_graph->degrees[v] - 1)
					cout << p_graph->arrayOfArrays[v][i] << ",";
				else
					cout << p_graph->arrayOfArrays[v][i] << "]" << endl;
			}
		}
	}
}

int main()
{
	int V = 27;
	int E = 0; // yet to update
	struct Graph* graph = createGraph(V);
	ifstream ip("./input/new.csv");

	string node_1;
	string node_2;

	int src = 0;
	int dest = 0;

	int line = 0;

	while (ip.good()) {
		getline(ip, node_1, ',');
		getline(ip, node_2, '\n');

		stringstream start_node(node_1), end_node(node_2);

		start_node >> src;
		end_node >> dest;

		if (line < 2)
			line++;

		if (line == 2 && !node_1.empty()) {
			addEdge(graph, src, dest);
			E++;
		}
	}

	printGraph(graph);
	cout << "#Edges: " << E << endl;

	for (int i = 0; i < V; ++i)
	{
		free(graph->arrayOfArrays[i]);
	}
	free (graph->degrees);
	free (graph->capacities);
	free (graph->arrayOfArrays);
	free (graph);

	ip.close();
	exit(0);
}
