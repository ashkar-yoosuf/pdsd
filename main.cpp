#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>

using namespace std;

#define NO_OF_VERTICES 4039
#define EPSILON 0.5
#define approxFactor(EPSILON) (2 + 2 * EPSILON)

struct Graph
{
	int V;
	int E;
	int** arrayOfArrays;
	int* degrees;
	int* capacities;
};

struct GraphTilde
{
	int V{};
	int E{};
	int V_arr[NO_OF_VERTICES] = {};
};

void push(int **arr, int index, int value, int *size, int *capacity)
{
	if (*capacity != 0) {
		*arr = (int *) realloc(*arr, ((*size + 1) * sizeof(int)));
	} else {
		*arr = (int *) malloc(sizeof(int));
	}
	*capacity = (*size + 1);
	(*arr)[index] = value;
	*size = *size + 1;
}

void addEdge(struct Graph* graph, int src, int dest)
{
	push(&(graph->arrayOfArrays[src]), graph->degrees[src], dest, &(graph->degrees[src]), &(graph->capacities[src]));

	push(&(graph->arrayOfArrays[dest]), graph->degrees[dest], src, &(graph->degrees[dest]), &(graph->capacities[dest]));

}

struct Graph* createGraph()
{
	auto* graph =
			(struct Graph*) malloc(sizeof(struct Graph));

	graph->V = NO_OF_VERTICES;

	graph->E = 0;

	graph->arrayOfArrays =
			(int **) malloc(NO_OF_VERTICES * sizeof(*(graph->arrayOfArrays)));

	graph->degrees =
			(int *) malloc(NO_OF_VERTICES * sizeof(*(graph->degrees)));

	graph->capacities =
			(int *) malloc(NO_OF_VERTICES * sizeof(*(graph->capacities)));

	int i;
	for (i = 0; i < NO_OF_VERTICES; ++i) {
		graph->degrees[i] = 0;
		graph->capacities[i] = 0;
	}
	return graph;
}

void printGraph(struct Graph* graph)
{
	struct Graph* p_graph = graph;
	for (int v = 0; v < NO_OF_VERTICES; ++v)
	{
		if (p_graph->capacities[v] != 0)
		{
			cout << "[" << v << "] --> [";
			for (int i = 0; i < p_graph->capacities[v]; ++i)
			{
				if (i != p_graph->capacities[v] - 1)
					cout << p_graph->arrayOfArrays[v][i] << ",";
				else
					cout << p_graph->arrayOfArrays[v][i] << "] | degree: " << p_graph->degrees[v] << endl;
			}
		}
	}
}

inline void calcDensity(int num_edges, int num_vertices, double* density)
{
		*density = (double) num_edges/  (double) num_vertices;
}

inline void assignToTilde(struct Graph* graph, struct GraphTilde* graph_tilde)
{
	graph_tilde->V = graph->V;
	graph_tilde->E = graph->E;
}

void deallocateGraph(struct Graph* graph)
{
	for (int i = 0; i < NO_OF_VERTICES; ++i)
	{
		free(graph->arrayOfArrays[i]);
	}
	free (graph->degrees);
	free (graph->capacities);
	free (graph->arrayOfArrays);
	free (graph);
}

void densestComponent(struct Graph* graph, double rho_init)
{
	auto* graph_tilde = (struct GraphTilde *) malloc(sizeof(struct GraphTilde));
	graph_tilde->V = NO_OF_VERTICES;
	graph_tilde->E = graph->E;

	int target_element;
	double current_graph_rho = rho_init, current_graphTilde_rho = rho_init;
	bool isTildeChanged = false;

	while (graph->V > 0) {
		for (int v = 0; v < NO_OF_VERTICES; v++) {
			if ((graph_tilde->V_arr[v] != -1) && (graph->degrees[v] <= approxFactor(EPSILON) * current_graph_rho)) {
				graph_tilde->V_arr[v] = -1;
				graph->V--;
				for (int v_v = 0; v_v < graph->capacities[v]; v_v++) {
					target_element = graph->arrayOfArrays[v][v_v];
					if (target_element != -1) {
						for (int u_u = 0; u_u < graph->capacities[target_element]; u_u++) {
							if (graph->arrayOfArrays[target_element][u_u] == v) {
								graph->arrayOfArrays[target_element][u_u] = -1;
								graph->E--;
								break;
							}
						}
						graph->degrees[target_element]--;
						graph->arrayOfArrays[v][v_v] = -1;
					}
				}
				graph->degrees[v] = 0;
			}
		}
		calcDensity(graph->E, graph->V, &current_graph_rho);

		if (isTildeChanged) {
			calcDensity(graph_tilde->E, graph_tilde->V, &current_graphTilde_rho);
			isTildeChanged = false;
		}

		if (current_graph_rho > current_graphTilde_rho) {
			assignToTilde(graph, graph_tilde);
			isTildeChanged = true;
		}
	}

	cout << "Density of the densest component: " << current_graphTilde_rho << endl;

	free(graph_tilde);
}

int main()
{
	double rho_init;
	struct Graph* graph = createGraph();
	ifstream ip("./input/facebook_combined.csv");

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
			graph->E++;
		}
	}

	calcDensity(graph->E, NO_OF_VERTICES, &rho_init);

//	printGraph(graph);

	cout << "#Edges:" << graph->E << endl;
	cout << "Initial Density:" << rho_init << endl;

	auto start = chrono::steady_clock::now();
	densestComponent(graph, rho_init);
	auto end = chrono::steady_clock::now();

	cout << "Elapsed time in nanoseconds : "
			 << chrono::duration_cast<chrono::nanoseconds>(end - start).count()
			 << " ns" << endl;

	deallocateGraph(graph);

	ip.close();
	exit(0);
}
