#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>

using namespace std;

#define NO_OF_VERTICES 4039
#define EPSILON 0.1
#define approxFactor(EPSILON) (2 + 2 * EPSILON)

struct Graph
{
	int V;
	int E;
	int flag;
	int** flags;
	int** arrayOfArrays;
	int* degrees;
	int* capacities;
	int* failed_V_arr;
};

struct GraphTilde
{
	int V;
	int E;
};

inline void edgeMemAlloc(int **arr, int index, int value, int *degree, int *capacity) {

	if (*capacity != 0) {
		*arr = (int *) realloc(*arr, ((*capacity + 1) * sizeof(int)));
	} else {
		*arr = (int *) malloc(sizeof(int));
	}

	*capacity = (*capacity + 1);
	(*arr)[index] = value;
	*degree = *degree + 1;

}

inline void push(int **arr, int index, int value, int *degree, int *capacity, int same_node, bool first) {

	if (same_node) {

		if (first) {
			edgeMemAlloc(arr, index, value, degree, capacity);
		} else {
			*degree = *degree + 1;
		}

	} else {
		edgeMemAlloc(arr, index, value, degree, capacity);
	}
}

void addEdge(struct Graph* graph, int src, int dest) {

	bool same_node { false };

	if (src == dest)
		same_node = true;

	push(&(graph->arrayOfArrays[src]), graph->capacities[src], dest, &(graph->degrees[src]),
			&(graph->capacities[src]), same_node, true);
	push(&(graph->arrayOfArrays[dest]), graph->capacities[dest], src, &(graph->degrees[dest]),
			&(graph->capacities[dest]), same_node, false);

}

struct Graph* createGraph() {

	auto* graph =
			(struct Graph*) malloc(sizeof(struct Graph));

	graph->V = NO_OF_VERTICES;

	graph->E = 0;

	graph->flag = -1;

	graph->arrayOfArrays =
			(int **) malloc(NO_OF_VERTICES * sizeof(*(graph->arrayOfArrays)));

	graph->flags =
			(int **) malloc(NO_OF_VERTICES * sizeof(*(graph->flags)));

	graph->degrees =
			(int *) malloc(NO_OF_VERTICES * sizeof(*(graph->degrees)));

	graph->capacities =
			(int *) malloc(NO_OF_VERTICES * sizeof(*(graph->capacities)));


	graph->failed_V_arr =
			(int *) malloc(NO_OF_VERTICES* sizeof(int));

	int i;
	for (i = 0; i < NO_OF_VERTICES; ++i) {
		graph->degrees[i] = 0;
		graph->capacities[i] = 0;
		graph->failed_V_arr[i] = 0;
	}
	return graph;

}

struct GraphTilde* createGraphTilde(struct Graph* graph) {

	auto* graph_tilde =
			(struct GraphTilde*) malloc(sizeof(struct GraphTilde));

	graph_tilde->V = NO_OF_VERTICES;

	graph_tilde->E = graph->E;

	return graph_tilde;

}

void printGraph(struct Graph* graph) {

	for (int v = 0; v < NO_OF_VERTICES; ++v)
	{
		if (graph->capacities[v] != 0)
		{
			cout << "[" << v << "] --> [";
			for (int i = 0; i < graph->capacities[v]; ++i)
			{
				if (i != graph->capacities[v] - 1)
					cout << graph->arrayOfArrays[v][i] << ",";
				else
					cout << graph->arrayOfArrays[v][i] << "] | degree: " << graph->degrees[v] << endl;
			}
		}
	}
}

inline void calcDensity(int num_edges, int num_vertices, float* density) {

	*density = (double) num_edges/  (double) num_vertices;

}

inline void assignToTilde(struct Graph* graph, struct GraphTilde* graph_tilde) {

	graph_tilde->V = graph->V;
	graph_tilde->E = graph->E;

}

void deallocateGraph(struct Graph* graph) {

	for (int i = 0; i < NO_OF_VERTICES; ++i)
	{
		free(graph->arrayOfArrays[i]);
	}
	free (graph->degrees);
	free (graph->capacities);
	free (graph->arrayOfArrays);
	free (graph->flags);
	free (graph->failed_V_arr);
	free (graph);

}

inline void deallocateGraphTilde(struct GraphTilde* graph_tilde) {

	free (graph_tilde);

}

void init_vertexFlags(struct Graph* graph) {

	for (int v = 0; v < NO_OF_VERTICES; v++) {
		graph->flags[v] = new int[graph->capacities[v]];

		for (int i = 0; i < graph->capacities[v]; i++) {
			graph->flags[v][i] = 0;
		}
	}
}

inline void removeEdge(struct Graph* graph, int degree_reduction, int target_element, int flagOf_target_element, int v, int v_v) {

	if (flagOf_target_element == 0) {
		for (int u_u = 0; u_u < graph->capacities[target_element]; u_u++) {
			if (graph->arrayOfArrays[target_element][u_u] == v) {
				graph->flags[target_element][u_u] = graph->flag;
				graph->E--;
				break;
			}
		}
		graph->degrees[target_element] -= degree_reduction;
		graph->flags[v][v_v] = graph->flag;
	}

}

float maxDensity(struct Graph* graph, double rho_init) {

	struct GraphTilde* graph_tilde = createGraphTilde(graph);

	int target_element, flagOf_target_element;
	float current_graph_rho = rho_init, current_graphTilde_rho = rho_init;
	bool isTildeChanged {false};

	while (graph->V > 0) {
		for (int v = 0; v < NO_OF_VERTICES; v++) {
			if ((graph->failed_V_arr[v] == 0) && (graph->degrees[v] <= approxFactor(EPSILON) * current_graph_rho)) {
				graph->failed_V_arr[v] = graph->flag;
				graph->V--;

				for (int v_v = 0; v_v < graph->capacities[v]; v_v++) {

					target_element = graph->arrayOfArrays[v][v_v];
					flagOf_target_element = graph->flags[v][v_v];

					if (target_element != v) {
						removeEdge(graph, 1, target_element, flagOf_target_element, v, v_v);
					} else {
						removeEdge(graph, 2, target_element, flagOf_target_element, v, v_v);
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
			graph->flag--;
			isTildeChanged = true;
		}
	}

	deallocateGraphTilde(graph_tilde);

	return current_graphTilde_rho;

}

void densestComponent(struct Graph* graph) {

	for (int v = 0; v < NO_OF_VERTICES; v++) {
		if (graph->failed_V_arr[v] == graph->flag) {
			cout << "[" << v << "] --> |";
			for (int i = 0; i < graph->capacities[v]; ++i) {
				if (graph->flags[v][i] == graph->flag) {
					cout << graph->arrayOfArrays[v][i] << "|";
				}
			}
			cout << endl;
		}
	}

}

int densestComponentVertices(struct Graph* graph) {

	int dense_nodes = 0;
	for (int v = 0; v < NO_OF_VERTICES; v++) {
		if (graph->failed_V_arr[v] == graph->flag) {
			dense_nodes++;
			cout << v << endl;
		}
	}

	return dense_nodes;

}

int main() {

	freopen("./output/facebook_combined(0.1).txt", "w", stdout);

	double avg_elapsed_time = 0;
	float rho_init, max_density = 0;
	int dense_nodes = 0;

	for (int turn = 0; turn < 5; turn++) {

		ifstream ip("./input/facebook_combined.csv");

		struct Graph* graph = createGraph();

		string node_1;
		string node_2;

		int src = 0, dest = 0;

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

//	printGraph(graph);
		calcDensity(graph->E, NO_OF_VERTICES, &rho_init);
		init_vertexFlags(graph);

		auto start = chrono::steady_clock::now();
		max_density = maxDensity(graph, rho_init);
		auto end = chrono::steady_clock::now();

		avg_elapsed_time += chrono::duration_cast<chrono::microseconds>(end - start).count();

		if (turn == 0) {
			cout << "Initial Density: " << rho_init << endl
					 << "-----------------------------" << endl
					 << "VERTICES OF DENSEST COMPONENT" << endl
					 << "-----------------------------" << endl;
			dense_nodes = densestComponentVertices(graph);
			cout << "-----------------" << endl
					 << "DENSEST COMPONENT" << endl
					 << "-----------------" << endl;
			densestComponent(graph);
		}
		deallocateGraph(graph);

		ip.close();
	}

	cout << "---------------------------------------" << endl
			 << "Number of nodes in densest component: " << dense_nodes << endl
			 << "Density of densest component: " << max_density << endl
	     << "Average Elapsed time in microseconds : " << avg_elapsed_time << " μs" << endl
			 << "---------------------------------------" << endl;

	exit(0);

}
