#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <unordered_map>
#include <map>
#include <iterator>

using namespace std;

#define EPSILON 0.5
#define approxFactor(EPSILON) (2 + 2 * EPSILON)
#define TURNS 1


struct Graph
{
	typedef array<int, 2> ar;
	typedef unordered_map<int, ar*> map_type;
	typedef unordered_map<int, map_type*> super_map;
	super_map* sm;
	unordered_map<int, int>* degrees;
	int V;
	int E;
	int flag;
};

struct GraphTilde
{
	int V;
	int E;
};

inline void push(struct Graph* graph, int pos, int value) {

	auto it = graph->sm->find(pos);

	if (it != graph->sm->end()) {

		auto itt = graph->sm->at(pos)->find(value);

		if (itt != graph->sm->at(pos)->end()) {

			Graph::ar* flag_subDeg = graph->sm->at(pos)->at(value);
			(*flag_subDeg)[1]++;

		} else {

			auto* ar = new Graph::ar();
			(*ar)[0] = 0; (*ar)[1] = 1;
			graph->sm->at(pos)->emplace(make_pair(value, ar));

		}

		graph->degrees->at(pos)++;

	} else {

		auto* mt = new Graph::map_type();
		auto* ar = new Graph::ar();
		(*ar)[0] = 0; (*ar)[1] = 1;
		graph->sm->emplace(make_pair(pos, mt));
		graph->sm->at(pos)->emplace(make_pair(value, ar));
		graph->degrees->insert({pos, 1});

	}
}

//inline void push(struct Graph* graph, int pos, int value) {
//
//	auto it = graph->sm.find(pos);
//
//	if (it != graph->sm.end()) {
//
//		graph->sm[pos]->insert({value, 0});
//		graph->degrees[pos]++;
//
//	} else {
//
//		auto* mt = new Graph::map_type();
//		graph->sm.emplace(make_pair(pos, mt));
//		graph->sm[pos]->emplace(make_pair(value, 0));
//		graph->degrees.insert({pos, 1});
//
//	}
//}

void addEdge(struct Graph* graph, int src, int dest) {

	push(graph, src,  dest);
	push(graph, dest, src);

}

struct Graph* createGraph() {

	static auto* graph = new Graph();
	graph->V = 0;
	graph->E = 0;
	graph->flag = -1;
	graph->sm = new Graph::super_map();
	graph->degrees = new unordered_map<int, int>();

	return graph;
}

struct GraphTilde* createGraphTilde(struct Graph* graph) {

	auto* graph_tilde = new GraphTilde();
	graph_tilde->V = 0;
	graph_tilde->E = graph->E;

	return graph_tilde;

}

//void printGraph(struct Graph* graph) {
//
//	for (int v = 0; v < NO_OF_VERTICES; ++v)
//	{
//		if (graph->capacities[v] != 0)
//		{
//			cout << "[" << v << "] --> [";
//			for (int i = 0; i < graph->capacities[v]; ++i)
//			{
//				if (i != graph->capacities[v] - 1)
//					cout << graph->arrayOfArrays[v][i] << ",";
//				else
//					cout << graph->arrayOfArrays[v][i] << "] | degree: " << graph->degrees[v] << endl;
//			}
//		}
//	}
//}

inline void calcDensity(int num_edges, int num_vertices, float* density) {

	*density = (double) num_edges/  (double) num_vertices;

}

inline void assignToTilde(struct Graph* graph, struct GraphTilde* graph_tilde) {

	graph_tilde->V = graph->V;
	graph_tilde->E = graph->E;

}

void deallocateGraph(struct Graph* graph) {

	for (auto it = graph->sm->begin(); it != graph->sm->end();) {

		Graph::map_type temp = *(graph->sm->at(it->first));

		for (auto itt = temp.begin(); itt != temp.end();){

			delete temp[itt->first];
			itt = temp.erase(itt);

		}

		delete graph->sm->at(it->first);
		it = graph->sm->erase(it);

	}

	delete graph->sm;
	delete graph->degrees;
	delete graph;

}

inline void deallocateGraphTilde(struct GraphTilde* graph_tilde) {

	delete graph_tilde;

}

float maxDensity(struct Graph* graph, struct GraphTilde* graph_tilde, unordered_map<int, bool> active , double rho_init) {

	int edge_reduction = 0;
	float current_graph_rho = rho_init, current_graphTilde_rho = rho_init;
	bool isTildeChanged {false};

	unordered_map <int, bool> current_failed;

	while (graph->V > 0) {

		for (auto it = active.begin(); it != active.end();) {

			if (graph->degrees->at(it->first) <= approxFactor(EPSILON) * current_graph_rho) {

				current_failed.insert({it->first, true});
				graph->V--;
				it = active.erase(it);

			} else {

				it++;

			}
		}

		for (auto& it: current_failed) {

			Graph::map_type temp = *(graph->sm->at(it.first));

			for (auto& itt: temp) {

				Graph::ar* temp_1 = graph->sm->at(it.first)->at(itt.first);
				Graph::ar* temp_2 = graph->sm->at(itt.first)->at(it.first);

				if ((*temp_1)[0] == 0 && (*temp_2)[0] == 0) {

					(*temp_1)[0] = graph->flag;
					(*temp_2)[0] = graph->flag;
					edge_reduction = (*temp_1)[1];
					graph->degrees->at(it.first) -= edge_reduction;
					graph->degrees->at(itt.first) -= edge_reduction;
					graph->E -= edge_reduction;

				}
			}
		}

		current_failed.clear();

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

	return current_graphTilde_rho;

}

//void densestComponent(struct Graph* graph) {
//
//	for (int v = 0; v < NO_OF_VERTICES; v++) {
//		if (graph->failed_V_arr[v] == graph->flag) {
//			cout << "[" << v << "] --> |";
//			for (int i = 0; i < graph->capacities[v]; ++i) {
//				if (graph->flags[v][i] == graph->flag) {
//					cout << graph->arrayOfArrays[v][i] << "|";
//				}
//			}
//			cout << endl;
//		}
//	}
//
//}

//int densestComponentVertices(struct Graph* graph) {
//
//	int dense_nodes = 0;
//	for (int v = 0; v < NO_OF_VERTICES; v++) {
//		if (graph->failed_V_arr[v] == graph->flag) {
//			dense_nodes++;
//			cout << v << endl;
//		}
//	}
//
//	return dense_nodes;
//
//}

int main() {

	freopen("/home/ashkar/Documents/Graph_Output/slj1_SERIAL.txt", "w", stdout);

	double avg_elapsed_time = 0;
	float rho_init, max_density = 0;
	int dense_nodes = 0;

	unordered_map<int, bool> active;

	for (int turn = 0; turn < TURNS; turn++) {

//		ifstream ip("./input/as-skitter_reFormatted.csv");
		ifstream ip("/home/ashkar/Documents/Graph_Input/soc-LiveJournal1.csv");

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

			if (line == 2 && !node_1.empty() && src != dest) {

				active.insert({src, true});
				active.insert({dest, true});
				addEdge(graph, src, dest);
				graph->E++;

			}
		}

		struct GraphTilde* graph_tilde = createGraphTilde(graph);
		graph->V = active.size();

//		printGraph(graph);
		calcDensity(graph->E, graph->V, &rho_init);

		auto start = chrono::steady_clock::now();
		max_density = maxDensity(graph, graph_tilde, active, rho_init);
		auto end = chrono::steady_clock::now();

		avg_elapsed_time += chrono::duration_cast<chrono::microseconds>(end - start).count();

		if (turn == 0) {

			cout << "Initial Density: " << rho_init << endl
					 << "-----------------------------" << endl
					 << "VERTICES OF DENSEST COMPONENT" << endl
					 << "-----------------------------" << endl;
//			dense_nodes = densestComponentVertices(graph);
			cout << "-----------------" << endl
					 << "DENSEST COMPONENT" << endl
					 << "-----------------" << endl;
//			densestComponent(graph);

		}

		deallocateGraph(graph);
		deallocateGraphTilde(graph_tilde);

		ip.close();
	}

	cout << "---------------------------------------" << endl
			 << "Number of nodes in densest component: " << dense_nodes << endl
			 << "Density of densest component: " << max_density << endl
			 << "Average Elapsed time in microseconds : " << avg_elapsed_time / (float) TURNS << " μs" << endl
			 << "---------------------------------------" << endl;

	exit(0);

}
