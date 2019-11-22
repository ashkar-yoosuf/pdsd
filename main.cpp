#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <iterator>
#include <omp.h>

using namespace std;

#define EPSILON 0.5
#define approxFactor(EPSILON) (2 + 2 * EPSILON)
#define TURNS 1
#define NUM_THREADS 4

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

void addEdge(struct Graph* graph, int src, int dest) {

	push(graph, src,  dest);
	push(graph, dest, src);

}

inline void partitionGraph(unordered_map<int, bool>* active_I, unordered_map<int, bool>* active_II, unordered_map<int, bool>* active_III, unordered_map<int, bool>* active_IV, unordered_set<int>* active, int src, int dest) {

	unordered_set<int> active_deref = *active;
	auto it_end = active_deref.end();

	for (auto it = active_deref.begin(); it != it_end;) {
		active_I->insert({*it, true});
		it = active_deref.erase(it);
		if (it != it_end) {
			active_II->insert({*it, true});
			it = active_deref.erase(it);
			if (it != it_end) {
				active_III->insert({*it, true});
				it = active_deref.erase(it);
				if (it != it_end) {
					active_IV->insert({*it, true});
					it = active_deref.erase(it);
				}
			}
		}
	}

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
	graph_tilde->V = graph->V;
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

float maxDensity(struct Graph* graph, struct GraphTilde* graph_tilde, unordered_map<int, bool> active_I, unordered_map<int, bool> active_II, unordered_map<int, bool> active_III, unordered_map<int, bool> active_IV , double rho_init) {

	int edge_reduction = 0;
	float current_graph_rho = rho_init, current_graphTilde_rho = rho_init;
	bool isTildeChanged {false};

	unordered_map <int, bool> current_failed_I;
	unordered_map <int, bool> current_failed_II;
	unordered_map <int, bool> current_failed_III;
	unordered_map <int, bool> current_failed_IV;

	int edge_dec[4] = {0, 0, 0, 0};
	int vertex_dec[4] = {0, 0, 0, 0};

	while (graph->V > 0) {

#pragma omp parallel
		{
#pragma omp single
			{
#pragma omp task
				{
					for (auto it = active_I.begin(); it != active_I.end();) {
						if (graph->degrees->at(it->first) <= approxFactor(EPSILON) * current_graph_rho) {
							current_failed_I.insert({it->first, true});
							vertex_dec[0]--;
							it = active_I.erase(it);
						} else {
							it++;
						}
					}
				}

#pragma omp task
				{
					for (auto it = active_II.begin(); it != active_II.end();) {
						if (graph->degrees->at(it->first) <= approxFactor(EPSILON) * current_graph_rho) {
							current_failed_II.insert({it->first, true});
							vertex_dec[1]--;
							it = active_II.erase(it);
						} else {
							it++;
						}
					}
				}

#pragma omp task
				{
					for (auto it = active_III.begin(); it != active_III.end();) {
						if (graph->degrees->at(it->first) <= approxFactor(EPSILON) * current_graph_rho) {
							current_failed_III.insert({it->first, true});
							vertex_dec[2]--;
							it = active_III.erase(it);
						} else {
							it++;
						}
					}
				}

				for (auto it = active_IV.begin(); it != active_IV.end();) {
					if (graph->degrees->at(it->first) <= approxFactor(EPSILON) * current_graph_rho) {
						current_failed_IV.insert({it->first, true});
						vertex_dec[3]--;
						it = active_IV.erase(it);
					} else {
						it++;
					}
				}

#pragma omp taskwait

#pragma omp task default(shared) firstprivate(edge_reduction)
				{
					for (auto& it: current_failed_I) {
						Graph::map_type temp = *(graph->sm->at(it.first));
						for (auto& itt: temp) {
							Graph::ar* temp_1 = graph->sm->at(it.first)->at(itt.first);
							Graph::ar* temp_2 = graph->sm->at(itt.first)->at(it.first);
							if ((*temp_1)[0] == 0 && (*temp_2)[0] == 0) {
								(*temp_1)[0] = graph->flag;
								(*temp_2)[0] = graph->flag;
								edge_reduction = (*temp_1)[1];
								graph->degrees->at(it.first) -= edge_reduction;
#pragma omp atomic
								graph->degrees->at(itt.first) -= edge_reduction;
								edge_dec[0] -= edge_reduction;
							}
						}
					}
				}

#pragma omp task default(shared) firstprivate(edge_reduction)
				{
					for (auto& it: current_failed_II) {
						Graph::map_type temp = *(graph->sm->at(it.first));
						for (auto& itt: temp) {
							Graph::ar* temp_1 = graph->sm->at(it.first)->at(itt.first);
							Graph::ar* temp_2 = graph->sm->at(itt.first)->at(it.first);
							if ((*temp_1)[0] == 0 && (*temp_2)[0] == 0) {
								(*temp_1)[0] = graph->flag;
								(*temp_2)[0] = graph->flag;
								edge_reduction = (*temp_1)[1];
								graph->degrees->at(it.first) -= edge_reduction;
#pragma omp atomic
								graph->degrees->at(itt.first) -= edge_reduction;
								edge_dec[1] -= edge_reduction;
							}
						}
					}
				}

#pragma omp task default(shared) firstprivate(edge_reduction)
				{
					for (auto& it: current_failed_III) {
						Graph::map_type temp = *(graph->sm->at(it.first));
						for (auto& itt: temp) {
							Graph::ar* temp_1 = graph->sm->at(it.first)->at(itt.first);
							Graph::ar* temp_2 = graph->sm->at(itt.first)->at(it.first);
							if ((*temp_1)[0] == 0 && (*temp_2)[0] == 0) {
								(*temp_1)[0] = graph->flag;
								(*temp_2)[0] = graph->flag;
								edge_reduction = (*temp_1)[1];
								graph->degrees->at(it.first) -= edge_reduction;
#pragma omp atomic
								graph->degrees->at(itt.first) -= edge_reduction;
								edge_dec[2] -= edge_reduction;
							}
						}
					}
				}

				for (auto& it: current_failed_IV) {
					Graph::map_type temp = *(graph->sm->at(it.first));
					for (auto& itt: temp) {
						Graph::ar* temp_1 = graph->sm->at(it.first)->at(itt.first);
						Graph::ar* temp_2 = graph->sm->at(itt.first)->at(it.first);
						if ((*temp_1)[0] == 0 && (*temp_2)[0] == 0) {
							(*temp_1)[0] = graph->flag;
							(*temp_2)[0] = graph->flag;
							edge_reduction = (*temp_1)[1];
							graph->degrees->at(it.first) -= edge_reduction;
#pragma omp atomic
							graph->degrees->at(itt.first) -= edge_reduction;
							edge_dec[3] -= edge_reduction;
						}
					}
				}


			}
		}

		current_failed_I.clear();
		current_failed_II.clear();

		graph->V += vertex_dec[0] + vertex_dec[1] + vertex_dec[2] + vertex_dec[3];
		graph->E += edge_dec[0] + edge_dec[1] + edge_dec[2] + edge_dec[3];

		vertex_dec[0] = 0, vertex_dec[1] = 0, vertex_dec[2] = 0, vertex_dec[3] = 0;
		edge_dec[0] = 0, edge_dec[1] = 0, edge_dec[2] = 0, edge_dec[3] = 0;

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

	freopen("./output_NEW/[3]fc_4.txt", "w", stdout);

	double avg_elapsed_time = 0;
	float rho_init, max_density = 0;
	int dense_nodes = 0;

	unordered_set<int> active;
	unordered_map<int, bool> active_I;
	unordered_map<int, bool> active_II;
	unordered_map<int, bool> active_III;
	unordered_map<int, bool> active_IV;

	omp_set_num_threads(NUM_THREADS);

	for (int turn = 0; turn < TURNS; turn++) {

		ifstream ip("./input_NEW/[3]facebook_combined.csv");

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

				active.insert(src);
				active.insert(dest);

				addEdge(graph, src, dest);
				graph->E++;

			}
		}

		partitionGraph(&active_I, &active_II, &active_III, &active_IV, &active, src, dest);

		active.clear();

		graph->V = active_I.size() + active_II.size() + active_III.size() + active_IV.size();

		struct GraphTilde* graph_tilde = createGraphTilde(graph);

//		printGraph(graph);
		calcDensity(graph->E, graph->V, &rho_init);

		auto start = chrono::steady_clock::now();
		max_density = maxDensity(graph, graph_tilde, active_I, active_II, active_III, active_IV, rho_init);
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
