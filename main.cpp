#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <iterator>
#include <omp.h>

using namespace std;

#define EPSILON 0.5
#define approxFactor(EPSILON) (2 + 2 * EPSILON)
#define NUM_THREADS 1

typedef unordered_map<unsigned int, unsigned int> map_type;
typedef unordered_map<unsigned int, map_type*> super_map;
typedef unordered_set<unsigned int> set_type;
typedef array<set_type*, NUM_THREADS> nodes_track;

struct Graph
{
	super_map* sm;
	map_type* degrees;
	unsigned int V;
	unsigned long long E;
};

struct GraphTilde
{
	unsigned int V;
	unsigned long long E;
};

inline void partitionNodes(set_type* active, nodes_track* _actives, nodes_track* _fails)
{
	for (int i = 0; i < NUM_THREADS; i++) {
		(*_actives)[i] = new ::set_type();
		(*_fails)[i] = new ::set_type();
	}

	set_type active_deref = *active;
	auto it_end = active_deref.end();

	for (auto it = active_deref.begin(); it != it_end;) {
		for (int i = 0; i < NUM_THREADS; i++) {
			if (it != it_end) {
				(*_actives)[i]->insert(*it);
				it = active_deref.erase(it);
			}
		}
	}
}

inline void push(struct Graph* graph, unsigned int pos, unsigned int value) {

	auto it = graph->sm->find(pos);

	if (it != graph->sm->end()) {

		auto itt = graph->sm->at(pos)->find(value);

		if (itt != graph->sm->at(pos)->end()) {

			graph->sm->at(pos)->at(value)++;

		} else {

			graph->sm->at(pos)->emplace(make_pair(value, 1));

		}

		graph->degrees->at(pos)++;

	} else {

		auto* mt = new ::map_type();
		graph->sm->emplace(make_pair(pos, mt));
		graph->sm->at(pos)->emplace(make_pair(value, 1));
		graph->degrees->insert({pos, 1});

	}
}

void addEdge(struct Graph* graph, unsigned int src, unsigned int dest) {

	push(graph, src,  dest);
	push(graph, dest, src);

}

struct Graph* createGraph() {

	static auto* graph = new Graph();
	graph->V = 0;
	graph->E = 0;
	graph->sm = new ::super_map();
	graph->degrees = new ::map_type();

	return graph;
}

struct GraphTilde* createGraphTilde(struct Graph* graph) {

	auto* graph_tilde = new GraphTilde();
	graph_tilde->V = 0;
	graph_tilde->E = graph->E;

	return graph_tilde;

}

inline void calcDensity(unsigned long long num_edges, unsigned int num_vertices, float* density) {

	*density = (double) num_edges/  (double) num_vertices;

}

inline void assignToTilde(struct Graph* graph, struct GraphTilde* graph_tilde) {

	graph_tilde->V = graph->V;
	graph_tilde->E = graph->E;

}

void deallocateGraph(struct Graph* graph) {

	for (auto it = graph->sm->begin(); it != graph->sm->end();) {

		map_type temp = *(graph->sm->at(it->first));

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

inline void deallocateMinor(nodes_track* _actives, nodes_track* _fails) {

	for (int i = 0; i < NUM_THREADS; i++) {
		delete (*_actives)[i];
		delete (*_fails)[i];
	}
	delete _actives;
	delete _fails;
}

float maxDensity(struct Graph* graph, struct GraphTilde* graph_tilde, nodes_track* _actives, nodes_track* _fails, double rho_init) {

	unsigned int edge_reduction = 0, id = 0, iid = 0;
	unsigned int* temp_1, * temp_2;
	map_type temp;
	float current_graph_rho = rho_init, current_graphTilde_rho = rho_init;
	bool isTildeChanged {false};

	auto* edge_dec = new array<int, NUM_THREADS>();
	auto* vertex_dec = new array<int, NUM_THREADS>();

	while (graph->V > 0) {

#pragma omp parallel
		{
#pragma omp task default(shared)
			{
				int _id;
#pragma omp critical
				{
					_id = id;
					id++;
				}
				for (auto it = (*_actives)[_id]->begin(); it != (*_actives)[_id]->end();) {

					if (graph->degrees->at(*it) <= approxFactor(EPSILON) * current_graph_rho) {

						(*_fails)[_id]->insert(*it);
						(*vertex_dec)[_id]--;
						it = (*_actives)[_id]->erase(it);

					} else {

						it++;

					}
				}
			}
		}

#pragma omp parallel
		{
#pragma omp task default(shared) firstprivate(edge_reduction, temp, temp_1, temp_2)
			{
				int _id;
#pragma omp critical
				{
					_id = iid;
					iid++;
				}

				set_type st = *(*_fails)[_id];
				for (auto& it: st) {

					temp = *(graph->sm->at(it));

					for (auto& itt: temp) {

						temp_1 = &(graph->sm->at(it)->at(itt.first));
						temp_2 = &(graph->sm->at(itt.first)->at(it));

						if (*temp_1 != 0 && *temp_2 != 0) {

							*temp_2 = 0;
							edge_reduction = temp.at(itt.first);
							(*edge_dec)[_id] -= edge_reduction;
							__sync_sub_and_fetch(&(graph->degrees->at(itt.first)), edge_reduction);

						}
					}
				}
			}
		}

		id = 0, iid = 0;

		for (int i = 0; i < NUM_THREADS; i++) {
			graph->V += (*vertex_dec)[i];
			(*vertex_dec)[i] = 0;
			graph->E += (*edge_dec)[i];
			(*edge_dec)[i] = 0;
		}

		if (graph->V != 0) {
			for (int i = 0; i < NUM_THREADS; i++) {
				(*_fails)[i]->clear();
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

	delete edge_dec;
	delete vertex_dec;

	return current_graphTilde_rho;

}

int densestComponent(nodes_track* _fails) {

	unsigned int num_nodes = 0;
	for (int i = 0; i < NUM_THREADS; i++) {
		num_nodes += (*_fails)[i]->size();
		for (auto& it: *(*_fails)[i]) {
			cout << it << endl;
		}
	}

	return num_nodes;
}

int main() {

	unsigned int max_threads = omp_get_max_threads();

	cout << "max threads: " << max_threads << endl;
	cout << "requested threads: " << NUM_THREADS << endl;

	if (NUM_THREADS <= max_threads) {

		double elapsed_time = 0;
		float rho_init, max_density = 0;
		unsigned int dense_nodes = 0;

		omp_set_num_threads(NUM_THREADS);

		set_type active;
		auto* fails = new nodes_track();
		auto* active_ar = new nodes_track();

		ifstream ip("./input/[3]facebook_combined.csv");

		struct Graph* graph = createGraph();

		string node_1;
		string node_2;

		unsigned int src = 0, dest = 0;

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

		partitionNodes(&active, active_ar, fails);

		struct GraphTilde* graph_tilde = createGraphTilde(graph);

		graph->V = active.size();

		active.clear();

		calcDensity(graph->E, graph->V, &rho_init);

		auto start = chrono::steady_clock::now();
		max_density = maxDensity(graph, graph_tilde, active_ar, fails, rho_init);
		auto end = chrono::steady_clock::now();

		elapsed_time += chrono::duration_cast<chrono::microseconds>(end - start).count();

		cout << "Initial Density: " << rho_init << endl;
	 	cout << "-----------------" << endl
				 << "DENSEST COMPONENT" << endl
				 << "-----------------" << endl;
		dense_nodes = densestComponent(fails);

		deallocateMinor(active_ar, fails);
		deallocateGraph(graph);
		deallocateGraphTilde(graph_tilde);

		ip.close();

		cout << "---------------------------------------" << endl
				 << "Number of nodes in densest component: " << dense_nodes << endl
				 << "Density of densest component: " << max_density << endl
				 << "Elapsed time in microseconds : " << elapsed_time << " μs" << endl
				 << "---------------------------------------" << endl;
	} else
		cout << "Thread limit exceeded" << endl;

	exit(0);
}
