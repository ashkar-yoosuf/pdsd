//
//  ParAlgo.cpp
//
//
//  Created by Dinuka Manohara De Zoysa on 11/30/19.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <iterator>
#include <vector>
#include <deque>
#include <string>
#include <boost/tokenizer.hpp>
#include <omp.h>

using namespace std;

#define NUM_THREADS 4

typedef vector<long> long_vector;
typedef unordered_map<long, long> long_dict;
typedef unordered_map<long, long_dict*> l_d_dict;
typedef unordered_set<long> long_set;

struct Graph
{
    l_d_dict* sm;
    long_dict* deg;
    long_vector* vert;
    long V;
    long E;
};

inline void push(struct Graph* graph, int pos, int value) {

    auto it = graph->sm->find(pos);

    if (it != graph->sm->end()) {
        auto itt = graph->sm->at(pos)->find(value);
		if (itt != graph->sm->at(pos)->end()) {
			graph->sm->at(pos)->at(value)++;
		} else {
			graph->sm->at(pos)->emplace(make_pair(value, 1));
		}
        graph->deg->at(pos)++;
    } else {
        auto* mt = new ::long_dict();
        graph->sm->emplace(make_pair(pos, mt));
        graph->sm->at(pos)->emplace(make_pair(value, 1));
        graph->deg->insert({pos, 1});
        graph->vert->push_back(pos);
    }
}

void addEdge(struct Graph* graph, int src, int dest) {

    push(graph, src,  dest);

}

struct Graph* createGraph() {

    static auto* graph = new Graph();
    graph->V = 0;
    graph->E = 0;
    graph->sm = new ::l_d_dict();
    graph->deg = new ::long_dict();
    graph->vert = new ::long_vector();

    return graph;
}

inline void calcDensity(long num_edges, long num_vertices, float* density) {

    *density = (double) num_edges/  (double) num_vertices;

}

void deallocateGraph(struct Graph* graph) {

    for (auto it = graph->sm->begin(); it != graph->sm->end();) {

        auto temp = *(graph->sm->at(it->first));

        delete graph->sm->at(it->first);
        it = graph->sm->erase(it);

    }

    delete graph->sm;
    delete graph->deg;
    delete graph->vert;
    delete graph;

}

float maxDensity(struct Graph* graph) {
   
    float mden = 0;
    long mv = 0;
    long me = 0;
    float den = 0;
    long mcore = 0;
    long deleted = 0;
    long temp = 0;
    long visited = 0;
    long additional_e = 0;
    long_vector *legit_buff = new ::long_vector();

    #pragma omp parallel
    {
        long k = 0;
        long v = 0;
        long my_deleted = 0;
        long my_temp = 0;
        float den = 0;

        long *buff = (long *)malloc( (graph->V*sizeof(long)) / NUM_THREADS );

        long start = 0, end = 0;

        while( visited < graph->V ) {

    #pragma omp for schedule(static) 
            for(long i = 0; i < graph->V; i++)  {
                if( graph->deg->at(graph->vert->at(i)) == k ) {
                    buff[end] = graph->vert->at(i);
                    end ++;
                }
            }

            //Get work from curr queue and also add work after the current size
            while( start < end ) {

                v = buff[start];
                start ++;
                my_deleted += k;
                for(auto p = graph->sm->at(v)->begin() ; p != graph->sm->at(v)->end(); p++) {
                    auto u = p->first;
                    for(int i = 0; i < p->second; i++) {
                        long deg_u = graph->deg->at(u);

                        if( deg_u > k ) {
                            int du = __sync_fetch_and_sub(&graph->deg->at(u), 1);

                            if( du == (k+1) ) {
                                buff[end] = u;
                                end ++;
                            }

                            if( du <= k ) {
                                __sync_fetch_and_add(&graph->deg->at(u), 1);
                            }

                            if (du > k) {
                                my_temp += 1;
                            } else {
                                my_deleted += 1;
                                my_temp -= 1;
                            }
                        } //deg_u > k
                    }
                } //visit adjacencies
            }  //end of while loop

            __sync_fetch_and_add(&visited, end);
            __sync_fetch_and_add(&deleted, my_deleted);
            __sync_fetch_and_add(&temp, my_temp);

    #pragma omp barrier 
            start = 0;
            end = 0;
            my_deleted = 0;
            my_temp = 0;
            k += 1;


        #pragma omp single 
        {
            calcDensity(graph->E - ((deleted + temp) / 2), graph->V - visited, &den);

            if (den > mden) {
                mden = den;
                me = graph->E - ((deleted + temp) / 2);
                mv = graph->V - visited;
                mcore = k - 1;
            }
        }
            
        } //end of main while loop
        
        //free(buff);

        long legits = 0;
        long legit_v;
        //long *buff = (long *)malloc( (graph->V*sizeof(long)) / NUM_THREADS );

#pragma omp for schedule(static) 
        for(long i = 0; i < graph->V; i++)  {
            if((graph->deg->at(graph->vert->at(i)) < mcore) && (graph->deg->at(graph->vert->at(i)) > mden)) {
                buff[end] = graph->vert->at(i);
                end ++;
            }
        }

        while( start < end ) {
            legits = 0;
            v = buff[start];
            start ++;
            for(auto p = graph->sm->at(v)->begin() ; p != graph->sm->at(v)->end(); p++) {
                
                auto u = p->first;

                if (graph->deg->at(u) >= mcore) {
                    legits += p->second;
                }
            } //visit adjacencies

            if (legits > mden) {
                #pragma omp critical
                {
                    additional_e += legits;
                    legit_buff->push_back(v);
                }
            }
        } //end of while loop

        free(buff);

#pragma omp barrier 

#pragma omp for reduction(+: additional_e) schedule(dynamic)
        for (long i = 0; i < legit_buff->size() - 1; i++) {
            legit_v = legit_buff->at(i);
            graph->deg->at(legit_v) = mcore; // to make it easy to retrieve densest subgraph nodes
            for (long j = i + 1; j < legit_buff->size(); j++) {
                if (graph->sm->at(legit_v)->find(legit_buff->at(j)) != graph->sm->at(legit_v)->end()) {
                    additional_e += graph->sm->at(legit_v)->at(legit_buff->at(j));
                }
            }
        }

#pragma omp barrier
    }//end of parallel region

    me += additional_e;
    mv += legit_buff->size();
    calcDensity(me, mv, &mden); 

    delete legit_buff;

    return mden;
}

int main() {

    int max_threads = omp_get_max_threads();

    cout << "max threads: " << max_threads << endl;
    cout << "requested threads: " << NUM_THREADS << endl;

    if (NUM_THREADS <= max_threads) {

        double elapsed_time = 0;
        float rho_init, max_density = 0;
        int dense_nodes = 0;

        omp_set_num_threads(NUM_THREADS);

        //ifstream ip("as-skitter.csv");

        struct Graph* graph = createGraph();

        string node_1;
		string node_2;

		ifstream ip;
		ip.open("as-skitter.csv");
		string line;
		unsigned long long line_id = 0;
		unsigned long long src = 0, dest = 0;
		int count = 0;

		if (ip.is_open()) {
			while ((!ip.eof())) {
				getline(ip, line);
				boost::tokenizer<separator_type> tokenizer(line, separator_type(" "));
				auto it = tokenizer.begin();

				while (it != tokenizer.end()) {

					if (count == 0) {
						if (it == tokenizer.begin()) {

							try {
								stoull(*it);
								count++;
								src = ++line_id - 1;
							} catch (exception& e) {
								it = tokenizer.end();
							}

						}
					} else {
						if (line_id == 1) {

							graph->V = stoull(*it++);
							graph->E = stoull(*it);
							it = tokenizer.end();

						} else if (line_id > 1) {

							dest = stoull(*it++);
							if (src != dest) {
								active.insert(src);
								addEdge(graph, src, dest);
							} else graph->E--;
						}
					}
				}
				count = 0;
			}
		}
        

        cout << "Here, after File Reading" << endl;

        auto start = chrono::steady_clock::now();
        

        float mden = maxDensity(graph);

        
        auto end = chrono::steady_clock::now();

        elapsed_time += chrono::duration_cast<chrono::microseconds>(end - start).count();

        cout << "Maximum Density: " << mden << endl
             << "Elapsed time in seconds : " << elapsed_time / 1000000 << " s" << endl;

        /*
        cout << "Initial Density: " << rho_init << endl;
        cout << "-----------------" << endl
             << "DENSEST COMPONENT" << endl
             << "-----------------" << endl;
        */

        deallocateGraph(graph);

        ip.close();

        /*
        cout << "---------------------------------------" << endl
             << "Number of nodes in densest component: " << dense_nodes << endl
             << "Density of densest component: " << max_density << endl
             << "Elapsed time in microseconds : " << elapsed_time << " μs" << endl
             << "---------------------------------------" << endl;
        */
    } else
        cout << "Thread limit exceeded" << endl;

    exit(0);
}

