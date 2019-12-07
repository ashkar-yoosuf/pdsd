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
#include <omp.h>

using namespace std;

#define NUM_THREADS 1

typedef vector<long> long_vector;
typedef unordered_map<long, long> long_dict;
typedef unordered_map<long, long_vector*> l_v_dict;
typedef unordered_set<long> long_set;
typedef unordered_map<long, deque<long>*> l_ld_dict;
//typedef unordered_map<long, long_set*> bucket;

struct Graph
{
    l_v_dict* sm;
    long_dict* deg;
    long_vector* vert;
    int V;
    int E;
};

inline void push(struct Graph* graph, int pos, int value) {

    auto it = graph->sm->find(pos);

    if (it != graph->sm->end()) {
        graph->sm->at(pos)->push_back(value);
        graph->deg->at(pos)++;
    } else {
        auto* mt = new ::long_vector();
        graph->sm->emplace(make_pair(pos, mt));
        graph->sm->at(pos)->push_back(value);
        graph->deg->insert({pos, 1});
        graph->vert->push_back(pos);
        graph->V++;
    }
}

void addEdge(struct Graph* graph, int src, int dest) {

    push(graph, src,  dest);
    push(graph, dest, src);

}

struct Graph* createGraph() {

    static auto* graph = new Graph();
    graph->V = 0;
    graph->E = 0;
    graph->sm = new ::l_v_dict();
    graph->deg = new ::long_dict();
    graph->vert = new ::long_vector();

    return graph;
}

inline void calcDensity(long num_edges, long num_vertices, float* density) {

    *density = (double) num_edges/  (double) num_vertices;

}

void deallocateGraph(struct Graph* graph) {

    for (auto it = graph->sm->begin(); it != graph->sm->end();) {

        long_vector temp = *(graph->sm->at(it->first));

        delete graph->sm->at(it->first);
        it = graph->sm->erase(it);

    }

    delete graph->sm;
    delete graph->deg;
    delete graph->vert;
    delete graph;

}

float maxDensity(struct Graph* graph) {
    long md = 0;

    cout << graph->V << endl;

    #pragma omp parallel for reduction(max: md)
    for (long i = 0; i < graph->V; i++) {
        long deg = graph->deg->at(graph->vert->at(i));
        if (deg > md) {
            md = deg;
        }
    }

    //cout << "Here, after Max Degree" << endl;

    long deg_bucket[md + 3] = {0};
    #pragma omp parallel for reduction(+: deg_bucket) schedule(dynamic)
    for (long i = 0; i < graph->V; i++) {
        long deg = graph->deg->at(graph->vert->at(i));
        for (long j = deg + 2; j < md + 2; j++) {
            deg_bucket[j]++;
        }
    }

    //cout << "Here, after Bucket Init" << endl;

    long vert[graph->V];

    #pragma omp parallel for
    for (long i = 0; i < graph->V; i++) {
        long v = graph->vert->at(i);
        long d;
        #pragma omp critical
        {
            d = graph->deg->at(v);
            vert[deg_bucket[d + 1]] = v;
            deg_bucket[d + 1]++;
        }
        
    }

    //cout << "Here, after Categorizing" << endl;

    long k = 0;
    int b_count = 0;
    long edges = graph->E;
    long vertices = graph->V;
    //long e = 0;
    long ver = 0;
    float mden = 0;
    long mv = 0;
    long me = 0;
    float den = 0;
    long deleted = 0;
    long temp = 0;
    l_ld_dict k_buffers;

    while (k < md) {
        //edges = edges - e;
        //vertices = vertices - n;
        cout << edges - ((deleted + temp) / 2) << ", " << vertices - ver << endl;
        calcDensity(edges - ((deleted + temp) / 2), vertices - ver, &den);

        if (den > mden) {
            mden = den;
            me = edges - ((deleted + temp) / 2);
            mv = vertices;
        }

        //e = 0;
        //n = 0;

        #pragma omp parallel for reduction(+: deleted, temp, ver)
        for (long i = deg_bucket[k]; i < deg_bucket[k + 1]; i++) {
            deque<long> buff;
            long v = vert[i];
            long d;
            deleted += k;
            for (auto nei = graph->sm->at(v)->begin(); nei != graph->sm->at(v)->end(); nei++) {
                if (graph->deg->at(*nei) > k) {
                    d = __sync_fetch_and_sub(&(graph->deg->at(*nei)), 1);
                    
                    if (vert[k] != *nei) {
                        #pragma omp critical
                        {
                            long *a = &(vert[k]);
                            long b = *nei;
                            vert[k] = *nei;
                        }
                    }
                    if (d == k + 1) {
                        buff.push_back(*nei);
                    }
                    if (d <= k) {
                        __sync_fetch_and_add(&(graph->deg->at(*nei)), 1);
                    }
                    temp += 1;
                }else {
                    deleted += 1;
                    temp -= 1;
                }
            }

            while (!buff.empty()) {
                v = buff.front();
                buff.pop_front();
                deleted += k;
                for (auto nei = graph->sm->at(v)->begin(); nei != graph->sm->at(v)->end(); nei++) {
                    if (graph->deg->at(*nei) > k) {
                        d = __sync_fetch_and_sub(&(graph->deg->at(*nei)), 1);
                        if (d == k + 1) {
                            buff.push_back(*nei);
                        }
                        if (d <= k) {
                            __sync_fetch_and_add(&(graph->deg->at(*nei)), 1);
                        }
                        temp += 1;
                    }else {
                        deleted += 1;
                        temp -= 1;
                    }
                }
                ver += 1;
            }
        }

        v += deg_bucket[k + 1] - deg_bucket

        k++;
    }

    //cout << "ME : " << me << ", MV : " << mv << endl;
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

        ifstream ip("ca-GrQc-new.csv");

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
                addEdge(graph, src, dest);
                graph->E++;
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

