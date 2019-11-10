#include <iostream>
#include "maxflow.cpp"
using namespace std;

/*
 * Main Contains Menu
 */
int main(){
    int graph[V][V] = { {0, 16, 13, 0, 0, 0},
                        {0, 0, 10, 12, 0, 0},
                        {0, 4, 0, 0, 14, 0},
                        {0, 0, 9, 0, 0, 20},
                        {0, 0, 0, 7, 0, 4},
                        {0, 0, 0, 0, 0, 0}
    };
    minCut(graph, 0, 5);
    return 0;
}