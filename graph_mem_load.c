/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include<stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>
//#include<conio.h>
#include <curses.h>
#include <string.h>
#define LEN 256
#define INITIAL_CAPACITY 0

double start, end;

FILE * fr, * fp;
int prev_token;

// Global variables related to the Graph
int E = 0; // number of edges
int V = 4039; // initial number of vertices
int V_count = 4039;
float c_density = 0.0; // graph density
float max_density = 0.0; // maximum density found
float epsilon = 0.25;
float factor = 0;

int d_del, v_del;
omp_lock_t* vertex_lock;

int *queue;
int queueSize = 0;

// A structure to represent an adjacency list node
//struct AdjListNode
//{
//    int dest;
//    struct AdjListNode* next;
//};

// A structure to represent an adjacency list
struct AdjList
{
    int *p_arr;
};

struct NodeDegree
{
    int degree;
};

struct Capacity
{
    int capacity;
};

struct Status
{
    int state;
};

struct Graph
{
    int V;
    struct AdjList* array;
    struct NodeDegree* degrees;
    struct Capacity* capacities;
    struct Status* states;
};

void push(int **arr, int index, int value, int *size, int *capacity){
    
    if((*size + 1) > *capacity){
        *arr = realloc(*arr, ((*size + 1) * sizeof(int)));
        *capacity = (*capacity + 1);
    }
    (*arr)[index] = value;
    *size = *size + 1;
    printf ("size: %d, capacity: %d\n", *size, *capacity);
}

void queuePush(int **arr, int value, int *size){
    
    *arr = realloc(*arr, ((*size + 1) * sizeof(int)));
    *size = (*size + 1);
    
    (*arr)[*size] = value;
    printf ("Queue Size: %d", *size);
}

void addEdge(struct Graph* graph, int src, int dest)
{
    printf ("source: %d , destination: %d\n", src, dest);
    push(&(graph->array[src].p_arr), graph->degrees[src].degree, dest, &(graph->degrees[src].degree), &(graph->capacities[src].capacity));
    
    printf ("source: %d , destination: %d\n", dest, src);
    push(&(graph->array[dest].p_arr), graph->degrees[dest].degree, src, &(graph->degrees[dest].degree), &(graph->capacities[dest].capacity));
    
    E = E + 1;
}

struct Graph* createGraph(int V)
{
    vertex_lock = (omp_lock_t*) malloc(V * sizeof(omp_lock_t));
    
    struct Graph* graph =
    (struct Graph*) malloc(sizeof(struct Graph));
    graph->V = V;
    
    graph->array =
    (struct AdjList*) malloc(V * sizeof(struct AdjList));
    
    graph->degrees =
    (struct NodeDegree*) malloc(V * sizeof(struct NodeDegree));

    graph->capacities =
    (struct Capacity*) malloc(V * sizeof(struct Capacity));
    
    graph->states =
    (struct Status*) malloc(V * sizeof(struct Status));
    
    int i;
    for (i = 0; i < V; ++i) {
        graph->array[i].p_arr = malloc(INITIAL_CAPACITY * sizeof(graph->array[i].p_arr[0]));
        graph->degrees[i].degree = 0;
        graph->capacities[i].capacity = INITIAL_CAPACITY;
        graph->states[i].state = 1;
        omp_init_lock(&vertex_lock[i]);
    }
    return graph;
}

void deleteNode(struct Graph* graph, int node)
{
    //    int linked_node;
    //
    //    struct AdjListNode currentNode = graph->array[node].head;
    //
    //    struct AdjList* currentList = graph->array[node];
    //
    //    while (currentNode){
    //        linked_node = currentNode->dest;
    //        currentNode = currentNode->next;
    //    }
}

void printGraph(struct Graph* graph)
{
    int v;
    int i;
    struct Graph* p_graph = graph;
    for (v = 0; v < p_graph->V; ++v)
    {
        if (p_graph->degrees[v].degree != 0)
        {
            printf("[%d] --> [", v);
            for (i = 0; i < p_graph->degrees[v].degree; ++i)
            {
                if (i != p_graph->degrees[v].degree - 1)
                    printf("%d,", p_graph->array[v].p_arr[i]);
                else
                    printf("%d]\n", p_graph->array[v].p_arr[i]);
            }
        }
    }
}

void getData(struct Graph* graph, char *buff);

int main()
{
    //    fr = fopen ("./results.txt","w");
    fp = fopen("./input/facebook_combined.csv","r");
    //    fp = fopen("./test_file_end.csv","r");
    //    fp = fopen("./artist_edges.csv","r");
    int V = 4039;
    struct Graph* graph = createGraph(V);
    
    factor = 2 * (1 + epsilon);
    
    int buff_size = 20;
    char buff[buff_size];
    int count=0;
    fgets(buff, buff_size, (FILE*)fp);
    
    do
    {
        if(count != 0)
        {
            getData(graph, buff);
        }
        count++;
        
        getch();// check the outputs while debugging
    } while(fgets (buff, buff_size, fp)!= NULL );
    
    fclose (fp);
    fclose (fr);
    
    //printf ("i'm here 1\n");
    
    printGraph(graph);
    
    max_density = (float)V / (float)E;
    
    start = omp_get_wtime();
    
    while(V_count > 0) {
        //printf ("i'm here 2\n");
#pragma omp parallel default(shared)
        {
            // parallel process
            int i;
            //degree = 0;
            //v_del = 0;
            
            for (i = 0; i < V; i++) {
                //printf ("i'm here 3\n");
                //#ifdef DEBUG
                printf ("i value: %d\n", i);
                //#endif
                //omp_set_lock(&vertex_lock[i]);
                int my_state = graph->states[i].state;
                //omp_unset_lock(&vertex_lock[i]);
                
                if (graph->degrees[i].degree <= factor * c_density && my_state == 1) {
                    queuePush(&(queue), i, &(queueSize));
                }
            }
            
#pragma omp for reduction(+: d_del, v_del)
            for (i = 0; i < queueSize; i++) {
                int v = queue[i];
                int count = graph->capacities[v].capacity;
                //#ifdef DEBUG
                printf ("\tcount: %d\n", count);
                //#endif
                // parallel process
                int j;
                for (j = 0; j < count; j++) {
                    // vertex deletion
                    //printf ("i'm here 5\n");
#pragma omp task firstprivate(j)
                    {
                        int check_vertex = graph->array[v].p_arr[j];
                        //#ifdef DEBUG
                        printf ("\tcheck vertex: %d\n", check_vertex);
                        //#endif
                        if (check_vertex != -1) {
                            int cvn_count = graph->capacities[check_vertex].capacity;
                            //#ifdef DEBUG
                            printf ("\tcvn_count: %d\n", cvn_count);
                            //#endif
                            int k;
                            for (k = 0; k < cvn_count; k++) {
                                //printf ("i'm here 6\n");
                                //#ifdef DEBUG
                                printf ("\t\tnow checking: %d\n", graph->array[check_vertex].p_arr[k]);
                                //#endif
                                omp_set_lock(&vertex_lock[check_vertex]);
                                //printf ("i'm here 7\n");
                                //printf ("now checking: %d\n", graph->array[check_vertex].p_arr[k]);
                                //printf ("i'm here 7_x\n");
                                if (graph->array[check_vertex].p_arr[k] == v && graph->array[check_vertex].p_arr[k] != -1 && graph->states[check_vertex].state == 1) {
                                    //printf ("i'm here 8\n");
                                    graph->array[check_vertex].p_arr[k] = -1;
                                    graph->degrees[check_vertex].degree -= 1;
                                    //omp_unset_lock(&vertex_lock[check_vertex]);
                                    //degree += 1;
                                }
                                omp_unset_lock(&vertex_lock[check_vertex]);
                            }
                        }
                    }
                }
                
                omp_set_lock(&vertex_lock[v]);
                //printf ("i'm here 9\n");
                printf ("Freeing node: %d\n", v);
                free(graph->array[v].p_arr);
                graph->states[v].state = 0;
                //printf ("i'm here 10\n");
                omp_unset_lock(&vertex_lock[v]);
                
                v_del += 1;
                d_del += graph->degrees[v].degree;
                
            }
        }
        
#pragma omp master
        {
            printf ("Deleted edges: %d\n", d_del);
            printf ("Deleted vertices: %d\n", v_del);
            E = E - d_del;
            V_count = V_count - v_del;
            
            if (V_count > 0) {
                c_density = (float)E / V_count;
            }
            
            printf ("C density: %f\n", c_density);
            if (c_density > max_density){
                max_density = c_density;
                printf("\nCurrent maximum density = %f\n\n", max_density);
            }
        }
    }
    
    end = omp_get_wtime();
    
    printf("Elapsed time = %f\n", end - start);
}

void getData(struct Graph* graph, char buff[])
{
    char *token = strtok(buff,",");
    int counter=0;
    fr = fopen ("./results.txt","w");
    
    while( token != NULL )
    {
        counter++;
        
        if (counter == 1) prev_token = atoi(token);
        else if (counter == 2) {addEdge(graph, prev_token, atoi(token));}
        
        token = strtok(NULL,",");
    }
}
