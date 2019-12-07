void PKC_org(graph_t *g, int *deg) {
    long n = g->n;
    long visited = 0;

#pragma omp parallel
{
    int level = 0;

    vid_t *buff = (vid_t *)malloc( (n*sizeof(vid_t)) / NUM_THREADS );
    assert( buff != NULL);

    long start = 0, end = 0;

    #pragma omp for schedule(static)
    for (long i = 0; i < n; i++)  {
        deg[i] = (g->num_edges[i+1] - g->num_edges[i]);
    }

    while( visited < n ) {

        #pragma omp for schedule(static)
        for(long i = 0; i < n; i++)  {

            if( deg[i] == level ) {
                buff[end] = i;
                end ++;
            }
        }

        //Get work from curr queue and also add work after the current size
        while( start < end ) {

            vid_t v = buff [start];
            start ++;

            for(eid_t j = g->num_edges[v]; j < g->num_edges[v+1]; j++) {
                vid_t u = g->adj[j];

                    int deg_u = deg[u];

                    if( deg_u > level ) {
                        int du = __sync_fetch_and_sub(&deg[u], 1);

                        if( du == (level+1) ) {
                            buff[end] = u;
                            end ++;
                        }

                        if( du <= level ) {
                            __sync_fetch_and_add(&deg[u], 1);
                        }
                    }  //deg_u > level

            } //visit adjacencies
        }  //end of while loop

        __sync_fetch_and_add(&visited, end);

#pragma omp barrier
        start = 0;
        end = 0;
        level = level + 1;

    }   //end of #visited < n

    free( buff );

}  //#end of parallel region

}
