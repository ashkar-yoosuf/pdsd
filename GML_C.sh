#!/bin/sh

gcc -fopenmp -o graph_mem_load.out graph_mem_load.c -lncurses

if [[ $? -eq 0 ]]; then
        ./graph_mem_load.out
fi


