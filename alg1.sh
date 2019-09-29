#!/bin/sh

gcc -o main.out main.cpp

if [[ $? -eq 0 ]]; then
        ./main.out
fi
