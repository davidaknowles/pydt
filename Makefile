#!/bin/bash 
#g++ main.cpp -g -I ~/boost_1_47_0/ -o testps.bin
g++ -O3 main.cpp -I ~/scratch/boost_1_51_0/ -o ts
#./testps.bin
#python plot_ddt.py train1.txt tree.txt