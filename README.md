pydt
====

Pitman Yor diffusion tree sampler

To compile you will need the Boost libraries. e.g. if boost is in ~/boost_1_51_0/ then run: 

g++ -O3 main.cpp -I ~/boost_1_51_0/ -o pydt

To run on the CCLE data provided here you can then run for example

./pydt 0 1.5 sensitivity.txt cell_names.txt 1000 1 1

Just running pydt with no args will explain what these inputs are. 