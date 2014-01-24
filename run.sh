#!/bin/sh
#
# the next line is a "magic" comment that tells codine to use bash
#$ -S /bin/bash
#
# and now for some real work

./ts 0 0.0 sensitivity.txt cell_names.txt 10000 0 $PBS_ARRAYID
#./ts 0 -1.0 sensitivity.txt cell_names.txt 10000 0 $PBS_ARRAYID
#./ts 1 -1.0 sensitivity.txt cell_names.txt 50000 0 $PBS_ARRAYID
