#!/bin/sh

read -p 'Clusters: ' k
read -p 'Points per cluster: ' pk
read -p 'MPI Nodes: ' n
# read -p 'Input file: ' f
# read -p 'KMeans iterations: ' iters
echo 'running...\n'

tp=$(($k * $pk))

# Compile
mpicc points_par.c -o points_par;
# Generate points
python3 ../points_generator/randomclustergen/generaterawdata.py -c $k -p $pk -o ../points_generator/input/cluster.csv 
# Run
mpiexec -n $n ./points_par -c $k -t $tp -i ../points_generator/input/cluster.csv