#!/bin/sh

read -p 'Clusters: ' k
read -p 'Strands per cluster: ' p
read -p 'Strand lengths: ' l
read -p 'MPI nodes: ' n
# read -p 'Input file: ' f
# read -p 'KMeans iterations: ' iters
echo 'running...\n'

tp=`expr $k \* $p`

# Generate random data based on input
python3 ../dna_generator/data_generator.py -c $k -p $p -l $l;
# Compile program 
mpicc dna_par.c -o dna_par;
# Run
mpiexec -n $n ./dna_par -c $k -t $tp -l $l -i ../dna_generator/dna_data.csv