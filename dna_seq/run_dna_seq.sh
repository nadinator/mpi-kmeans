#!/bin/sh

read -p 'Clusters: ' k
read -p 'Strands per cluster: ' p
read -p 'Strand lengths: ' l
# read -p 'Input file: ' f
# read -p 'KMeans iterations: ' iters
echo 'running...\n'

tp=`expr $k \* $p` # Total points (for dna_seq.c)

# Generate random data based on input
python3 ../dna_generator/data_generator.py -c $k -p $p -l $l;
# Compile program
gcc dna_seq.c -o dna_seq;
# Run program based on the same input
./dna_seq -c $k -t $tp -l $l -i ../dna_generator/dna_data.csv