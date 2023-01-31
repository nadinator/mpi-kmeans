#!/bin/sh

read -p 'Clusters: ' k
read -p 'Total Points: ' tp
# read -p 'Input file: ' f
# read -p 'KMeans iterations: ' iters

gcc points_seq.c -o points_seq;
python3 ../points_generator/randomclustergen/generaterawdata.py -c $k -p $tp -o ../points_generator/input/cluster.csv # Generate new random data
./points_seq -c $k -t $tp -i ../points_generator/input/cluster.csv