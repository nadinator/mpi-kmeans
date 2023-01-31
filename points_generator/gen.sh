#!/bin/sh

read -p 'Clusters: ' k
read -p 'Points per cluster: ' b
read -p 'Output file name: ' output

echo ***GENERATING $b INPUT POINTS EACH IN $k CLUSTERS***
./points_generator.py -c $k -p $b -o $output
