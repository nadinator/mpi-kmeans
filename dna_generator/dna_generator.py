# Author: Nadim Bou Alwan
# Date: Mon Oct 10 2022

""" 
DESCRIPTION
This script generates a file containing (clusters*strands) DNA strands,
where each strand is centered around some cluster and is a maximum of DIST 
distance away from that cluster's centroid.

USAGE
DON'T run this file on its own. Must be run from the /dna_seq directory to work. See run_dna_seq.sh
"""


import sys
import getopt
from random import randint

# BUG: Gets stuck for really large values

def getRandomBase():
  dna_types = ['A', 'G', 'T', 'C']
  return dna_types[randint(0, 3)]

def randomStrand(l):
  strand = ""
  for i in range(l):
    strand += getRandomBase()
  return strand

def distance (s1, s2):
  assert(len(s1) == len(s2))
  d = 0
  for i in range(len(s1)):
    if s1[i] != s2[i]:
      d += 1
  return d

def usage():
    print('$> python3 data_generator.py <required args>\n' + \
        '\t-c <#>\t\tNumber of clusters to generate (min 0)\n' + \
        '\t-p <#>\t\tNumber of strands per cluster (min 0)\n' + \
        '\t-l <#>\t\tLength of DNA strands (min 1)\n')  

def run():
  k = -1
  pk = -1
  length = -1

  # Get command-line arguments
  try:
    optlist, _ = getopt.getopt(sys.argv[1:], 'c:p:l:')
  except getopt.GetoptError as err:
      usage()
      print(str(err))
      sys.exit(2)

  # Check arguments
  for key, val in optlist:
    if key == '-c':
      k = int(val)
    elif key == '-p':
      pk = int(val)
    elif key == '-l':
      length = int(val)

  # check required arguments were inputted  
  if k < 0 or pk < 0 or length < 1:
    usage()
    sys.exit()  

  # The list of strands
  strands = []

  # Generate strands around a centroid, $c times
  for i in range(k):
    # Create a random centroid
    centroid = randomStrand(length);
    strands_c = []
    # Create a set of strands uniformly distributed around this centroid
    while (len(strands_c) < pk):
      strand = randomStrand(length)
      if (distance(strand, centroid) < length/2): # < (length/(k % length))): # Problem if k>=length
        strand += '\n' # newline-terminate the string
        strands_c.append(strand)
    # Add the centroid-distributed strands to the list of all strands
    strands += strands_c


  # Generate & write the DNA strands to a file
  with open("../dna_generator/dna_data.csv", "w+") as f: 
    f.writelines(strands)

  return


run()