#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <limits.h>
#include <mpi.h>
#include <stddef.h>

#define ROOT 0

int l; // Length of DNA strands
int k; // Clusters 
int ts; // Total number of points
int iters = 100000; // Max iterations of KMeans algorithm
char* input_file = NULL; // Input file to read data from
int rank, world_size;

typedef struct strands {
	int cluster; // Cluster strand belongs to
	char* s;
} strand;

typedef struct centroids {
	int members; // Number of strands in cluster
	char* s;
} centroid;

strand* strands = NULL;
centroid* centroids = NULL;


/* Checks flags of input arguments, returning 1 on error and 0 otherwise */
int check_flags(int k_FLAG, int total_FLAG, int iters_FLAG, int file_FLAG, int l_FLAG) {
	if (!k_FLAG) {
		fprintf(stderr, "Error: option -c must be included\n");
		return 1;
	} else if (k < 1) {
		fprintf(stderr, "Error: clusters cannot be less than one.\n");
		return 1;
	} else if (!total_FLAG) {
		fprintf(stderr, "Error: option -t must be included\n");
		return 1;
	} else if (ts < 1) {
		fprintf(stderr, "Error: total points cannot be less than one\n");
		return 1;
	} else if (!file_FLAG) {
		fprintf(stderr, "Error: option -i must be included\n");
		return 1;
	} else if (!l_FLAG) {
		fprintf(stderr, "Error: option -l must be included\n");
		return 1;
	} else if (l < 1) {
		fprintf(stderr, "Error: length cannot be less than one\n");
		return 1;
	} else if (iters_FLAG && iters < 1) {
		fprintf(stderr, "Error: iterations cannot be less than one\n");
		return 1;
	} else if (ts < world_size) {
		fprintf(stderr, "Error: Can't have fewer strands than processes\n");
		return 1;
	}

	return 0;
}


/* Handle the input arguments and assign global variables 
 * Returns 1 on error, 0 otherwise 
 */
int handle_args(int argc, char** argv) {
	int o;
	int file_FLAG = 0;
	int k_FLAG = 0;
	int total_FLAG = 0;
	int iters_FLAG = 0;
	int l_FLAG = 0;

	// Get option values
	while ((o = getopt(argc, argv, "c:t:i:n:l:")) != -1) {
		switch (o) {
			case 'c':
				if (optarg && optarg[0] != '-') {
					k = atoi(optarg);
					k_FLAG = 1;
				}
				break;
			case 't':
				if (optarg && optarg[0] != '-') {
					ts = atoi(optarg);
					total_FLAG = 1;
				}
				break;
			case 'i':
				if (optarg && optarg[0] != '-') {
					input_file = optarg;
					file_FLAG = 1;
				}
				break;
			case 'n':
				if (optarg && optarg[0] != '-') {
					iters = atoi(optarg); // 1 if unitialized (for some reason)
					iters_FLAG = 1;
				}
				break;
			case 'l':
				if (optarg && optarg[0] != '-') {
					l = atoi(optarg);
					l_FLAG = 1;
				}
				break;
			case '?':
				if (optarg == NULL && strrchr("ctinl", optopt))
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
		}
	}

	return check_flags(k_FLAG, total_FLAG, iters_FLAG, file_FLAG, l_FLAG);
}


// Frees the strands pointer and all strand strings
void free_strands() {
  int i;
  for (i = 0; i < ts; i++) free(strands[i].s);
  free(strands);
}


// Frees the centroids pointer and all centroid strings
void free_centroids() {
  int i;
  for (i = 0; i < k; i++) free(centroids[i].s);
  free(centroids);
}


/* Debugging tool for showing all the points */
void print_strands() {
	int i;
  strand s;
  for (i = 0; i < ts; i++) {
    s = strands[i];
    printf("Strand %i: %s cluster: %i\n", i+1, s.s, s.cluster);
  }
}


/* Debugging tool for showing all the centroids */
void print_centroids() {
	centroid c;
	int i;
	for (i = 0; i < k; i++) {
		c = centroids[i];
		printf("Centroid %i: %s members: %i\n", i, c.s, c.members);
	}
}


/* Calculate the distance between two strands by measuring the 
 * the number of changes necessary to go from one strand to the other
 * REQUIRES: strlen(s1) == strlen(s2)
 */
int distance(char* s1, char* s2) {
  int i;
	int d = 0;

	for (i = 0; i < strlen(s1); i++) {
		if (s1[i] != s2[i]) {
			d++;
		}
	}

	return d;
}


/* Returns the max of two numbers */
int max(int x, int y) {
	if (x > y) {
		return x;
	} else {
		return y;
	}
}


// Creates the STRAND datatype for passing around strands
void create_strand_dt(MPI_Datatype* dt) {
	MPI_Datatype tmp_STRAND;
	MPI_Aint lb, extent;
	
	const int count = 2; // s.cluster, s.s
	int block_lengths[count] = {1, l+1};
	MPI_Aint elem_displacements[count] = {offsetof(strand, cluster), offsetof(strand, s)};
	MPI_Datatype types[count] = {MPI_INT, MPI_CHAR};
	
	MPI_Type_create_struct(count, block_lengths, elem_displacements, types, &tmp_STRAND);
	MPI_Type_get_extent(tmp_STRAND, &lb, &extent);
	MPI_Type_create_resized(tmp_STRAND, lb, extent, dt);
	MPI_Type_commit(dt);
}


// Iniitialize points array in root process
void init_strands() {
  int i;
  strand str;
	char strand_buf[l+1];

	// Allocate strands memory
	strands = (strand*) malloc(sizeof(strand) * ts);
	if (!strands) {
		fprintf(stderr, "Couldn't allocate memory for strands\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// Read input file and populate strands list with its data
	FILE* fp = fopen(input_file, "r");
	for (i = 0; i < ts; i++) {	
    // Read a single line 
		if (fscanf(fp, "%s\n", strand_buf) != 1) {
			fprintf(stderr, "Reading failure.\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

    // Set new strand's values
    str.s = (char*) malloc(sizeof(char) * (l+1));
    if (!str.s) {
		  fprintf(stderr, "Couldn't allocate memory for a strand's string\n");
		  MPI_Abort(MPI_COMM_WORLD, 1);
	  }
    strcpy(str.s, strand_buf);
    str.cluster = 0; // All strands initially belong to cluster 0

    // Add strand to list
		strands[i] = str;
	}

	fclose(fp);
}


// Initialize centroids array in root process
void init_centroids() {
	int i;
	strand tmp_strand;
	centroid tmp_centroid;
	srand(time(NULL));

	// Allocate memory for centroids array
	centroids = (centroid*) malloc(sizeof(centroid) * k);
	if (!centroids) {
		fprintf(stderr, "Couldn't allocate memory for centroids.\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// Set random initial centroids
	for (i = 0; i < k; i++) {		
    // Pick a random strand 
		tmp_strand = strands[((int) rand()) % ts];

    // Set new centroid's values
    tmp_centroid.s = (char*) malloc(sizeof(char) * (l+1));
    if (!tmp_centroid.s) {
		  fprintf(stderr, "Couldn't allocate memory for centroid's string.\n");
		  MPI_Abort(MPI_COMM_WORLD, 1);
	  }
		strcpy(tmp_centroid.s, tmp_strand.s);
		if (i == 0) { // All strands initially belong to cluster 0
			tmp_centroid.members = ts;
		} else {
			tmp_centroid.members = 0;
		}

    // Add centroid to list
		centroids[i] = tmp_centroid;
	}
}


// Writes results to a file called points_par_results.csv
void write_output_file() {
	int i;
	int size = 0;
	char* results = (char*) malloc(k * 256); // k lines of output

	// Create the output file data
	for (i = 0; i < k; i++) {
		size += sprintf(results + size, "%s, %i\n", 
										centroids[i].s, centroids[i].members); 
	}

	// Write results to file
	FILE* fp = fopen("dna_par_results.csv", "w+");
	fprintf(fp, "%s", results);
	
	fclose(fp);
	free(results);
}

// Prints results of running the KMeans algorithm
void print_results(int iter, double start_time) {
	iter = iter == 0 ? 0 : iter-1;
	double total_time = MPI_Wtime() - start_time;

	printf("========================== Final Results ==========================\n");
	printf("Iterations ran: %d\n", iter);
	printf("Time taken: %f seconds\n", total_time);
	printf("\n");
	print_centroids();
	printf("\n");
	printf("See output file dna_par_results for handout-formatted results.\n");
	printf("===================================================================\n");
}


// Sum over gathered membership numbers to find centroid[k].members for all k
// Runs only in root.
void collect_centroid_members(int* c_all, int world_size) {
	int centroid_local_members, total_members;
	int i, j;

	for (i = 0; i < k; i++) { // Pick a centroid
		total_members = 0;
		for (j = 0; j < world_size; j++) { // Jump over other processes
			centroid_local_members = c_all[j * k + i]; // Skip to centroid's value
			total_members += centroid_local_members;
		}
		centroids[i].members = total_members;
	}
}


/* Update the centroids' string values.
 * Runs only in root.
 */
void update_centroid_values() {
	strand str;
  int i, j;
  int cluster;
	int a_count, c_count, g_count, t_count, highest_count;
	char centroid_s_new[l+1];
	int table[k][4][l]; // Table containing frequency values for each centroid
	memset(table, 0, k * 4 * l * sizeof(int));

	// Create table, where table[k] is centroid k's frequency table
	for (i = 0; i < ts; i++) {
    str = strands[i];
		cluster = str.cluster;
    for (j = 0; j < l; j++) {
      // Update table[cluster][base][j]
      switch (str.s[j]) {
        case 'A':
          table[cluster][0][j]++;
          break;
        case 'C':
          table[cluster][1][j]++;
          break;
        case 'G':
          table[cluster][2][j]++;
          break;
        case 'T':
          table[cluster][3][j]++;
          break;
      }
    }
	}

	// Create the new centroid string by finding the highest frequency base  
	// of each string index and setting that base as centroid_s_new[index]
  for (i = 0; i < k; i++) { // Repeat for each centroid
    for (j = 0; j < l; j++) {
      a_count = table[i][0][j];
      c_count = table[i][1][j];
      g_count = table[i][2][j];
      t_count = table[i][3][j];
      highest_count = max(a_count, max(c_count, max(g_count, t_count)));
      if (highest_count == a_count) {
        centroid_s_new[j] = 'A';
      } else if (highest_count == c_count) {
        centroid_s_new[j] = 'C';
      } else if (highest_count == g_count) {
        centroid_s_new[j] = 'G';
      } else {
        centroid_s_new[j] = 'T';
      }
    }
    centroid_s_new[l] = '\0';
    strcpy(centroids[i].s, centroid_s_new); 
  }
}


/* Update cluster membership of the points by calculating distance to each centroid.
 * Return whether any of the clusters changed. 
 * Runs only in root.
 */
bool cluster() {
	bool clusters_changed = false;
	strand str;
	centroid centroid_old, centroid_new;
	int cluster_old, cluster_new;
	double d_old, d_new;
	int i;

	// For each point, calculate its distance to each cluster centroid
	for (i = 0; i < ts; i++) {
		// Get the point's current values
		str = strands[i];
		cluster_old = str.cluster;
		centroid_old = centroids[cluster_old]; 
		d_old = distance(str.s, centroid_old.s);
		
		// Compare point's distance to each centroid different than its own cluster's
		for (cluster_new = 0; cluster_new < k; cluster_new++) {
			if (cluster_new != cluster_old) {
				centroid_new = centroids[cluster_new];
				d_new = distance(str.s, centroid_new.s);
				
				// If the new distance is less than the old, then update the point's cluster
				if (d_new < d_old) {
					strands[i].cluster = cluster_new; // Membership update
					centroids[cluster_old].members--;
					centroids[cluster_new].members++;
					cluster_old = cluster_new;
					centroid_old = centroid_new;
					d_old = d_new;
					clusters_changed = true;
				}
			}
		}
	}
	
	// Calculate centroids based on new members
	if (clusters_changed) update_centroid_values();

	return clusters_changed;
}

/* 
 * MAIN
 */
int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  bool clusters_changed;
	int i;
	double start_time; 

  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Master node initializes everything
  if (rank == ROOT) {
		// Start timer 
		start_time = MPI_Wtime();
    // Read commandline input and populate global variables
    if (handle_args(argc, argv)) MPI_Abort(MPI_COMM_WORLD, 0);; // Exit on error
    // Initialize array of points by reading from input file
		init_strands();
    print_strands();
		// Initialize centroids array by picking random points
		init_centroids();
    print_centroids();
		// Do initial clustering
		clusters_changed = cluster();
  }

  // Create the STRAND datatype for message-passing
  MPI_Datatype STRAND;
  create_strand_dt(&STRAND);

  // Broadcast k, ts, iters, clusters_changed, and strand length
	MPI_Bcast(&k, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&l, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	MPI_Bcast(&ts, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	MPI_Bcast(&iters, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	MPI_Bcast(&clusters_changed, 1, MPI_C_BOOL, ROOT, MPI_COMM_WORLD);

  // Set up variables for MPI_Scatter
	int* send_counts = NULL; 
	int* displacements = NULL;
	int num_local_strands = ts / world_size; // Number of strands given to each process

  // Start KMeans iterations
	int iter = 0;
	while (iter < iters && clusters_changed) {
		// Receive points for processing cluster membership
		strand* local_strands = (strand*) calloc(num_local_strands, sizeof(strand));
		MPI_Scatter(strands, num_local_strands * sizeof(strand), MPI_BYTE, 
								local_strands, num_local_strands * sizeof(strand), MPI_BYTE, 
								ROOT, MPI_COMM_WORLD);

		// Find how many points in the subset belong to each cluster
		int* local_members = (int*) calloc(k, sizeof(int)); // local_members[i] is cluster i's membership count
		for (i = 0; i < num_local_strands; i++) {
			local_members[local_strands[i].cluster]++; // Populate list with centroid memberships
		}
		
		// Gather list of members from all processes to the root
		int* c_all = NULL; // List of centroid members from each process
		if (rank == ROOT) c_all = (int*) calloc(k * world_size, sizeof(int));
		MPI_Gather(local_members, k, MPI_INT, c_all, k, MPI_INT, ROOT, MPI_COMM_WORLD);
		
		// Update centroid members in the root and run clustering algorithm
		if (rank == ROOT) {
			collect_centroid_members(c_all, world_size);
			clusters_changed = cluster();
		}

		iter++;
		MPI_Bcast(&clusters_changed, 1, MPI_C_BOOL, ROOT, MPI_COMM_WORLD);
		free(local_strands);
		free(local_members);
		if (c_all != NULL) free(c_all);
	}

  // Write cluster data to an output file and print results
	if (rank == ROOT) {	
		write_output_file();
		print_results(iter, start_time);
	}

  // Free pointers
	free(send_counts);
	free(displacements);
	if (rank == ROOT) {
		free_strands(); 
		free_centroids();
	}
	
  // Finalize MPI environment and return
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
  return 0;
}
