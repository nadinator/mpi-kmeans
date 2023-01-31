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

int k; // Clusters 
int tp; // Total number of points
int iters = 100000; // Max iterations of KMeans algorithm
char* input_file = NULL; // Input file to read data from
int world_size;
int rank;

typedef struct Point {
	int cluster; // Cluster point belongs to. 0 <= cluster < k
	double x;
	double y;
} point;

typedef struct Centroid {
	int members; // Number of points in cluster
	double x;
	double y;
} centroid;

point* points = NULL; // List of points
centroid* centroids = NULL; // List of centroids


/* Check the argument flags. Return 1 on error, 0 otherwise */
int check_flags(int k_FLAG, int total_FLAG, int iters_FLAG, int file_FLAG) {
	if (!k_FLAG) {
		fprintf(stderr, "Error: option -c must be included\n");
		return 1;
	} else if (k < 1) { // And tp?
		fprintf(stderr, "Error: clusters cannot be less than one\n");
		return 1;
	} else if (!total_FLAG) {
		fprintf(stderr, "Error: option -t must be included\n");
		return 1;
	} else if (tp < 1) {
		fprintf(stderr, "Error: total points cannot be less than one\n");
		return 1;
	} else if (!file_FLAG) {
		fprintf(stderr, "Error: option -i must be included\n");
		return 1;
	} else if (iters_FLAG && iters < 1) {
		fprintf(stderr, "Error: iterations cannot be less than one\n");
		return 1;
	} else if (tp < world_size) {
		fprintf(stderr, "Error: Can't have fewer points than processes\n");
		return 1;
	}

	return 0;
}


/* Handle input arguments. Return 1 on error, 0 otherwise */
int handle_args(int argc, char** argv) {
	int o;
	int file_FLAG = 0;
	int k_FLAG = 0;
	int total_FLAG = 0;
	int iters_FLAG = 0;

	while ((o = getopt(argc, argv, "c:t:i:n:")) != -1) {
		switch (o) {
			case 'c':
				if (optarg != NULL && optarg[0] != '-') {
					k = atoi(optarg);
					k_FLAG = 1;
				}
				break;
			case 't':
				if (optarg != NULL && optarg[0] != '-') {
					tp = atoi(optarg);
					total_FLAG = 1;
				}
				break;
			case 'i':
				if (optarg != NULL && optarg[0] != '-') {
					input_file = optarg;
					file_FLAG = 1;
				}
				break;
			case 'n':
				if (optarg != NULL && optarg[0] != '-') {
					iters = atoi(optarg); 
					iters_FLAG = 1;
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

	return check_flags(k_FLAG, total_FLAG, iters_FLAG, file_FLAG);
}


/* Debugging tool for showing all the points */
void print_points() {
	point p;
	int i;

	printf("--------- Points --------\n");
	for (i = 0; i < tp; i++) {
		p = points[i];
		printf("Point %i: x:%.15f y:%.15f cluster: %i\n", i+1, p.x, p.y, p.cluster);	
	}
	printf("-------- End Points -------\n");
}


/* Debugging tool for showing all the centroids */
void print_centroids() {
	centroid c;
	int i;

	for (i = 0; i < k; i++) {
		c = centroids[i];
		printf("Centroid %i: x:%.15f y:%.15f members: %i\n", i, c.x, c.y, c.members);
	}
}


/* Calculate the Euclidean distance between two points */
double distance(double x1, double y1, double x2, double y2) {
	return sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
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


/* Update the value of centroid[c] (repeat this k times to update all of them) 
 * Runs only in root.
 */
void update_centroid_values() {
	int i;
	double* sums_x = calloc(k, sizeof(double));
	double* sums_y = calloc(k, sizeof(double));

	// Sum over members' values
	for (i = 0; i < tp; i++) {
		int cluster = points[i].cluster;
		sums_x[cluster] += points[i].x;
		sums_y[cluster] += points[i].y;
	}

	// New values are the mean of all members
  for (i = 0; i < k; i++) {
		int members = centroids[i].members;
		if (members != 0) {
			centroids[i].x = sums_x[i] / members;
			centroids[i].y = sums_y[i] / members;
		}
	}

	free(sums_x);
	free(sums_y);
}


/* Update cluster membership of the points by calculating distance to each centroid.
 * Return whether any of the clusters changed. 
 * Runs only in root.
 */
bool cluster() {
	bool clusters_changed = false;
	point p;
	centroid centroid_old, centroid_new;
	int cluster_old, cluster_new;
	double d_old, d_new;
	int i;

	// For each point, calculate its distance to each cluster centroid
	for (i = 0; i < tp; i++) {
		// Get the point's current values
		p = points[i];
		cluster_old = p.cluster;
		centroid_old = centroids[cluster_old]; 
		d_old = distance(p.x, p.y, centroid_old.x, centroid_old.y);
		
		// Compare point's distance to each centroid different than its own cluster's
		for (cluster_new = 0; cluster_new < k; cluster_new++) {
			if (cluster_new != cluster_old) {
				centroid_new = centroids[cluster_new];
				d_new = distance(p.x, p.y, centroid_new.x, centroid_new.y);
				
				// If the new distance is less than the old, then update the point's cluster
				if (d_new < d_old) {
					points[i].cluster = cluster_new; // Membership update
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


// Creates the POINT datatype for passing around points
void create_point_dt(MPI_Datatype* dt) {
	MPI_Datatype tmp_POINT;
	MPI_Aint lb, extent;
	
	const int count = 3; // p.cluster, p.x, p.y
	int block_lengths[count] = {1, 1, 1};
	MPI_Aint elem_displacements[count] = {offsetof(point, cluster), offsetof(point, x), offsetof(point, y)};
	MPI_Datatype types[count] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE};
	
	MPI_Type_create_struct(count, block_lengths, elem_displacements, types, &tmp_POINT);
	MPI_Type_get_extent(tmp_POINT, &lb, &extent);
	MPI_Type_create_resized(tmp_POINT, lb, extent, dt);
	MPI_Type_commit(dt);
}


// Iniitialize points array in root process
void init_points() {
	int i = 0;
	double x = 0.0;
	double y = 0.0;

	// Allocate points memory
	points = (point*) malloc(sizeof(point) * tp);
	if (!points) {
		fprintf(stderr, "Couldn't allocate memory for points\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// Read input data
	FILE* fp = fopen(input_file, "r");
	for (i = 0; i < tp; i++) {	
		if (fscanf(fp, "%lf,%lf", &x, &y) != 2) {
			fprintf(stderr, "Reading failure.\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		point p;
		p.x = x;
		p.y = y;
		p.cluster = 0;
		points[i] = p;
	}
	fclose(fp);
}


// Initialize centroids array in root process
void init_centroids() {
	int i;
	point tmp_point;
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
		tmp_point = points[((int) rand()) % tp];
		tmp_centroid.x = tmp_point.x;
		tmp_centroid.y = tmp_point.y;
		if (i == 0) { // All points initially belong to cluster 0
			tmp_centroid.members = tp;
		} else {
			tmp_centroid.members = 0;
		}
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
		size += sprintf(results + size, "%f, %f, %i\n", 
										centroids[i].x, centroids[i].y, centroids[i].members); 
	}

	// Write results to file
	FILE* fp = fopen("points_par_results.csv", "w+");
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
	printf("See output file points_par_results for handout-formatted results.\n");
	printf("===================================================================\n");
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
		init_points();
		// Initialize centroids array by picking random points
		init_centroids();
		// Do initial clustering
		clusters_changed = cluster();
  }

	// Create the datatype for sending points around
	MPI_Datatype POINT;
	create_point_dt(&POINT);
	
	// Broadcast k, tp, iters, and clusters_changed
	MPI_Bcast(&k, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	MPI_Bcast(&tp, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	MPI_Bcast(&iters, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	MPI_Bcast(&clusters_changed, 1, MPI_C_BOOL, ROOT, MPI_COMM_WORLD);

	// Set up variables for MPI_Scatterv
	int* send_counts = NULL; 
	int* displacements = NULL;
	int num_local_points = tp / world_size; // Number of points given to each process
	int num_rm = tp % world_size; // Possible remainder points
	int displacement_num = num_local_points; // Number of points to displace the sendbuf by
	// Populate send_counts
	send_counts = (int*) malloc(sizeof(int) * world_size);
	for (i = 0; i < world_size; i++) {send_counts[i] = num_local_points;} 
	if (num_rm != 0) {
		// Give remainder points to last process
		send_counts[world_size-1] = num_local_points + num_rm;
		if (rank == world_size -1) 
			num_local_points = num_local_points + num_rm;
	}
	// Populate displacements
	displacements = (int*) malloc(sizeof(int) * world_size);
	for (i = 0; i < world_size; i++) displacements[i] = i * displacement_num; 
	
	// Start KMeans iterations
	int iter = 0;
	while (iter < iters && clusters_changed) {
		// Receive points for processing cluster membership
		point* local_points = (point*) calloc(send_counts[rank], sizeof(point));
		MPI_Scatterv(points, send_counts, displacements, POINT, 
							   local_points, send_counts[rank], POINT, 
								 ROOT, MPI_COMM_WORLD);

		// Find how many points in the subset belong to each cluster
		int* local_members = (int*) calloc(k, sizeof(int)); // local_members[i] is cluster i's membership count
		for (i = 0; i < num_local_points; i++) {
			local_members[local_points[i].cluster]++; // Populate list with centroid memberships
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
		free(local_points);
		free(local_members);
		if (c_all != NULL) free(c_all);
	}

	// Write cluster data to an output file and print results
	if (rank == ROOT) {	
		write_output_file();
		print_results(iter, start_time);
	}
	
	// Free pointers
	if (points != NULL) free(points); 
	if (centroids != NULL) free(centroids);
	free(send_counts);
	free(displacements);
	
	// Finalize MPI environment and return
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
