#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <stdbool.h>
#include <limits.h>


int l; // Length of DNA strands
int k; // Clusters 
int ts; // Total number of points
int iters = 100000; // Max iterations of KMeans algorithm
char* input_file = NULL; // Input file to read data from

typedef struct strands {
	int cluster; // Cluster strand belongs to
	char* s;
	struct strands* next;
} strand;

typedef struct centroids {
	int cluster; // Cluster this is a centroid of
	int members; // Number of strands in cluster
	char* s;
	struct centroids* next;
} centroid;

// Pointers to list of datapoints
strand* strands_head = NULL;
strand* strands_tail = NULL;

// Pointers to list of centroids
centroid* cent_head = NULL;
centroid* cent_tail = NULL;

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


/* Recursively frees all strands */
void free_strands(strand* sp) {
	if (sp) {
		free_strands(sp->next);
		free(sp->s);
		free(sp);
	}
}


/* Recursively frees all centroids */
void free_centroids(centroid* c) {
	if (c) {
		free_centroids(c->next);
		free(c->s);
		free(c);
	}
}


/* Debugging tool for showing all the points */
void print_strands() {
	strand* sp = strands_head;
	int i = 1;
	while (sp) {
		printf("Strand %i: %s\n", i, sp->s);
		sp = sp->next;
		i++;
	}
}


/* Debugging tool for showing all the centroids */
void print_centroids() {
	centroid* c = cent_head;
	int i = 0;
	while (c) {
		printf("Centroid %i: s: %s c: %i m: %i\n", i, c->s, c->cluster, c->members);
		c = c->next;
		i++;
	}
}


/* Debugging tool for showing all the variables */
void print_all() {
	printf("Clusters: %i\n", k);
	printf("Total strands: %i\n", ts);
	printf("K Means Iterations: %i\n", iters);
	printf("Input file: %s\n", input_file);
	print_strands();
	print_centroids();
}


/* Return the centroid of a given cluster k */
centroid* get_centroid(int k) {
	centroid* tmp = cent_head;
	for (int j = 0; j < k && tmp; j++) {
		tmp = tmp->next;
	}
	return tmp;
}


/* Return the i_th strand from the list of strands */
strand* get_strand(int i) {
	strand* tmp = strands_head;
	for (int j = 0; j < i && tmp; j++) {
		tmp = tmp->next;
	}
	return tmp;
}


/* Adds a DNA strand read in from the file to the linked list */
void add_strand(char* dna_strand) {
	// Remove the newline character at the end of the strand
	dna_strand[strcspn(dna_strand, "\n")] = '\0';

	// Create the strand object
	strand* new_strand = (strand*) malloc(sizeof(strand));
	new_strand->s = (char*) malloc(sizeof(char) * strlen(dna_strand)); 
	memcpy(new_strand->s, dna_strand, strlen(dna_strand));
	new_strand->cluster = 0;
	new_strand->next = NULL;

	// Add strand to list by setting head/tail pointers
	if (strands_head == NULL) {
		strands_head = new_strand;
		strands_tail = new_strand;
	} else {
		strands_tail->next = new_strand;
		strands_tail = new_strand;
	}
}


/* Calculate the distance between two strands by measuring the 
 * the number of changes necessary to go from one strand to the other
 * REQUIRES: strlen(s1) == strlen(s2)
 */
int distance(char* s1, char* s2) {
	int d = 0;

	for (int i = 0; i < strlen(s1); i++) {
		if (s1[i] != s2[i]) {
			d++;
		}
	}

	return d;
}


/* Create k random centroids */
void intialize_centroids() {
	strand* tmp;

	// Assign a random point to be a centroid, k times
	for (int i = 0; i < k; i++) {
		// Move to selected point
		tmp = get_strand(((int) random()) % ts);

		// Set centroid values to strand values
		centroid* new_c = (centroid*) malloc(sizeof(centroid));
		new_c->cluster = i;
		if (i == 0) { // All strands initially assigned to cluster 0
			new_c->members = ts;
		} else {
			new_c->members = 0;
		}
		new_c->next = NULL;
		new_c->s = (char*) malloc(sizeof(char) * (l+1));
		strcpy(new_c->s, tmp->s);
		
		// Update list of centroids
		if (cent_head == NULL) {
			cent_head = new_c;
			cent_tail = new_c;
		} else {
			cent_tail->next = new_c;
			cent_tail = new_c;
		}
	}
}


/* Returns the max of two numbers */
int max(int x, int y) {
	if (x > y) {
		return x;
	} else {
		return y;
	}
}


/* Update the value of a single centroid (repeat this k times to update all of them) */
void update_centroid(int c) {
	strand* sp;
	int a_count, c_count, g_count, t_count, highest_count;
	char* cent_new = (char*) malloc((l+1) * sizeof(char)); // The new centroid strand
	int table[4][l]; // Table containing frequency values
	memset(table, 0, 4 * l * sizeof(int));
	centroid* cent = get_centroid(c);

	// Get table values 
	for (int i = 0; i < ts; i++) {
		sp = get_strand(i);
		if (sp->cluster == c) {
			for (int j = 0; j < l; j++) {
				// Update table[b][j] according to the base
				switch (sp->s[j]) {
					case 'A':
						table[0][j]++;
						break;
					case 'C':
						table[1][j]++;
						break;
					case 'G':
						table[2][j]++;
						break;
					case 'T':
						table[3][j]++;
						break;
				}
			}
		}
	}

	// Create the new centroid by finding the highest frequency base of each 
	// index and setting that base as cent_new[index]
	for (int j = 0; j < l; j++) {
		a_count = table[0][j];
		c_count = table[1][j];
		g_count = table[2][j];
		t_count = table[3][j];
		highest_count = max(a_count, max(c_count, max(g_count, t_count)));
		if (highest_count == a_count) {
			cent_new[j] = 'A';
		} else if (highest_count == c_count) {
			cent_new[j] = 'C';
		} else if (highest_count == g_count) {
			cent_new[j] = 'G';
		} else {
			cent_new[j] = 'T';
		}
	}
	cent_new[l] = '\0';
	free(cent->s); // free old centroid pointer
	cent->s = cent_new; /* New centroid assignment */
}


/* Update cluster membership of the points then change the centroids.
 * Returns if any of the clusters changed
 */
bool update_clusters() {
	bool clusters_changed = false;

	// For each point, calculate its distance to each cluster centroid
	for (int i = 0; i < ts; i++) {
		// Get the point's current values
		strand* p = get_strand(i);
		int cluster_old = p->cluster;
		centroid* cent_old = get_centroid(cluster_old);
		int d_old = distance(p->s, cent_old->s);
		
		// Compare point's distance to each centroid different than its own cluster's
		for (int cluster_new = 0; cluster_new < k; cluster_new++) {
			if (cluster_new != cluster_old) {
				centroid* cent_new = get_centroid(cluster_new);
				int d_new = distance(p->s, cent_new->s);

				// If the new distance is less than the old, then update the strand's cluster
				if (d_new < d_old) { 
					p->cluster = cluster_new; /* Membership update */
					cent_old->members--;
					cent_new->members++;
					cluster_old = cluster_new;
					cent_old = cent_new;
					d_old = d_new;
					clusters_changed = true;
				}
			}
		}
	}

	// Calculate new centroids
	for (int i = 0; i < k; i++) 
		update_centroid(i);

	return clusters_changed;
}


/* Runs the KMeans algorithm on the points and returns the data in handout-specified form */
char* k_means() {
	char* data = (char*) malloc(k * 256); // Really do need to malloc (vars get freed after block)
	int size = 0;

	// Set random initial centroids
	intialize_centroids();

	// Cluster until we reach desired iterations or the clusters stop changing
	for (int i = 0; i < iters && update_clusters(); i++) {}

	// Create the output file data
	for (int i = 0; i < k; i++) {
		centroid* cent = get_centroid(i);
		size += sprintf(data + size, "%s, %i\n", cent->s, cent->members); 
	}
	sprintf(data + size, ""); // Null-terminate data

	return data;
}



/* 
 * MAIN
 */
int main(int argc, char** argv) { // argc is 1 when no params passed in
	// Read in commandline arguments
	if (handle_args(argc, argv)) return 1; // Exit on error

	// Read the input file to get the strands
	char* buf = (char*) malloc(l+2); // TODO: Check if strands equal input file (also for points.c)
	FILE* fp = fopen(input_file, "r");
	while (fgets(buf, l+2, fp) != NULL) {add_strand(buf);}
	fclose(fp);
	free(buf);

	// Run kmeans algorithm and produce results
	char* cluster_data = k_means();

	// Write results to an output file
	fp = fopen("dna_seq_results.csv", "w+");
	fprintf(fp, "%s", cluster_data);
	fclose(fp);

	// Free malloc'ed structures and return
	free(cluster_data);
	free_strands(strands_head);
	free_centroids(cent_head);
	return 0;
}
