#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <stdbool.h>
#include <limits.h>

int k; // Clusters 
int tp; // Total number of points
int iters = 100000; // Max iterations of KMeans algorithm
char* input_file = NULL; // Input file to read data from


typedef struct points {
	int cluster; // Cluster point belongs to
	double x;
	double y;
	struct points* next;
} point;

typedef struct centroids {
	int cluster; // Cluster this is a centroid of
	int members; // Number of points in cluster
	double x;
	double y;
	struct centroids* next;
} centroid;

// Pointers to list of datapoints
point* points_head = NULL;
point* points_tail = NULL;

// Pointers to list of centroids
centroid* cent_head = NULL;
centroid* cent_tail = NULL;


/* Recursively frees all points */
void free_points(point* p) {
	if (p) {
		free_points(p->next);
		free(p);
	}
}


/* Recursively frees all centroids */
void free_centroids(centroid* c) {
	if (c) {
		free_centroids(c->next);
		free(c);
	}
}


/* Debugging tool for showing all the points */
void print_points() {
	point* p = points_head;
	int i = 1;
	while (p) {
		printf("Point %i: x:%.15f y:%.15f\n", i, p->x, p->y);
		p = p->next;
		i++;
	}
}


/* Debugging tool for showing all the centroids */
void print_centroids() {
	centroid* c = cent_head;
	int i = 1;
	while (c) {
		printf("Centroid %i: x:%.15f y:%.15f\n", i, c->x, c->y);
		c = c->next;
		i++;
	}
}


/* Return the centroid of a given cluster i */
centroid* get_centroid(int k) {
	centroid* tmp = cent_head;
	for (int j = 0; j < k && tmp; j++) {
		tmp = tmp->next;
	}
	return tmp;
}


/* Return the ith point from the list of points */
point* get_point(int i) {
	point* tmp = points_head;
	for (int j = 0; j < i && tmp; j++) {
		tmp = tmp->next;
	}
	return tmp;
}


/* Converts a line of the input file to a list containing two doubles */
double* string_to_point(char* str) {
	double* points = (double*) malloc(sizeof(double) * 2);
	double pt1, pt2;
	int length = strlen(str);
	int i = 0;
	int j = 0;
	char* buf1 = (char*) malloc(256);
	char* buf2 = (char*) malloc(256);

	// Read in first point
	while (i < length && str[i] != ',') {
		buf1[i] = str[i];
		i++;
	}
	
	// Convert buf1 to double for pt1
	buf1[i] = '\0';
	pt1 = strtod(buf1, NULL);
	points[0] = pt1;
	free(buf1);

	// Read in second point
	i++;
	while (i < length && str[i] != '\n') {
		buf2[j] = str[i];
		i++;
		j++;
	}

	// Convert buf2 to double for pt2
	buf2[j] = '\0';
	pt2 = strtod(buf2, NULL);
	points[1] = pt2;
	free(buf2);

	return points;
}


/* Adds a point read in from the file to the linked list */
void add_point(double* points) {
	point* newPoint = (point*) malloc(sizeof(point));
	newPoint->x = points[0];
	newPoint->y = points[1];
	free(points); 
	newPoint->cluster = 0;
	newPoint->next = NULL;

	if (points_head == NULL) {
		points_head = newPoint;
		points_tail = newPoint;
	} else {
		points_tail->next = newPoint;
		points_tail = newPoint;
	}
}


/* Calculate the Euclidean distance between two points */
double distance(double x1, double y1, double x2, double y2) {
	return sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
}


/* Create k random centroids */
void intialize_centroids() {
	point* tmp;

	// Assign a random point to be a centroid, k times
	for (int i = 0; i < k; i++) {
		// Move to selected point
		tmp = get_point(((int) random()) % tp);

		// Set centroid values to point values
		centroid* new_c = (centroid*) malloc(sizeof(centroid));
		new_c->cluster = i;
		if (i == 0) { // All points initially assigned to cluster 0
			new_c->members = tp;
		} else {
			new_c->members = 0;
		}
		new_c->x = tmp->x;
		new_c->y = tmp->y;
		new_c->next = NULL;
		
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


/* Update the value of a single centroid (repeat this k times to update all of them) 
 * Returns whether a centroid's value has changed from the previous iteration
 */ 
void update_centroid(int c) {
	double sum_x = 0.0; 
	double sum_y = 0.0;
	centroid* cent = get_centroid(c);

	// Sum over members' values
	for (int i = 0; i < tp; i++) {
		point* p = get_point(i);
		if (p->cluster == c) {
			sum_x += p->x;
			sum_y += p->y;
		}
	}

	// New values are the mean of all members
	if (cent->members != 0) {
		cent->x = sum_x / cent->members;
		cent->y = sum_y / cent->members;
	}
}


/* Update cluster membership of the points and change the centroids */
bool update_clusters() {
	bool clusters_changed = false;

	// For each point, calculate its distance to each cluster centroid
	for (int i = 0; i < tp; i++) {
		// Get the point's current values
		point* p = get_point(i);
		int cluster_old = p->cluster;
		centroid* cent_old = get_centroid(cluster_old);
		double d_old = distance(p->x, p->y, cent_old->x, cent_old->y);
		
		// Compare point's distance to each centroid different than its own cluster's
		for (int cluster_new = 0; cluster_new < k; cluster_new++) {
			if (cluster_new != cluster_old) {
				centroid* cent_new = get_centroid(cluster_new);
				double d_new = distance(p->x, p->y, cent_new->x, cent_new->y);
				
				// If the new distance is less than the old then update the point's cluster
				if (d_new < d_old) { 
					p->cluster = cluster_new; // Membership update
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
	for (int i = 0; i < k; i++) update_centroid(i);

	return clusters_changed;
}


/* Runs the KMeans algorithm on the points and returns the data in handout-specified form */
char* k_means() {
	char* data = (char*) malloc(k * 256);
	int size = 0;
	
	// Set random initial centroids
	intialize_centroids();
	
	// Cluster until we reach desired iterations or the clusters stop changing
	for (int i = 0; i < iters && update_clusters(); i++) {}

	// Create the output file data
	for (int i = 0; i < k; i++) {
		centroid* cent = get_centroid(i);
		size += sprintf(data + size, "%f, %f, %i\n", cent->x, cent->y, cent->members); 
	}
	sprintf(data + size, "");

	return data;
}


/* Check the argument flags. Return 1 on error, 0 otherwise */
int check_flags(int k_FLAG, int total_FLAG, int iters_FLAG, int file_FLAG) {
	if (!k_FLAG) {
		fprintf(stderr, "Error: option -c must be included\n");
		return 1;
	} else if (k < 0) { // And tp?
		fprintf(stderr, "Error: clusters cannot be less than zero\n");
		return 1;
	} else if (!total_FLAG) {
		fprintf(stderr, "Error: option -t must be included\n");
		return 1;
	} else if (tp < 0) {
		fprintf(stderr, "Error: total points cannot be less than zero\n");
		return 1;
	} else if (!file_FLAG) {
		fprintf(stderr, "Error: option -i must be included\n");
		return 1;
	} else if (iters_FLAG && iters < 0) {
		fprintf(stderr, "Error: iterations cannot be less than zero\n");
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
					iters = atoi(optarg); // 1 if unitialized (for some reason)
					iters_FLAG = 1;
				}
				break;
			case '?':
				printf("here");
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


/* 
 * MAIN
 */
int main(int argc, char** argv) { // argc is 1 when no params passed in
	if (handle_args(argc, argv)) return 1; // Exit on error

	// Read the input file to get the data points
	char* buf = (char*) malloc(256);
	FILE* fp = fopen(input_file, "r");
	while (fgets(buf, 256, fp) != NULL) add_point(string_to_point(buf));
	fclose(fp);
	free(buf);

	// Now I have a list of points, and the real Kmeans algorithm begins.
	char* cluster_data = k_means();

	// Write cluster data to an output file
	fp = fopen("points_seq_results.csv", "w+");
	fprintf(fp, "%s", cluster_data);
	fclose(fp);

	free(cluster_data);
	free_points(points_head);
	free_centroids(cent_head);
	return 0;
}
