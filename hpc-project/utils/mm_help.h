#include <stdio.h>
#include <stdlib.h>
#include "../utils/mmio.h"


void utils_open_matrix_files(FILE **inFile, FILE **outFile, int argc, char *argv[]);
void utils_get_matrix_size(FILE *inFile, int *N);
void utils_read_matrix(FILE *inFile, int N, double matrix[N][N]);
void utils_print_matrix(int N, double matrix[N][N]);
void utils_print_array(int N, double arr[N]);
void utils_fprint_matrix(int N, double matrix[N][N], FILE *outFile);
void utils_fprint_array(int N, double arr[N], FILE *outFile);
void utils_get_split(int argc, char *argv[], int *p, int *p1, int *p2);