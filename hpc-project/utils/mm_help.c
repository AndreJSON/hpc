#include "mm_help.h"

void utils_open_matrix_files(FILE **inFile, FILE **outFile, int argc, char *argv[]) {
	if(argc < 7) {
		printf("Too few arguments");
		exit(-1);
	}
	if(((*inFile = fopen(argv[4], "r")) == NULL) || ((*outFile = fopen(argv[6], "w")) == NULL)) {
		printf("Matrix file does not exist!\n");
		exit(-1);
	}
}

void utils_get_matrix_size(FILE *inFile, int *N) {
	MM_typecode matCode;
	int tmp1;
	if(mm_read_banner(inFile, &matCode) != 0) {
		printf("Problem reading banner!\n");
		exit(-1);
	}
	if(mm_read_mtx_array_size(inFile, N, &tmp1) != 0) {
		printf("Problem getting matrix size!\n");
		exit(-1);
	}
}

void utils_read_matrix(FILE *inFile, int N, double matrix[N][N]) {
	int i, j;
	double tmp3;
	for(i = 0; i < N; i++) {
		for(j = 0; j < N; j++) {
			fscanf(inFile, "%lg\n",&tmp3);
			matrix[i][j] = tmp3;
		}
	}
}

void utils_print_matrix(int N, double matrix[N][N]) {
	int i, j;
	for(i = 0; i < N; i++) {
		for(j = 0; j < N; j++) {
			printf("%.3f  ", matrix[i][j]);
		}
		printf("\n");
	}
}

void utils_print_array(int N, double arr[N]) {
	int i;
	for(i = 0; i < N; i++) {
		printf("%.3f  ", arr[i]);
	}
	printf("\n");
}


void utils_fprint_matrix(int N, double matrix[N][N], FILE *outFile) {
	int i, j;
	fprintf(outFile,"%%%%MatrixMarket matrix array real general\n");
	mm_write_mtx_array_size(outFile, N, N);
	for(i = 0; i < N; i++) {
		for(j = 0; j < N; j++) {
			if(matrix[i][j] >= 0)
				fprintf(outFile," ");
			fprintf(outFile, "%.17f\n", matrix[i][j]);
		}
	}
}

void utils_fprint_array(int N, double arr[N], FILE *outFile) {
	int i, j;
	fprintf(outFile,"%%%%MatrixMarket matrix array real general\n");
	mm_write_mtx_array_size(outFile, N, 1);
	for(i = 0; i < N; i++) {
		if(arr[i] >= 0)
			fprintf(outFile," ");
		fprintf(outFile, "%.17f\n", arr[i]);
	}
}

void utils_get_split(int argc, char *argv[], int *p, int *p1, int *p2) {
	int tmp = argv[2][1];
	if(tmp > 9) { //first digit is only 1 char.
		*p1 = argv[2][0] - '0';
	} else {
		*p1 = (argv[2][0] - '0') * 10 + (argv[2][1] - '0');
	}
	*p2 = (*p)/(*p1);
	if((*p1) * (*p2) != *p) {
		printf("Incorrect p1 or p2!\n");
		exit(-1);
	}
}