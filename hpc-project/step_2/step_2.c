#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "../utils/mm_help.h"

int localPivotIndex(int N, double matrix[N][N], int k, int rank, int np) { //index of locally largest pivot element.
	int i, index;
	double val = 0.0;
	for(i = rank; i < N; i+=np) {
		if((i >= k) && (fabs(matrix[i][k]) >= val)) {
			val = fabs(matrix[i][k]);
			index = i;
		}
	}
	if(val != 0.0) {
		return index;
	} else {
		return -1;
	}
}

void copyToBuf(int N, double matrix[N][N], double b[N], double buf[N+1], int k) {
	int i;
	for(i = 0; i < N; i++) {
		buf[i] = matrix[k][i];
	}
	buf[N] = b[k];
}

void copyFromBuf(int N, double matrix[N][N], double b[N], double buf[N+1], int k) {
	int i;
	for(i = 0; i < N; i++) {
		matrix[k][i] = buf[i];
	}
	b[k] = buf[N];
}

void exchangeRow(int N, double matrix[N][N], double b[N], double buf[N+1], int index, int k) {
	int i;
	for(i = 0; i < N; i++) {
		matrix[index][i] = matrix[k][i];
		matrix[k][i] = buf[i];
	}
	b[index] = b[k];
	b[k] = buf[N];
}

void exchangeBufRow(int N, double matrix[N][N], double b[N], double buf[N+1], int index, int k) {
	int i;
	double tmp;
	for(i = 0; i < N; i++) {
		tmp = matrix[index][i];
		matrix[index][i] = buf[i];
		buf[i] = tmp;
	}
	tmp = b[index];
	b[index] = buf[N];
	buf[N] = tmp;
}

void gatherMatrix(int N, double matrix[N][N], double b[N], int rank, int np) {
	MPI_Status status;
	int i;
	for(i = 0; i < N; i++) {
		if(rank == 0 && i % np != 0) {
			MPI_Recv(&matrix[i][0], N, MPI_DOUBLE, i % np, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&b[i], 1, MPI_DOUBLE, i % np, 0, MPI_COMM_WORLD, &status);
		}
		else if(rank == i % np) {
			MPI_Send(&matrix[i][0], N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			MPI_Send(&b[i], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	}
}

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Status status;
	int rank, np, N, i, j, k, index;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	FILE *inFile;
	FILE *outFile;
	if(rank == 0) { //read header
		printf("Running step_2 using %d MPI processes.\n", np);
		utils_open_matrix_files(&inFile, &outFile, argc, argv);
		utils_get_matrix_size(inFile, &N);
	}

	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	double matrix[N][N];
	double buf[N+1]; //Used for temporary storage of pivot row.
	double l[N];
	double b[N];
	double x[N];
	for(i = 0; i < N; i++) {
		b[i] = 1;
		x[i] = 0;
	} 

	if(rank == 0) { //read elements
		utils_read_matrix(inFile, N, matrix);
		int dest;
		for(i = 0; i < N; i++) {
			dest = i % np;
			if(dest == 0)
				continue;
			MPI_Send(&matrix[i][0], N, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
		}
	} else { //If not rank 0, receive your elements through MPI communication.
		for(i = rank; i < N; i+=np) {
			MPI_Recv(&matrix[i][0], N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		}
	}
	struct {double val; int index; long rank;} lp, gp; //local and global pivot in structs.
	for(k = 0; k < N-1; k++) { //Forward elimination.
		index = localPivotIndex(N, matrix, k, rank, np);
		if(index != -1) {
			lp.val = fabs(matrix[index][k]);
			lp.index = index;
		} else {
			lp.val = 0.0;
			lp.index = -1;
		}
		//Determine global pivot from all local pivots.
		MPI_Allreduce(&lp, &gp, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
		gp.rank = gp.index % np;
		if(gp.rank == (k % np)) { //Pivot element is under supervision of owner of row k.
			if(rank == gp.rank) { //This process owns the pivot element.
				copyToBuf(N, matrix, b, buf, gp.index);
				exchangeRow(N, matrix, b, buf, index, k);
			}
		} else { //Pivot row and row k are not owned by the same process.
			if(rank == (k % np)) { //If this process owns the current row.
				copyToBuf(N, matrix, b, buf, k);
				MPI_Send(&buf, N+1, MPI_DOUBLE, gp.rank, 0, MPI_COMM_WORLD);
			} else if (rank == gp.rank){ //If this process owns the pivot row.
				MPI_Recv(&buf, N+1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
				exchangeBufRow(N, matrix, b, buf, gp.index, k);
			}
		}
		MPI_Bcast(&buf[0],N+1, MPI_DOUBLE, gp.rank, MPI_COMM_WORLD);
		if((k % np) == rank && rank != gp.rank) {
			copyFromBuf(N, matrix, b, buf, k);
		}
		i = k+1;
		while (i % np != rank) {
			i++;
		}
		for(; i<N; i+=np) {
			l[i] = matrix[i][k] / buf[k];
			for(j = 0; j < N; j++) {
				matrix[i][j] = matrix[i][j] - l[i] * buf[j];
			}
			b[i] = b[i] - l[i] * buf[N];
		}
	}
	gatherMatrix(N, matrix, b, rank, np);
	double sum;
	if(rank == 0) {
		for(k=N-1; k >= 0; k--) { //Backwards substitution.
			sum = 0.0;
			for(j=k+1; j < N; j++) {
				sum += matrix[k][j] * x[j];
			}
			x[k] = (b[k] - sum) / matrix[k][k];
		}
	}
	if(rank == 0) { //print elements to file
		utils_fprint_array(N, x, outFile);
	}

	MPI_Finalize();
}