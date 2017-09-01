//MATRIX MUL

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define DEBUG_RESULT 0
#define DEBUG_CODE 0
#define DEBUG_MATRIX 0

#define BLOCK_DIST 1

#define BUFFER_SIZE 1000000

#define MIN(a,b) (((a)<(b))?(a):(b))

void init(size_t N, size_t M, int matrix[N][M], int vector[M]) {
	int i, j;
	for(i = 0; i < N; i++) {
		for(j = 0; j < M; j++) {
			matrix[i][j] = i + j;
			if(i == 0) {
				vector[j] = j;
			}
		}
	}
}

void scatterBlock(size_t N, size_t M, int matrix[N][M], int vector[M], int np, int rank, int distribution[np+1]) {
	MPI_Status status;
	int i;
	
	for(i = 0; i < np+1; i++) {
		if(i == np) {
			distribution[i] = N;
		} else {
			distribution[i] = (N / np) * i + MIN(N%np,i);
		}
		
		if(DEBUG_CODE && rank == 0) {
			printf("rank %d starts at row %d\n", i, distribution[i]);
		}
	}

	if(rank == 0) {
		for(i = 1; i < np; i++) {
			MPI_Bsend(&matrix[distribution[i]][0], M*(distribution[i+1] - distribution[i]), MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Bsend(&vector[0], M, MPI_INT, i, 0, MPI_COMM_WORLD);
		}
	} else {
		MPI_Recv(&matrix[distribution[rank]][0], M*(distribution[rank+1] - distribution[rank]), MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&vector[0], M, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	}
}

void scatterCyclic(size_t N, size_t M, int matrix[N][M], int vector[M], int np, int rank) {
	MPI_Status status;
	int i;

	if(rank == 0) {
		for(i = 0; i < N; i++) {
			if(i%np == 0) //0 dont send to 0
				continue;
			if(i < np) //Only do the first lap
				MPI_Bsend(&vector[0], M, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Bsend(&matrix[i][0], M, MPI_INT, i%np, 0, MPI_COMM_WORLD);
		}
	} else {
		MPI_Recv(&vector[0], M, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		for(i = rank; i < N; i+=np) {
			MPI_Recv(&matrix[i][0], M, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		}
	}
}

void multiplyBlock(size_t N, size_t M, int matrix[N][M], int vector[M], int np, int rank, int distribution[np+1], int result[N]) {
	int i, j;
	for(i = distribution[rank]; i < distribution[rank+1]; i++) {
		result[i] = 0;
		for(j = 0; j < M; j++) {
			result[i] += matrix[i][j] * vector[j];
		}
	}
}

void multiplyCyclic(size_t N, size_t M, int matrix[N][M], int vector[M], int np, int rank, int distribution[np+1], int result[N]) {
	int i, j;
	for(i = rank; i < N; i+=np) {
		result[i] = 0;
		for(j = 0; j < M; j++) {
			result[i] += matrix[i][j] * vector[j];
		}
	}
}

void shareResultBlock(size_t N, int np, int distribution[np+1], int result[N]) {
	int i;
	for(i = 0; i < np; i++) {
		MPI_Bcast(&result[distribution[i]],distribution[i+1]-distribution[i],MPI_INT,i,MPI_COMM_WORLD);
	}
}

void shareResultCyclic(size_t N, int np, int distribution[np+1], int result[N]) {
	int i;
	for(i = 0; i < N; i++) {
		MPI_Bcast(&result[i],1,MPI_INT,i%np,MPI_COMM_WORLD);
	}
}

int main(int argc, char *argv[]) {
	MPI_Init(NULL, NULL);
	MPI_Buffer_attach(malloc(BUFFER_SIZE), BUFFER_SIZE);
	int np, rank, M, N, i, j;
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	N = strtol(argv[1],NULL,10); //Number of rows
	M = strtol(argv[2],NULL,10); //Length of rows

	int matrix[N][M], vector[M], distribution[np+1], result[N];

	//Initialize matrix
	if(rank == 0) {
		init(N, M, matrix, vector);
		//Print initial matrix
		if(DEBUG_MATRIX) {
			for(i = 0; i < N; i++) {
				for(j = 0; j < M; j++) {
					printf("%d",matrix[i][j]);
				}
				printf("\n");
			}
		}
	}

	if(BLOCK_DIST) {
		//Scatter matrix rows
		scatterBlock(N, M, matrix, vector, np, rank, distribution);
		//Calculate your part of the solution
		multiplyBlock(N, M, matrix, vector, np, rank, distribution, result);
		//Share your part of the solution with all other processes and get their parts.
		shareResultBlock(N, np, distribution, result);
	} else {
		//Scatter matrix rows
		scatterCyclic(N, M, matrix, vector, np, rank);
		//Calculate your part of the solution
		multiplyCyclic(N, M, matrix, vector, np, rank, distribution, result);
		//Share your part of the solution with all other processes and get their parts.
		shareResultCyclic(N, np, distribution, result);
	}

	if(DEBUG_CODE && rank == 1) {
		for(i = 0; i < M; i++) {
			printf("%d \n", vector[i]);
		}
	}

	if(DEBUG_RESULT && rank == 1) {
		for(i = 0; i < N; i++) {
			printf("%d \n", result[i]);
		}
	}

	MPI_Finalize();
}