//ring mesh hyper

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define DEBUG_RESULT 0

#define BUFF_SIZE 1000000

void broadcast_ring(void *buffer, int count, int root, MPI_Comm comm) {
	int i, rank, np;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	if(rank != 0) {
		MPI_Recv(buffer, count, MPI_INT, rank-1, 0, comm, &status);
	}
	if(rank != np-1) {
		MPI_Send(buffer, count, MPI_INT, rank+1, 0, comm);
	}
	MPI_Barrier(comm);
}

void broadcast_mesh(void *buffer, int count, int root, MPI_Comm comm) {
	int i, rank, np, xSize;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	xSize = (int)sqrt(np);

	if(xSize > 1) { //X direction communication should happen.
		if(rank > 0 && rank < xSize) {
			MPI_Recv(buffer, count, MPI_INT, rank-1, 0, comm, &status);
		}
		if(rank < xSize-1) {
			MPI_Send(buffer, count, MPI_INT, rank+1, 0, comm);
		}
	}
	//Y direction
	if(rank >= xSize) {
		MPI_Recv(buffer, count, MPI_INT, rank-xSize, 0, comm, &status);
	}
	if(rank < np-xSize) {
		MPI_Send(buffer, count, MPI_INT, rank+xSize, 0, comm);
	}
	MPI_Barrier(comm);
}

void broadcast_hypercube(void *buffer, int count, int root, MPI_Comm comm) {
	int i, rank, np, xSize;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	for(i = 1; i < np; i*=2) {
		if(rank < i) {
			MPI_Send(buffer, count, MPI_INT, rank+i, 0, comm);
		} else if (rank < i*2) {
			MPI_Recv(buffer, count, MPI_INT, rank-i, 0, comm, &status);
		}
	}
	MPI_Barrier(comm);
}

int main(int argc, char *argv[]) {
	if ( argc < 2 )
		exit(-1);
	int rank, msgLen, numTests, i;
	double startTime;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	msgLen = atoi(argv[1]);
	numTests = atoi(argv[2]);

	int message[msgLen];

	if (!rank) {
		for (i = 0; i < msgLen; i++) {
			message[i] = i;
		}
	}

	if(!rank) {
		startTime = MPI_Wtime();
	}
	for (i = 0; i < numTests; i++) {
		MPI_Bcast(message, msgLen, MPI_INT, 0, MPI_COMM_WORLD);
	}

	if(!rank) {
		printf("Regular   %.10f seconds\n", (MPI_Wtime()-startTime)/numTests);
		startTime = MPI_Wtime();
	}
	for (i = 0; i < numTests; i++) {
		broadcast_ring(message, msgLen, 0, MPI_COMM_WORLD);
	}

	if(!rank) {
		printf("Ring      %.10f seconds\n", (MPI_Wtime()-startTime)/numTests);
		startTime = MPI_Wtime();
	}
	for (i = 0; i < numTests; i++) {
		broadcast_mesh(message, msgLen, 0, MPI_COMM_WORLD);
	}

	if(!rank) {
		printf("Mesh      %.10f seconds\n", (MPI_Wtime()-startTime)/numTests);
		startTime = MPI_Wtime();
	}
	for (i = 0; i < numTests; i++) {
		broadcast_hypercube(message, msgLen, 0, MPI_COMM_WORLD);
	}

	if(!rank) {
		printf("Hypercube %.10f seconds\n", (MPI_Wtime()-startTime)/numTests);
	}

	if(DEBUG_RESULT) {
		for(i = 0; i < msgLen; i++) {
			if(message[i] != i)
				printf("FAULTY BROADCAST, rank %d expected %d but has value %d\n", rank, i, message[i]);
		}
	}

	MPI_Finalize();
}