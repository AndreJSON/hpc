//TORUS

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define DEBUG_CODE 0
#define DEBUG_CODE2 0
#define DEBUG_X 0
#define DEBUG_Y 0
#define DEBUG_RESULT 0

#define BUFF 100000000

void gather_ring_blocking( float *x, int blocksize, float *y, int p, int myrank)
{
    int i, succ, pred;
    int send_offs, recv_offs;
    int n = sqrtl((double)p);
    MPI_Status status;
    MPI_Buffer_attach(malloc(BUFF), BUFF);


    //Move blocks from vector x to correct position in local version of collection vector.
    for (i = 0; i < blocksize; i++)
        y[i + myrank*blocksize] = x[i];

    //////// X DIRECTION ///////////
    if(myrank % n == (n - 1))
        succ = myrank - n + 1;
    else
        succ = myrank + 1;
    if(myrank % n == 0)
        pred = myrank + n - 1;
    else
        pred = myrank - 1;

    if(DEBUG_X == 1 && DEBUG_CODE == 1) {
        printf("x-wise rank %d has succ %d and pred %d\n", myrank, succ, pred);
    }

    for (i = 0; i < n-1; i++) {
        send_offs = (myrank - (myrank%n) + ((myrank-i+n)%n))*blocksize;
        recv_offs = (myrank - (myrank%n) + ((myrank-1-i+n)%n))*blocksize;
        if(DEBUG_X == 1 && DEBUG_CODE2 == 1) {
            printf("x-wise,round %d rank %d has send_offs %d and recv_offs %d\n", i, myrank, send_offs, recv_offs);
        }

        MPI_Bsend(y + send_offs, blocksize, MPI_FLOAT, succ, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(y + recv_offs, blocksize, MPI_FLOAT, pred, 0,
                 MPI_COMM_WORLD, &status);
    }

    //////// Y DIRECTION ///////////
    if(myrank >= p-n)
        succ = myrank % n;
    else
        succ = myrank + n;
    if(myrank < n)
        pred = p - n + myrank;
    else
        pred = ((myrank - n)+p) % p;

    if(DEBUG_Y == 1 && DEBUG_CODE == 1) {
        printf("y-wise rank %d has succ %d and pred %d\n", myrank, succ, pred);
    }

    for (i = 0; i < n-1; i++) {
        send_offs = ((p + myrank - (myrank%n) - (i*n))%p)*blocksize;
        recv_offs = ((p + myrank - (myrank%n) - (i*n) - n)%p)*blocksize;
        if(DEBUG_Y == 1 && DEBUG_CODE2 == 1) {
            printf("y-wise,round %d rank %d has send_offs %d and recv_offs %d\n", i, myrank, send_offs, recv_offs);
        }

        MPI_Bsend(y + send_offs, blocksize*n, MPI_FLOAT, succ, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(y + recv_offs, blocksize*n, MPI_FLOAT, pred, 0,
                 MPI_COMM_WORLD, &status);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char *argv[])
{
    MPI_Init(NULL, NULL);

    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np); //Get number of processes.
    int blocksize = 1;
    if(argc > 1) {
        blocksize = strtol(argv[1],NULL,10);
    }
    float  x[blocksize];
    int i, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for(i = 0; i < blocksize; i++) {
        x[i] = rank + 0.01*i;
    }
    float  y[np*blocksize];

    int startTime = MPI_Wtime();
    gather_ring_blocking(x, blocksize, y, np, rank);

    if(rank == 0) {
        printf("Time %.10f\n", (MPI_Wtime()-startTime));
    }
    if(rank == 0 && DEBUG_RESULT) {
        for(i = 0; i < blocksize*np; i++) {
            printf("%f\n", y[i]);
        }
    }
    MPI_Finalize();
}