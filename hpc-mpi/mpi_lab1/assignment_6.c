//NONBLOCKING

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define DEBUG_CODE 0
#define DEBUG_RESULT 0

void gather_ring_blocking( float *x, int blocksize, float *y)
{
    int i, p, myrank, succ, pred;
    int send_offs, recv_offs;
    MPI_Status status;
    MPI_Request send_req, recv_req;

    //Get total number of processes.
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    //Get rank (id) of this process.
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


    //Move blocks from vector x to correct position in local version of collection vector.
    for (i = 0; i < blocksize; i++)
        y[i + myrank*blocksize] = x[i];

    succ = (myrank+1) % p; //Get successor process.
    pred = (myrank-1+p) % p; //Get predecessor process.

    for (i = 0; i < p-1; i++)
    {
        send_offs = ((myrank-i+p)%p)*blocksize; //Which index of y we will start sending from
        recv_offs = ((myrank-i-1+p)%p)*blocksize; //Which index of y we will start receiving from.

        MPI_Isend(y + send_offs, blocksize, MPI_FLOAT, succ, 0, MPI_COMM_WORLD, &send_req);
        MPI_Irecv(y + recv_offs, blocksize, MPI_FLOAT, pred, 0, MPI_COMM_WORLD, &recv_req);

        MPI_Wait(&send_req, &status);
        MPI_Wait(&recv_req, &status);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char *argv[]) {
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
    gather_ring_blocking(x, blocksize, y);

    if(rank == 0)
        printf("Time %.10f\n", (MPI_Wtime()-startTime));

    if(rank == 1 && DEBUG_RESULT) {
        for(i = 0; i < blocksize*np; i++) {
            printf("%f\n", y[i]);
        }
    }
    MPI_Finalize();
}