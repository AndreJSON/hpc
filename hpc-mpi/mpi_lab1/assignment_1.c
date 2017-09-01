#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define DEBUG_CODE 0
#define DEBUG_RESULT 0

#define BUFF 100000000

void gather_ring_blocking( float *x, int blocksize, float *y)
{
    int i, p, myrank, succ, pred;
    int send_offs, recv_offs;
    MPI_Status status;
    MPI_Buffer_attach(malloc(BUFF), BUFF);


    //Get total number of processes.
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    //Get rank (id) of this process.
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // Keep track of which process this process has sent and received from
    //int recieved[p-1];
    //int sent[p-1];


    //Move blocks from vector x to correct position in local version of collection vector.
    for (i = 0; i < blocksize; i++)
        y[i + myrank*blocksize] = x[i];

    succ = (myrank+1) % p; //Get successor process.
    pred = (myrank-1+p) % p; //Get predecessor process.

    for (i = 0; i < p-1; i++)
    {
        send_offs = ((myrank-i+p)%p)*blocksize; //Which index of y we will start sending from
        recv_offs = ((myrank-i-1+p)%p)*blocksize; //Which index of y we will start receiving from.

        MPI_Bsend(y + send_offs, blocksize, MPI_FLOAT, succ, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(y + recv_offs, blocksize, MPI_FLOAT, pred, 0,
                 MPI_COMM_WORLD, &status);
        //recieved[i] = ((myrank-i-1+p)%p);
        //sent[i] = ((myrank-i+p)%p);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    //Print received to stdout and sent to stderr, run with > received.txt 2> sent.txt
    /*for(i =0; i < p-1;i++){
    	printf("Rank %d recieved from: %d\n", myrank, recieved[i]);
    	fprintf(stderr,"Rank %d sent data belonging to: %d\n", myrank, sent[i]);
    }*/
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
    gather_ring_blocking(x, blocksize, y);

    if(rank == 0) {
        printf("Time %.10f\n", (MPI_Wtime()-startTime));
    }

    if(rank == 1 && DEBUG_RESULT) {
        for(i = 0; i < blocksize*np; i++) {
            printf("%f\n", y[i]);
        }
    }
    MPI_Finalize();
}