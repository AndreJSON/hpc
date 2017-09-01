#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

#define DEBUG_CODE 0
#define DEBUG_RESULT 0

void gather_ring_blocking( float *x, int * blocksizes, float *y, int total_blocksize)
{
    int i, p, myrank, succ, pred;
    int send_offs, recv_offs;
    MPI_Status status;


    //Get total number of processes.
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    //Get rank (id) of this process.
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int blockstart = 0;
    //Calculate where this ranks data starts in y.
    for(i = 0; i<myrank; i++){
    	blockstart += blocksizes[i];
    }
    //Move blocks from vector x to correct position in local version of collection vector.
    for (i = 0; i < blocksizes[myrank]; i++)
        y[i + blockstart] = x[i];

    succ = (myrank+1) % p; //Get successor process.
    pred = (myrank-1+p) % p; //Get predecessor process.
    recv_offs = blockstart; 
    for (i = 0; i < p-1; i++)
    {
        send_offs = recv_offs; //Which index of y we will start sending from
        recv_offs = (recv_offs - blocksizes[(myrank - i - 1 + p) % p] + total_blocksize) % total_blocksize; //Which index of y we will start receiving from.
        if(DEBUG_CODE == 1)
            printf("rank %i har so %i och ro %i\n", myrank, send_offs, recv_offs);

        MPI_Send(y + send_offs, blocksizes[(myrank - i + p) % p], MPI_FLOAT, succ, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(y + recv_offs, blocksizes[(myrank - i - 1 + p) % p], MPI_FLOAT, pred, 0,
                 MPI_COMM_WORLD, &status);
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(NULL, NULL);
    int i, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np); //Get number of processes.
    srand(time(NULL) + rank);
    int max_blocksize = 1;
    if(argc > 1) {
            max_blocksize = strtol(argv[1],NULL,10);
        }
    int blocksize = (rand() % max_blocksize) + 1; //Choose random blocksize between 1 and 4.

    float  x[blocksize];

//    printf("Rank %d has blocksize: %d\n", rank,blocksize);
    for(i = 0; i < blocksize; i++) {
        x[i] = rank + 0.01*i;
    }

    int blocksizes[np];
    blocksizes[rank] = blocksize;
    for(i =0; i < np; i++){
    	MPI_Bcast(&blocksizes[i],1,MPI_INT,i,MPI_COMM_WORLD); //Let every process broadcast their blocksize one after another.
    }


	//Communicate blocksize to all processes

    int total_blocksize = 0;
    for(i = 0; i < np; i++){
    	total_blocksize += blocksizes[i];
    }

    float  y[total_blocksize];
    gather_ring_blocking(x, blocksizes, y, total_blocksize);

    if(rank == 0 && DEBUG_RESULT) {
        for(i = 0; i < total_blocksize; i++) {
            printf("%f\n", y[i]);
        }
        printf("Blocksizes: ");
        for(i = 0; i<np;i++){
         	printf("%i,", blocksizes[i]);
         }
        printf("\n");
    }
    MPI_Finalize();
}
