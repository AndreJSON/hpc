#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

#define DEBUG_CODE 0
#define DEBUG_RESULT 0


void gather_ring_blocking( float *x, int blocksize, float **y, int max_blocksize, int* y_count)
{

    int i, p, myrank, succ, pred;
    MPI_Status status;

	
    //Get total number of processes.
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    //Get rank (id) of this process.
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    float* tmp = (float *) malloc(4*blocksize);
    //Move blocks from vector x to correct position in local version of collection vector.
    for (i = 0; i < blocksize; i++){
        tmp[i] = x[i];
	}
    y[myrank] = tmp;
	y_count[myrank] = blocksize;

    succ = (myrank+1) % p; //Get successor process.
    pred = (myrank-1+p) % p; //Get predecessor process.
    for (i = 0; i < p-1; i++)
    {
    	float recv_temp[max_blocksize];

    	//printf("%lu\n" ,sizeof(y[(myrank - i + p) % p])/sizeof(y[myrank][0]));
        MPI_Send(y[(myrank - i + p) % p], y_count[(myrank-i+p)%p], MPI_FLOAT, succ, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(recv_temp, max_blocksize, MPI_FLOAT, pred, 0,
                 MPI_COMM_WORLD, &status);
		int recieved_rank = (myrank-i-1+p)%p;
		MPI_Get_count(&status, MPI_INT, &y_count[recieved_rank]); //status lets us see how much data was actually sent.
		
        tmp = (float *) malloc(4*y_count[recieved_rank]);
        y[recieved_rank] = tmp;
		for(int j = 0; j< y_count[recieved_rank]; j++){
			y[recieved_rank][j] = recv_temp[j];
		}
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(NULL, NULL);

    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np); //Get number of processes.
	int i, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int max_blocksize = 1;
    if(argc > 1) {
            max_blocksize = strtol(argv[1],NULL,10);
        }
	srand(time(NULL) + rank);
    int blocksize = (rand() % max_blocksize) + 1;

    float  x[blocksize];
    for(i = 0; i < blocksize; i++) {
        x[i] = rank + 0.01*i;
    }
    float * y[np];
	int y_count[np];
    gather_ring_blocking(x, blocksize, y, max_blocksize,y_count);

	
    for(i = 0; i < np; i++) {
        if(rank == 0 && DEBUG_RESULT) {
    		for(int j = 0; j < y_count[i];j++){
                printf("%f\n", y[i][j]);
      	    }
        }
        free(y[i]);
	}

    MPI_Finalize();
}
