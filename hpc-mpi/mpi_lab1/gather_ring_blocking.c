
#include <mpi.h>

#define NP 4            // number of ranks
#define BLOCK_SIZE  1

void gather_ring_blocking( float *x, int blocksize, float *y)
{
    int i, p, myrank, succ, pred;
    int send_offs, recv_offs;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    for (i = 0; i < blocksize; i++)
        y[i + myrank*blocksize] = x[i];

    succ = (myrank+1) % p;
    pred = (myrank-1+p) % p;

    for (i = 0; i < p-1; i++)
    {
        send_offs = ((myrank-i+p)%p)*blocksize;
        recv_offs = ((myrank-i-1+p)%p)*blocksize;

        MPI_Send(y + send_offs, blocksize, MPI_FLOAT, succ, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(y + recv_offs, blocksize, MPI_FLOAT, pred, 0,
                 MPI_COMM_WORLD, &status);
    }
}

int main(int argc, char *argv[])
{
    float  x[BLOCK_SIZE];
    float  y[NP*BLOCK_SIZE];

    gather_ring_blocking(x, BLOCK_SIZE, y);

    return 0;
}
