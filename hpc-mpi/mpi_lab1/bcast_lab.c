
#include <stdlib.h>
#include <mpi.h>

void
broadcast_ring(void     *buffer,
               int       count,
               int       root,
               MPI_Comm  comm)
{
    MPI_Datatype datatype = MPI_INT; // Use this MPI function calls.

    // Implementation goes here.
}

void
broadcast_mesh(void     *buffer,
               int       count,
               int       root,
               MPI_Comm  comm)
{
    MPI_Datatype datatype = MPI_INT; // Use this MPI function calls.

    // Implementation goes here.
}

void
broadcast_hypercube(int      *buffer,
                    int       count,
                    int       root,
                    MPI_Comm  comm)
{
    MPI_Datatype datatype = MPI_INT; // Use this MPI function calls.

    // Implementation goes here.
}

int main(int argc, char *argv[])
{
    if ( argc < 2 )
        exit(-1);

    MPI_Init(&argc, &argv);

    int i, cur_msg_len, p, myrank;

    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int max_msg_len = atoi(argv[1]);
    int nr_tests    = atoi(argv[2]);

    int message[max_msg_len];

    if ( myrank == 0 )
       for (i = 0; i < max_msg_len; i++) message[i] = i;

    for (cur_msg_len = 1; cur_msg_len <= max_msg_len; cur_msg_len *= 2)
    {
        // MPI_Wtime() -> take start time on rank 0

        for (i = 0; i < nr_tests; i++)
        {
            MPI_Bcast(message, cur_msg_len, MPI_INT, 0, MPI_COMM_WORLD);
            // broadcast_ring(message, cur_msg_len, 0, MPI_COMM_WORLD);
            // broadcast_mesh(message, cur_msg_len, 0, MPI_COMM_WORLD);
            // broadcast_hypercube(message, cur_msg_len, 0, MPI_COMM_WORLD);
        }

        // MPI_Wtime() -> take end time on rank 0
        // Calculate on rank 0:
        // time = (end_time - start_time)/nr_tests;
        // print message length and time
    }

    MPI_Finalize();
    return 0;
}