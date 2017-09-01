.DEFAULT_GOAL = all
MPI_CC = mpicc
MPI_EXEC = mpiexec
MPI_CORE_SETTINGS = -n 4
OMP_CC = gcc
OMP_FLAGS = -fopenmp
UTILS = utils/mmio.c utils/mm_help.c -lm


all: step2

step2:
	$(MPI_CC) -o step_2/step_2.x step_2/step_2.c $(UTILS)

step3:
	$(MPI_CC) -o step_3/step_3.x step_3/step_3.c $(UTILS)

run: test3

test2: step2
	$(MPI_EXEC) $(MPI_CORE_SETTINGS) step_2/step_2.x -g 2x1 -i data_in/matrix_0512.mm -o test_result/result_0512.mm

test3: step3
	$(MPI_EXEC) $(MPI_CORE_SETTINGS) step_3/step_3.x -g 2x2 -i data_in/testin_b.mm -o test_result/result_b.mm