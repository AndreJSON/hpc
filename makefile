CC = g++
FLAGS = -fopenmp -Wall -Wextra

hello_world:
	$(CC) $(FLAGS) basic/hello_world.cpp
	@echo "--------------------- DONE -----------------------"

stream_benchmark:
	$(CC) $(FLAGS) basic/stream_benchmark.cpp
	@echo "--------------------- DONE -----------------------"

clean:
	rm -f *.out *.o