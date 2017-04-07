.DEFAULT_GOAL := hello_world
CC = g++
FLAGS = -fopenmp -std=c++11 -Wall -Wextra

hello_world:
	$(CC) $(FLAGS) basic/hello_world.cpp
	@echo "--------------------- DONE -----------------------"

clean:
	rm -f *.out *.o