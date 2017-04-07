#include "hello_world.hpp"
using std::cout; using std::endl;

int main(void) {
	#pragma omp parallel num_threads(4)
	{
		cout << "Hello"<<  endl;
		cout << "World!" << endl;
		cout << "From" << endl;
		cout << omp_get_thread_num() << endl;
	}
	return 0;
}