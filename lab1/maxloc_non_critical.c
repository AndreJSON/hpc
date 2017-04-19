#include "stdio.h" // printf
#include "stdlib.h" // rand for instance. 
#include "time.h"   // time(0) to get random seed
#include "omp.h"

#define N 1000000
double x[N];
int main(int argc, char* argv[]){
    // Array with doubles
    int MAX_THREADS = omp_get_max_threads();
    int padding = 128;
    double maxval[MAX_THREADS*padding]; int maxloc[MAX_THREADS*padding];
    int mloc;
    double mval;
    int i;
    for(i = 0; i < MAX_THREADS*padding; i++) {
    	maxval [i] = 0.0;
    	maxloc [i] = 0;
    }
    // fill the array
    srand(time(0)); // seed
    for(i=0; i < N;i++){
	       // Generate random number between 0 and 1
	       x[i] = ((double)(rand()) / RAND_MAX)*((double)(rand()) / RAND_MAX)*((double)(rand()) / RAND_MAX)*1000;
	}

    // calculate first serially to get values to checking the results obtianed with OpenMP
    for (i=0; i < N; i++){
		if (x[i] > mval) {
				mval = x[i]; 
				mloc = i; 
		}
	}
	double maxval_1 = mval; int maxloc_1 = mloc;

	// OpenMP Version
	double start_time, run_time;
	start_time = omp_get_wtime();
	
	#pragma omp parallel for
	for (i=0; i < N; i++){
		int threadid = omp_get_thread_num();
		if (x[i] > maxval[threadid*padding]) {
			maxval[threadid*padding] = x[i]; 
			maxloc[threadid*padding] = i;
		}
	}
	int maxloc2 = 0;
	double maxval2 = 0.0;
	for(i=0;i<MAX_THREADS*padding;i++){
		if(maxval[i] > maxval2){
			maxval2 = maxval[i];
			maxloc2 = maxloc[i];
		}
	}
	run_time = omp_get_wtime() - start_time;
    printf("maxloc computation in %f seconds\n",run_time);
    printf("maxval (omp)   = %f maxloc (omp)   = %d \n",maxval2, maxloc2);
    printf("maxval (s)     = %f maxloc (s)     = %d \n",maxval_1, maxloc_1);
    if (maxloc_1 != maxloc2)
    	printf("Test failed\n");

	return 0;
}