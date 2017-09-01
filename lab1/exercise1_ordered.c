#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main (void) {
	int nthreads, tid, i;
	nthreads = omp_get_max_threads();
	/* Fork a team of threads giving them their own copies of variables */

	#pragma omp parallel for ordered private(i, tid)
        for (i=0; i < nthreads;i++){
	/* Obtain thread number */
		tid = omp_get_thread_num();
		#pragma omp ordered
		printf("Hello from thread = %d\n", tid);

		/* Only master thread does this */
		if (tid == nthreads-1)
		{
		nthreads = omp_get_num_threads();
		printf("Number of threads = %d\n", nthreads);
		}

	}
			/* All threads join master thread and disband */

}

