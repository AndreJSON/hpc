#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "../utils/mm_help.h"

#define max(a,b) \
	({ __typeof__ (a) _a = (a); \
	__typeof__ (b) _b = (b); \
	_a > _b ? _a : _b; })

int p1,p2;

//Returns the column group number of the line k.
int Co(int k, int n) {
	int perGroup = n / p2;
	return k / perGroup;
}

//Returns the row group number of the line k.
int Ro(int k, int n) {
	int perGroup = n / p1;
	return k / perGroup;
}

//Returns the column group number of processor q.
int Cop(int q) {
	return q % p2;
}

//Returns the row group number of processor q.
int Rop(int q) {
	return q / p2;
}

//Returns true if me is part of the specified column group.
int cmember(int me, int group) {
	return Cop(me) == group;
}

//Returns true if me is part of the specified row group.
int rmember(int me, int group) {
	return Rop(me) == group;
}

int startRow(int me, int n) {
	return (me / p2) * (n / p1);
}

int endRow(int me, int n) {
	return ((me+p2) / p2) * (n / p1);
}

int startCol(int me, int n) {
	return Cop(me) * n / p2;
}

int endCol(int me, int n) {
	return (Cop(me) + 1) * n / p2;
}

int crank (int globalRank, int group) {
	return Rop(globalRank);
}

int rrank (int globalRank, int group) {
	return Cop(globalRank); //Column group number will be the rank in the row group.
}

int max_col_loc(int n, int me, double a[n][n], int col) {
	int i, index;
	double val = 0.0;
	for(i = startRow(me, n); i < endRow(me,n); i++) {
		if((i >= col) && (fabs(a[i][col]) >= val)) {
			val = fabs(a[i][col]);
			index = i;
		}
	}
	if(val != 0.0) {
		return index;
	} else {
		return -1;
	}
}

void exchange_row_loc(int n, double a[n][n], double b[n], int r, int k) {
	int i;
	double tmp;
	for(i = 0; i < n; i++) {
		tmp = a[r][i];
		a[r][i] = a[k][i];
		a[k][i] = tmp;
	}
	tmp = b[r];
	b[r] = b[k];
	b[k] = tmp;
}

void copy_row_loc(int n, double a[n][n], double b[n], int k, double buf[n+1]) {
	int i;
	for(i = 0; i < n; i++) {
		buf[i] = a[k][i];
	}
	buf[n] = b[k];
}

void exchange_row_buf(int n, double a[n][n], double b[n], int r, double buf[n+1], int k) {
	int i;
	double tmp;
	for(i = 0; i < n; i++) {
		tmp = a[r][i];
		a[r][i] = buf[i];
		buf[i] = tmp;
	}
	tmp = b[r];
	b[r] = buf[n];
	buf[n] = tmp;
}

int cgrp_leader(int group) {
	return group;
}

int compute_partner(int rowGroup, int rank) {
	return rowGroup * p2 + rank % p2;
}

int compute_size(int n, int k, int rowGroup) {
	if(startCol(rowGroup,n) > k) {
		return n / p2;
	} else if(k < endCol(rowGroup,n)) {
		endCol(rowGroup,n) - (k + 1);
	} else {
		return 0;
	}
}

void compute_elim_fact_loc(int n, double a[n][n], double b[n], int k, double buf[n+1], double elim_buf[n+1], int me) {
	int i;
	for(i=startRow(me,n);i<endRow(me,n);i++) {
		if(i>k) {
			elim_buf[i] = a[i][k] / buf[k];
		}
	}
}

void compute_local_entries(int n,double a[n][n],double b[n],int k,double elim_buf[n+1],double buf[n], int me) {
	int i, j;
	for(i = startRow(me,n); i < endRow(me,n); i++) {
		if(i>k) {
			for(j = startCol(me,n); j < endCol(me,n); j++) {
				a[i][j] = a[i][j] - elim_buf[i] * buf[j];
			}
			b[i] = b[i] - elim_buf[i] * buf[n];
		}
	}
}

void gatherMatrix(int n, double a[n][n], double b[n], int me, MPI_Comm cc) {
	MPI_Status status;
	int i;
	for(i = 0; i < n; i++) { //Up along the column groups
		if(crank(me, Cop(me)) == 0 && Ro(i,n) != Rop(me)) { //If this processor has colRank 0 and row doesnt belong to said processor.
			MPI_Recv(&a[i][startCol(me,n)],n/p2, MPI_DOUBLE, Ro(i,n), 0, cc, &status);
			MPI_Recv(&b[i],1, MPI_DOUBLE, Ro(i,n), 0, cc, &status);
		} else if(Ro(i,n) != Rop(me)) { //If this processor owns the row.
			MPI_Send(&a[i][startCol(me,n)],n/p2, MPI_DOUBLE, 0, 0, cc);
			MPI_Send(&b[i],1, MPI_DOUBLE, 0, 0, cc);
		}
	}
}

void backward_substitution(int n, double a[n][n], double b[n], double x[n]) {
	//TODO
}

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Status status;
	int me, p, n, i, j, k, r, q, q1, size, buf_size, elim_size, psz;
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	utils_get_split(argc, argv, &p, &p1, &p2);
	MPI_Comm cc,rc;
	MPI_Comm_split(MPI_COMM_WORLD,Cop(me),me,&cc);
	MPI_Comm_split(MPI_COMM_WORLD,Rop(me),me,&rc);

	FILE *inFile;
	FILE *outFile;
	if(me == 0) { //read header
		printf("Running step_2 using %d MPI processes.\n", p);
		utils_open_matrix_files(&inFile, &outFile, argc, argv);
		utils_get_matrix_size(inFile, &n);
	}

	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	double a[n][n];
	double buf[n+1]; //Used for temporary storage of pivot row.
	double elim_buf[n+1];
	double b[n];
	double x[n];
	for(i = 0; i < n; i++) {
		b[i] = 1;
		x[i] = 0;
	}
	/***** INITS DONE, NOW READ AND SHARE MATRIX ******/

	if(me == 0) { //read elements
		utils_read_matrix(inFile, n, a);
	}
	int dest;
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j+= n/p2) {
			dest = Co(j,n) + p2*Ro(i,n);
			if(dest == 0)
				continue;
			if(me == 0)
				MPI_Send(&a[i][j], n/p2, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
			else if(dest == me)
				MPI_Recv(&a[i][j], n/p2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		}
	}

	struct {double val; int pvtline;} z, y; //local and global pivot in structs.
	for(k = 0; k < n-1; k++) { //Forward elimination.
		if(cmember(me, Co(k,n))) {
			r = max_col_loc(n,me,a,k);
			z.pvtline = r; z.val = fabs(a[r][k]);
			MPI_Reduce(&z,&y,1,MPI_DOUBLE_INT,MPI_MAXLOC,cgrp_leader(Co(k,n)),cc);
		}
		MPI_Bcast(&y,1,MPI_DOUBLE_INT,cgrp_leader(Co(k,n)),MPI_COMM_WORLD);
		r = y.pvtline;
		if(Ro(k,n) == Ro(r,n)) { //pivot and current row belong to same row group.
			if(rmember(me,Ro(k,n))) {
				if(r != k) exchange_row_loc(n,a,b,r,k);
				copy_row_loc(n,a,b,k,buf);
			}
		}
		else { //pivot and current row are in different row groups.
			if(rmember(me,Ro(k,n))) {
				copy_row_loc(n,a,b,k,buf);
				q = compute_partner(Ro(r,n),me);
				psz = compute_size(n,k,Ro(k,n));
				MPI_Send(&buf[k],psz+1,MPI_DOUBLE,q,0,MPI_COMM_WORLD);
			}
			else if(rmember(me,Ro(r,n))) { //This processor has part of the pivot.
				q = compute_partner(Ro(k,n),me);
				psz = compute_size(n,k,Ro(r,n));
				MPI_Recv(&buf[k],psz+1,MPI_DOUBLE,q,0,MPI_COMM_WORLD,&status);
				exchange_row_buf(n,a,b,r,buf,k);
			}
		}
		for(q=0; q<p; q++) { //Broadcast pivot row.
			if(rmember(q,Ro(r,n)) && cmember(me,Cop(q))) {
				q1 = crank(q,Cop(q)); buf_size = compute_size(n,k,Ro(k,n));
				MPI_Bcast(&buf[k], buf_size, MPI_DOUBLE,q1,cc);
			}
		}
		if((Ro(k,n) != Ro(r,n)) && rmember(me,Ro(k,n))) { //If pivot wasn't already in the buf
			copy_row_loc(n,a,b,k,buf);
		}
		//if this processor is part of the column group containing col k.
		if(cmember(me,Co(k,n))) {
			compute_elim_fact_loc(n,a,b,k,buf,elim_buf,me);
		}
		for(q=0; q<p; q++) { //Broadcast elimination factors.
			if(cmember(q,Co(k,n)) && rmember(me,Rop(q))) {
				q1 = rrank(q,Rop(q)); elim_size = compute_size(n,k,Co(k,n));
				MPI_Bcast(&elim_buf[max(startRow(q,n),k+1)],elim_size,MPI_DOUBLE,q1,rc);
			}
		}
	}
	gatherMatrix(n,a,b,me,cc);
	backward_substitution(n,a,b,x);
	if(me == 0) { //print elements to file
		//får nog fel element nere i vänstra hörnet.
		utils_print_matrix(n, a);
		//utils_print_array(n,x);
	}

	MPI_Finalize();
}