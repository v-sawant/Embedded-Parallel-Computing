/*
 * HOST side program
 * QR factorization using Householder Transformation on 4x4 Matrix
 * Developed by Vinay Sawant
 */
 
/* Epiphany Host Application */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <math.h>

#include <e-loader.h>  // Comment this line if you are using old SDK ver 5.13.09.10

#include "e-hal.h"

#define FAIL(...) { fprintf(stderr, __VA_ARGS__); exit(1); }
#define SHM_OFFSET 0x01000000

#define IN_ROWS 4//row size of input matrix
#define IN_COLS 4//column size of input matrix

#define COREUSED 1//16//8//4//2//1//Maximum speedup at (unsigned int) IN_ROWS/2, thereafter no significant improvement in speedup.
#define NB_STAGES 11 //number of stages used (tentative)
#define Double float //double

typedef struct{
	Double flag;
	Double stage;
	Double test;
	Double H[IN_ROWS][IN_ROWS][IN_COLS];
	Double Q[IN_ROWS][IN_ROWS];
	Double R[IN_ROWS][IN_COLS];
	Double A[IN_ROWS][IN_COLS];
	unsigned long long total_cycles[COREUSED];
	unsigned long long stage_cycles[NB_STAGES];
}shm_t;

Double in[IN_ROWS][IN_COLS] = {
	{5.124575, 7.898948, 7.946213, 3.592056}, {2.522771,1.012173,8.776333, 6.716550}, {8.535061, 1.283802, 0.308171, 1.881512}, {3.747190, 1.811495, 6.834885, 9.154308}
};

Double I[IN_ROWS][IN_ROWS] = {
	{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}
};

void matrix_clear(Double* matrix, int m, int n)
{
	//zero all elements of a
	for(int i = 0; i < m; i++){
	    for(int j = 0; j < n; j++){
	        *((matrix+i*n)+j) = 0;
	    }
	}
    return;
}

void matrix_transfer2(Double* x, Double* y, unsigned int m, unsigned int n)
{
	//double y[x->m][x->n]; 
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			*((y+i*n)+j) = *((x+i*n)+j);
			//y[i][j] = x->v[i][j];
	return;
}

void matrix_show2(Double* y, unsigned int m, unsigned int n)
{
	for(int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			printf(" %8.5f", *((y+i*n)+j));
		}
		printf("\n");
	}
	printf("\n");
}

void vector_show(Double *v, int vsize)
{
	for (int j = 0; j < vsize; j++) {
		printf(" %8.3f", v[j]);
	}
	printf("\n");
}

// ||x|| 
Double vnorm(Double x[], int n)
{
	Double sum = 0;
	for (int i = 0; i < n; i++) sum += x[i] * x[i];
	return sqrt(sum);
}

// c = a + b * s 
void vmadd(Double a[], Double b[], Double s, Double c[], int n)
{
	for (int i = 0; i < n; i++)
		c[i] = a[i] + s * b[i];
	return;
}

// y = x / d 
void vdiv(Double x[], Double d, Double y[], int n)
{
	for (int i = 0; i < n; i++) 
	    y[i] = x[i] / d;
	return;
}

void _minor(Double x[][IN_COLS], Double y[][IN_COLS], int d)
{
	//Double m[IN_ROWS][IN_COLS];
	matrix_clear((Double*)y, IN_ROWS, IN_COLS);
	for(int i=0; i<d; i++)
		y[i][i] = 1;
	for(int i=d; i< IN_ROWS; i++)
	for(int j=d; j< IN_COLS; j++)
		y[i][j] = x[i][j];
	return;
}

/* take c-th column of m, put in v */
void getColumn(Double *x, Double *v, unsigned int col_nb)
{
	for (unsigned int i = 0; i < IN_ROWS; i++)
		v[i] = *((x+i*IN_COLS)+col_nb);
	return;
}

void getHouseholderMatrix(Double v[], Double h[][IN_ROWS])
{
    matrix_clear((Double*)h, IN_ROWS, IN_ROWS);
	for (int i = 0; i < IN_ROWS; i++)
		for (int j = 0; j < IN_ROWS; j++)
			h[i][j] = -2 * v[i] * v[j];
	for (int i = 0; i < IN_ROWS; i++)
		h[i][i] += 1;
 
	return;
}

void matrix_mul(Double* x, Double* y, Double* z, unsigned int xm, unsigned int xn, unsigned int ym, unsigned int yn)
{
	if (xn != ym) return;
	
    matrix_clear((Double*)z, xm, yn);	
	for (int i = 0; i < xm; i++)
		for (int j = 0; j < yn; j++)
			for (int k = 0; k < xn; k++)
				*((z+i*yn)+j) += *((x+i*xn)+k) * *((y+k*yn)+j);
	return;
}

void matrix_transpose(Double* matrix, unsigned int m, unsigned int n)
{
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < i; j++) {
			Double t = *((matrix+i*n)+j);
			*((matrix+i*n)+j) = *((matrix+j*n)+i);
			*((matrix+j*n)+i) = t;
		}
	}
}

//Main - host program starts from this function
int main(int argc, char *argv[])
{
fprintf(stderr, "Starting Epiphany Host Code....\n");

fprintf(stderr, "BP1\n");
	char *filename = "bin/testcode.srec";
	e_epiphany_t dev;
	e_mem_t      emem;
fprintf(stderr, "BP2\n");
	e_set_host_verbosity(H_D0);
	e_set_loader_verbosity(L_D0);
fprintf(stderr, "BP3\n");
	/* init board, reset epiphany, open workgroup */
	if(e_init(NULL) != E_OK)
		FAIL("Can't init!\n");
fprintf(stderr, "BP4\n");
	e_reset_system();
	if(e_open(&dev, 0, 0, 4, 4) != E_OK)
		FAIL("Can't open!\n");
fprintf(stderr, "BP5\n");
	if(e_alloc(&emem, SHM_OFFSET, sizeof(shm_t)) != E_OK)
		FAIL("Can't alloc!\n");

	shm_t shm;
printf("Size of shm:%llu\n", (unsigned long long) (sizeof shm));

	shm.flag = 0;
	shm.test = 11;
	shm.stage = 0;
	
	//Initiate stage cycles to zero
	for(unsigned int t=0; t< NB_STAGES; t++){
	shm.stage_cycles[t] = 0;
	}
/*//No need to  initiate any matrices in the host program	
	matrix_transfer2((Double*)I, (Double*)shm.Q, IN_ROWS, IN_ROWS);
	matrix_show2((Double*)shm.Q, IN_ROWS, IN_ROWS);
	matrix_transfer2((Double*)I, (Double*)shm.R, IN_ROWS, IN_COLS);
	matrix_show2((Double*)shm.R, IN_ROWS, IN_COLS);
	matrix_transfer2((Double*)in, (Double*)shm.A, IN_ROWS, IN_COLS);
	matrix_show2((Double*)shm.A, IN_ROWS, IN_COLS);
*///
	if(e_write(&emem, 0, 0, (off_t)0, &shm, sizeof(shm_t)) == E_ERR)
			FAIL("Can't clear shm!\n");
fprintf(stderr, "BP-for_write\n");
	
	// Load all required cores at once and start
		for(int pindx = 0; pindx<COREUSED; pindx++){
			int i = pindx/4;
			int j = pindx%4;
		/* load program */
		if(e_load(filename, &dev, i, j, E_TRUE) != E_OK)
			FAIL("Can't load!\n");
fprintf(stderr, "BP-for_load:core(%d,%d)\n", i, j);
		}
fprintf(stderr, "BP-stage:== %f\n",shm.stage);
fprintf(stderr, "BP-flag:== %f\n",shm.flag);
	while(1)
	{
		/* poll shm states for changes */
//		fprintf(stderr, "Polling shared memory.\n");

		/* read shm */
		if(e_read(&emem, 0, 0, (off_t)0, &shm, sizeof(shm_t)) == E_ERR)
				FAIL("Can't poll!\n");
//fprintf(stderr, "BP-for_read\n");
		//shm.flag = COREUSED;//TEST

fprintf(stderr, "BP-for_read, shm.test: %f\n", shm.test);// testing
fprintf(stderr, "BP-stage: %f\n",shm.stage);// testing
fprintf(stderr, "BP-flag: %f\n",shm.flag);// testing
		//Check for condition
		if (shm.flag == COREUSED)
		{
			for(int pnb = 0; pnb < COREUSED; pnb++)
			{
				
				printf("\nCore number:%d; Total cycles:%llu\n", pnb, shm.total_cycles[pnb]);
			}
			break;//break if all core programs are done processing
		}
	}
fprintf(stderr, "BP-WHILE1-END\n");

	puts("Product of Householder Matrices i.e. QT: \n"); matrix_show2((Double*)shm.H[0]+1, IN_ROWS, IN_ROWS);
	puts("Q=\n"); matrix_show2((Double*)shm.Q+1, IN_ROWS, IN_ROWS);
	puts("R=\n"); matrix_show2((Double*)shm.R+1, IN_ROWS, IN_COLS);
	puts("A=\n"); matrix_show2((Double*)shm.A+1, IN_ROWS, IN_COLS);// Q*R
	
	for(unsigned int t=0; t< NB_STAGES; t++){
	printf("stage-%u cycles: %llu\n", t, shm.stage_cycles[t]);//IMP
	}
	printf("shm.flag:%f\n", shm.flag); //testing
	printf("shm.stage:%f\n", shm.stage); //testing
	printf("shm.test:%f\n", shm.test); //testing
	printf("Size of shm:%llu\n", (unsigned long long) (sizeof shm)); //testing
	
	/* free shm buffer, close workgroup and finalize */
	if(e_free(&emem) != E_OK) FAIL("Can't free!\n");
	if(e_close(&dev) != E_OK) FAIL("Can't close!\n");
	if(e_finalize() != E_OK)  FAIL("Can't finalize!\n");

	printf("Program finished.\n");
	return(0);
}
