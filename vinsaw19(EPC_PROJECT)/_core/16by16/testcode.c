/*
 * Core side program
 * QR factorization using Householder Transformation on 16x16 Matrix
 * Developed by Vinay Sawant
 */
 
/* Epiphany test application */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include <e_lib.h>

//#include "testcode.h"//this is commented as mutex and barriers are not used here.

#define IN_ROWS 16 //row size of input matrix
#define IN_COLS 16 //column size of input matrix
//COREUSED macro can be defined in range 1 to 16
#define COREUSED 1//16//8//4//2//1

//#if (COREUSED != 1)
//#define MULTICORE
//#endif

#define NB_STAGES 11 //Number of stages used in implementation (tentative)
#define Double float //double


typedef struct{
	Double flag; // to update each core program completion status
	Double stage; // to monitor different stages
	Double test; // used for testing 
	Double H[IN_ROWS][IN_ROWS][IN_COLS]; // householder matrices
	Double Q[IN_ROWS][IN_ROWS]; // orthogonal matrix
	Double R[IN_ROWS][IN_COLS]; // upper triangular matrix
	Double A[IN_ROWS][IN_COLS]; // result of Q.R (i.e. input matrix)
	unsigned long long total_cycles[COREUSED]; // to store total cycles of each core used
	unsigned long long stage_cycles[NB_STAGES]; // to store stage-wise cycles
}shm_t;
//shm_t shm;

volatile shm_t shm SECTION(".shared_dram");

unsigned int row, col;
unsigned int timer_count;// used for timer0
unsigned long long total_cycles;

unsigned int timer1_count;// used for timer1

void delay()
{
	for(volatile int i = 0; i < 20000; i++)
		for(volatile int j = 0; j < 10000; j++)
			;
}

// Check for timer1's current value  and reset the timer if overflown
void chk_timer1_count()
{
	unsigned long long timer1_clk;
	
	timer1_clk = e_ctimer_get(E_CTIMER_1);
	
	if(timer1_clk <= 0)
	{
		timer1_count++;
		e_ctimer_set(E_CTIMER_1, E_CTIMER_MAX);
		e_ctimer_start(E_CTIMER_1, E_CTIMER_CLK);
	}
	
}

// Initiate timer 1 with value E_CTIMER_MAX
void init_timer1()
{
	timer1_count=0;
	e_ctimer_set(E_CTIMER_1, E_CTIMER_MAX);
	e_ctimer_start(E_CTIMER_1, E_CTIMER_CLK);
}

// Calculate total cycles of timer1
unsigned long long calc_time1()
{
	unsigned long long timer1_clk;
	
	timer1_clk = E_CTIMER_MAX - e_ctimer_get(E_CTIMER_1);
	
	unsigned long long stage_cycles = ((timer1_count*E_CTIMER_MAX)+timer1_clk);

	return stage_cycles;	
}


// Check for timer0's current value  and reset the timer if overflown
void chk_timer_count()
{
	unsigned long long timer_clk;
	
	timer_clk = e_ctimer_get(E_CTIMER_0);
	
	if(timer_clk <= 0)
	{
		timer_count++;
		e_ctimer_set(E_CTIMER_0, E_CTIMER_MAX);
		e_ctimer_start(E_CTIMER_0, E_CTIMER_CLK);
	}
	
}

// Initiate timer0 with value E_CTIMER_MAX
void init_timer()
{
	timer_count=0;
	e_ctimer_set(E_CTIMER_0, E_CTIMER_MAX);
	e_ctimer_start(E_CTIMER_0, E_CTIMER_CLK);
}

// Calculate total cycles of timer0
void calc_time()
{
	unsigned long long timer_clk;
	
	timer_clk = E_CTIMER_MAX - e_ctimer_get(E_CTIMER_0);
	
	total_cycles = ((timer_count*E_CTIMER_MAX)+timer_clk);
	
}

// Zero all elements of given matrix
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

// Transfer each element value from matrix x to matrix y
void matrix_transfer2(Double* x, Double* y, unsigned int m, unsigned int n)
{
	//double y[x->m][x->n]; 
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			*((y+i*n)+j) = *((x+i*n)+j);
			//y[i][j] = x->v[i][j];
	return;
}

// Print given matrix i.e. 2D array
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

// Print given vector i.e. 1D array
void vector_show(Double *v, int vsize)
{
	for (int j = 0; j < vsize; j++) {
		printf(" %8.3f", v[j]);
	}
	printf("\n");
}

// Calculate norm of a given vector 
Double vnorm(Double x[], int n)
{
	Double sum = 0;
	for (int i = 0; i < n; i++) sum += x[i] * x[i];
	return sqrt(sum);
}

// Addition of vectors
void vmadd(Double a[], Double b[], Double s, Double c[], int n)
{
	for (int i = 0; i < n; i++)
		c[i] = a[i] + s * b[i];
	return;
}

// Scaling of vector
void vdiv(Double x[], Double d, Double y[], int n)
{
	for (int i = 0; i < n; i++) 
	    y[i] = x[i] / d;
	return;
}

// Calculate minor of given matrices
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

// Get col_nb'th column of given matrix x and store it into vector v
void getColumn(Double *x, Double *v, unsigned int col_nb)
{
	for (unsigned int i = 0; i < IN_ROWS; i++)
		v[i] = *((x+i*IN_COLS)+col_nb);
	return;
}

// Get a householder matrix
void getHouseholderMatrix(Double v[], Double h[][IN_ROWS])
{
//    matrix_clear((Double*)h, IN_ROWS, IN_ROWS);
	for (int i = 0; i < IN_ROWS; i++)
		for (int j = 0; j < IN_ROWS; j++)
			h[i][j] = -2 * v[i] * v[j];
	for (int i = 0; i < IN_ROWS; i++)
		h[i][i] += 1;
 
	return;
}

// Matrix multiplication of x and y, result is stored in z
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

// Transpose of given matrix and update the same variable
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

Double in[IN_ROWS][IN_COLS] = {
	{5.124575, 7.898948, 7.946213, 3.592056, 6.081988,7.166724, 7.788865,9.798719,4.258254, 1.054377, 7.652938, 3.552302, 2.828004, 1.762015, 2.638748, 8.069129}, {2.522771,1.012173,8.776333, 6.716550,  4.709401,  1.210326, 1.280568, 3.150382, 9.152067, 2.220713, 4.314055, 2.801475, 0.355802, 4.751442, 2.452192, 8.104572}, {8.535061, 1.283802, 0.308171, 1.881512, 6.372757, 4.749593, 5.185006, 8.183881, 0.145697, 9.572844, 4.968377, 2.484793, 3.972274, 0.722540, 4.236591, 8.985795}, {3.747190, 1.811495, 6.834885, 9.154308, 5.894849, 0.924535, 4.470197, 0.745557, 2.168096, 6.291060, 6.742811, 0.373962, 4.722142, 6.529717, 9.824009, 7.356618}, {6.461459, 9.131959, 7.072122, 0.483971, 1.670396, 8.949241, 4.018287, 3.546607, 8.672633, 0.835284, 1.438314, 5.034810, 6.390177, 3.856355, 2.282590, 2.083286}, {9.260374, 6.927009, 4.018142, 0.385422, 1.099953, 3.724873, 9.097155, 9.176734, 6.069993, 8.241838, 9.264992, 1.061532, 3.459950, 2.737167, 2.078934, 3.678259}, {6.428624, 0.599322, 6.768367, 8.416201, 2.400816, 6.837505, 0.179495, 4.902166, 7.376048, 9.049492, 7.510052, 1.206717, 5.326238, 2.406306, 4.195184, 9.811339}, {8.330168, 8.243491, 5.997973, 8.827311, 5.079985, 6.432596, 3.651846, 9.114349, 1.086232, 7.299483, 4.097865, 8.284714, 0.099244, 9.042177, 1.308067, 3.394185}, {6.530458, 2.984305, 6.574254, 5.616988, 9.094666, 5.588702, 5.533086, 7.406141, 4.161861, 0.238815, 9.648564, 6.822000, 0.370360, 4.103983, 8.835959, 5.350417}, {6.909482, 7.122869, 4.538856, 0.413352, 0.563682, 8.213211, 3.641620, 3.405092, 1.848284, 4.093572, 8.746218, 8.198015, 4.431138, 8.115096, 3.752868, 9.309067}, {4.522346, 4.953636, 6.954018, 4.096635, 1.076255, 2.290696, 4.178736, 1.015898, 5.956050, 5.301290, 3.462898, 9.802859, 9.976361, 8.454410, 2.433835, 3.138024}, {6.622229, 8.626863, 4.287552, 5.543465, 5.243287, 7.461806, 2.883427, 3.698015, 4.793037, 1.159368, 1.629081, 2.267589, 0.671444, 9.046292, 4.284649, 2.955473}, {4.077998, 8.902470, 5.945147, 3.556053, 8.246637, 4.732425, 0.047543, 5.104174, 6.272343, 6.639598, 9.528616, 1.437121, 6.459073, 0.140771, 9.660500, 8.103308}, {7.141063, 7.821795, 4.805562, 5.816616, 8.541925, 0.016668, 2.015325, 4.703988, 5.374593, 0.995246, 5.829943, 5.670019, 8.307987, 6.370590, 9.468715, 1.161778}, {7.457226, 0.099300, 7.410219, 0.187684, 1.784756, 1.066501, 8.859654, 6.758221, 9.879946, 3.690862, 4.245235, 3.524249, 2.964446, 4.373680, 2.311288, 2.500651}, {4.752189, 0.592439, 2.946092, 0.001511, 6.350062, 6.068312, 2.986027, 0.368879, 3.349494, 5.183024, 7.386629, 0.150104, 9.933439, 1.343407, 0.192133, 6.214067}
};

Double I[IN_ROWS][IN_ROWS] = {
	{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}
};

Double z[IN_ROWS][IN_COLS], z1[IN_ROWS][IN_COLS];
Double q[IN_ROWS][IN_COLS];

// This function is used to compute all the householder matrices
int main1()
{
Double temp = 0;
//    printf("\n================================hhP1==================================\n");
    //Take input matrix from in and trnsfer it to z
    matrix_transfer2((Double*)in, (Double*)z, IN_ROWS, IN_COLS);
//    matrix_show2((Double*)z, IN_ROWS, IN_COLS);
    
    int k=0, l;
	for (k = 0; k < IN_COLS && k < IN_ROWS-1; k++) 
	{
//printf("\n----%d-th iteration-----\n", k);
	Double e[IN_ROWS], x[IN_ROWS], a;
    
    //matrix_minor
    _minor(z, z1, k);
//    puts("matrix_minor z1"); matrix_show2((Double*)z1, IN_ROWS, IN_COLS);

    //Transfer matrix z1 to matrix z
    matrix_transfer2((Double*) z1, (Double*)z, IN_ROWS, IN_COLS);
//    puts("updated z"); matrix_show2((Double*)z, IN_ROWS, IN_COLS);
    
    //Get column number k of z matrix
    getColumn((Double*)z, x, k);
//    printf("%d-th column of z, put in x\n",k); vector_show(x, IN_ROWS);
    
    //Get l2-norm of vector x
    a = vnorm(x, IN_ROWS);
//    printf("vnorm of x: %5.5f\n", a);
    
    //Sign change if diagonal element of input matrix is positive 
    if (in[k][k] > 0) a = -a;
//    printf("a: %5.5f\n", a); 
    //
    for (int i = 0; i < IN_ROWS; i++)
	{
    	e[i] = (i == k) ? 1 : 0;
	}
//	printf("e=>vector e:\n"); vector_show(e, IN_ROWS);
	
	//e = x + a.e
	vmadd(x, e, a, e, IN_ROWS);
//	printf("u=>vector e:\n"); vector_show(e, IN_ROWS);
	
	//e = e/norm(e)
	vdiv(e, vnorm(e, IN_ROWS), e, IN_ROWS);
//	printf("v=>vector e:\n"); vector_show(e, IN_ROWS);
	
	//Get kth householder matrix in q
	getHouseholderMatrix(e, q);
//	puts("H=>q[k]"); printf("q:\n"); matrix_show2((Double*)q, IN_ROWS, IN_COLS);
//temp = q[0][0];
	// Store kth householder matrix into shared memory
	matrix_transfer2((Double*)q, (Double*)(shm.H[k]), IN_ROWS, IN_COLS);
//	printf("H[%d]:\n", k); matrix_show2((Double*)shm.H[k], IN_ROWS, IN_COLS);
//if(k==0 || k==1) shm.H[k][0][0]=temp;	
	//z1 = q.z
	matrix_mul((Double*)q, (Double*)z, (Double*)z1, IN_ROWS, IN_COLS, IN_ROWS, IN_COLS);
//	puts("Ha=>z1"); matrix_show2((Double*)z1, IN_ROWS, IN_COLS);
	
	//Transfer matrix z1 to matrix z
    matrix_transfer2((Double*) z1, (Double*)z, IN_ROWS, IN_COLS);
//    puts("z"); matrix_show2((Double*)z, IN_ROWS, IN_COLS);
	
	}
    // Store identity matrix I as last element of householder matrix array	
    matrix_transfer2((Double*) I, (Double*)shm.H[k], IN_ROWS, IN_ROWS);
//    matrix_show2((Double*)shm.H[k], IN_ROWS, IN_ROWS);
    return 0;
}

// This function is used to compute product of all householder matrices i.e. Hn-1.Hn-2.....H1.H0
int main2(unsigned int last_index, unsigned int cid)
{
    if(last_index % 2 == 0){
    matrix_transfer2((Double*) I, (Double*)shm.H[last_index], IN_ROWS, IN_ROWS);
//    printf("H[%u]:\n", last_index); matrix_show2((Double*)shm.H[last_index], IN_ROWS, IN_ROWS);
    }
    int l=0;
    for(; l<last_index;)
    {
	//Each core have distinct pair of householder matrices
	if(cid == (l/2)%COREUSED){
        Double tempH[IN_ROWS][IN_ROWS];
        matrix_mul((Double*)shm.H[l+1], (Double*)shm.H[l], (Double*)tempH, IN_ROWS, IN_ROWS, IN_ROWS, IN_ROWS);
//	delay();//
	matrix_transfer2((Double*)tempH, (Double*)shm.H[l/2], IN_ROWS, IN_ROWS);
		}
        l += 2;//increment by 2
	chk_timer1_count();
    }
    return 0;    
}

//MAIN - core program starts from this function - QR decomposition using Householder tranformation
// In this function process runs from stage1 to succesive stages. 
//The execution cycles of each stage is calculated. Total cycles are also calculated; but that shows overflow value in most cases. So, we can consider sum of all stage cycles with delay cycles as total cycles. 
int main(int argc, char *argv[])
{
	init_timer ();

	init_timer1();	
	e_coords_from_coreid(e_get_coreid(), &row, &col);
	unsigned int cid = 4*row + col;
	
	shm.test = 99;
	shm.stage = 1;// initiate stage variable as 1 to start process from main1
	chk_timer1_count();
	shm.stage_cycles[0] = calc_time1();
	delay();
	if(shm.stage == 1 && cid == 0){
	init_timer1();
    // Stage 1: Calculation of Householder Matrices
	main1();
	shm.stage = 2;	
	chk_timer1_count();
	shm.stage_cycles[1] = calc_time1();
	delay();
	}
	while(shm.stage < 2);
	init_timer1();
    // Stage 2: Calculation of Qtranspose [PARALLELIZED]
    	unsigned int mul_count = IN_ROWS/2;
    	while(mul_count != 0)
    	{
       	main2(2*mul_count-1, cid);
       	mul_count /= 2;
   	}
   	shm.stage = 3;	
	chk_timer1_count();
	shm.stage_cycles[2] = calc_time1();
	delay();
	//Add condn for stage 3 onwards
	if(shm.stage == 3 && cid == 0){
	init_timer1();
//    	printf("H[%u]:\n", mul_count); matrix_show2((Double*)shm.H[mul_count], IN_ROWS, IN_ROWS); //here: mul_count=0
    
    // Stage 3: Calculation of R matrice
    	matrix_mul((Double*)shm.H[0], (Double*)in, (Double*)shm.R, IN_ROWS, IN_ROWS, IN_ROWS, IN_COLS);
//      puts("R:\n"); matrix_show2((Double*)shm.R, IN_ROWS, IN_COLS); 
	shm.stage = 4;
	chk_timer1_count();
	shm.stage_cycles[3] = calc_time1();
	delay();
	}
	
	if(shm.stage == 4 && cid == 0){
	init_timer1();
    // Stage 4: Calculation of Q matrice
    	matrix_transfer2((Double*) shm.H[0], (Double*)shm.Q, IN_ROWS, IN_ROWS);
    	matrix_transpose((Double*)shm.Q, IN_ROWS, IN_ROWS);
//    	puts("Q:\n"); matrix_show2((Double*)shm.Q, IN_ROWS, IN_ROWS); 
	shm.stage = 5;
	chk_timer1_count();
	shm.stage_cycles[4] = calc_time1();
	delay();
    }
	//Do stage 5 to see the result of Q*R
	if(shm.stage == 5 && cid == 0){
	init_timer1();
    // Stage 5: Calculation of A matrice i.e. Q*R
    	matrix_mul((Double*)shm.Q, (Double*)shm.R, (Double*)shm.A, IN_ROWS, IN_ROWS, IN_ROWS, IN_COLS);
//    	puts("A:\n"); matrix_show2((Double*)shm.A, IN_ROWS, IN_COLS); 
	shm.stage = 6;
	chk_timer1_count();
	shm.stage_cycles[5] = calc_time1();
	delay();
	}
	
    // Stage 6 is just to see last exiting core
	init_timer1();
	shm.test = cid;
	chk_timer1_count();
	shm.stage_cycles[6] = calc_time1();
	
    // Last stage to check cycles used in delay
	init_timer1();
	delay();
	chk_timer1_count();
	shm.stage_cycles[NB_STAGES-1] = calc_time1();
    // Compute total cycles	
	chk_timer_count();
	calc_time();
	shm.total_cycles[cid] = total_cycles;
	// increment the flag variable to indicate completion of process in cores
	shm.flag += 1; 
	while(1);
    return 0;
}
