/*************************************************************
 * Host side program
 * QR factorization using Householder Transformation on 16x16 matrix
 * Developed by Vinay Sawant
**************************************************************/

#include <stdio.h>
#include <math.h>

#define Double float //double//
#define IN_ROWS 8//row size of input matrix
#define IN_COLS 8//column size of input matrix

Double in[IN_ROWS][IN_COLS] = {
	{5.124575, 7.898948, 7.946213, 3.592056, 6.081988, 7.166724, 7.788865, 9.798719}, {2.522771,1.012173,8.776333, 6.716550, 4.709401, 1.210326, 1.280568, 3.150382}, {8.535061, 1.283802, 0.308171, 1.881512, 6.372757, 4.749593, 5.185006, 8.183881}, {3.747190, 1.811495, 6.834885, 9.154308, 5.894849, 0.924535, 4.470197, 0.745557}, {6.461459, 9.131959, 7.072122, 0.483971, 1.670396, 8.949241, 4.018287, 3.546607}, {9.260374, 6.927009, 4.018142, 0.385422, 1.099953, 3.724873, 9.097155, 9.176734}, {6.428624, 0.599322, 6.768367, 8.416201, 2.400816, 6.837505, 0.179495, 4.902166}, {8.330168, 8.243491, 5.997973, 8.827311, 5.079985, 6.432596, 3.651846, 9.114349}
};

Double I[IN_ROWS][IN_ROWS] = {
	{1,0,0,0,0,0,0,0},{0,1,0,0,0,0,0,0},{0,0,1,0,0,0,0,0},{0,0,0,1,0,0,0,0},{0,0,0,0,1,0,0,0},{0,0,0,0,0,1,0,0},{0,0,0,0,0,0,1,0},{0,0,0,0,0,0,0,1}
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

Double vnorm(Double x[], int n)
{
	Double sum = 0;
	for (int i = 0; i < n; i++) sum += x[i] * x[i];
	return sqrt(sum);
}

void vmadd(Double a[], Double b[], Double s, Double c[], int n)
{
	for (int i = 0; i < n; i++)
		c[i] = a[i] + s * b[i];
	return;
}

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


typedef struct{
	Double H[IN_ROWS][IN_ROWS][IN_COLS];
	Double Q[IN_ROWS][IN_ROWS]; // orthogonal matrix
	Double R[IN_ROWS][IN_COLS]; // upper triangular matrix
	Double A[IN_ROWS][IN_COLS]; // result of Q.R (i.e. input matrix)
}shm_t;
shm_t shm;

Double z[IN_ROWS][IN_COLS], z1[IN_ROWS][IN_COLS];
Double q[IN_ROWS][IN_COLS];

int main1()
{
    printf("\n================================hhP1==================================\n");
	//Double z[IN_ROWS][IN_COLS], z1[IN_ROWS][IN_COLS];
    //z=m
    matrix_transfer2((Double*) in, (Double*)z, IN_ROWS, IN_COLS);
//    matrix_show2((Double*)z, IN_ROWS, IN_COLS);
    
    int k=0, l;
	for (k = 0; k < IN_COLS && k < IN_ROWS-1; k++) 
	{
printf("\n----%d-th iteration-----\n", k);
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
	
	// Store kth householder matrix into shared memory
	matrix_transfer2((Double*)q, (Double*)shm.H[k], IN_ROWS, IN_COLS);
//	printf("H[%d]:\n", k); matrix_show2((Double*)shm.H[k], IN_ROWS, IN_COLS);
	
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

int main2(unsigned int last_index)
{
    if(last_index % 2 == 0){
    matrix_transfer2((Double*) I, (Double*)shm.H[last_index], IN_ROWS, IN_ROWS);
    printf("H[%u]:\n", last_index); matrix_show2((Double*)shm.H[last_index], IN_ROWS, IN_ROWS);
    }
    int l=0;
    for(; l<last_index;)
    {
        Double tempH[IN_ROWS][IN_ROWS];
        matrix_mul((Double*)shm.H[l+1], (Double*)shm.H[l], (Double*)tempH, IN_ROWS, IN_ROWS, IN_ROWS, IN_ROWS);
        matrix_transfer2((Double*)tempH, (Double*)shm.H[l/2], IN_ROWS, IN_ROWS);
        l += 2;
    }
    return 0;    
}

int main()
{
    // Stage 1: Calculation of Householder Matrices
    main1();

    // Stage 2: Calculation of Qtranspose
    unsigned int mul_count = IN_ROWS/2;
    while(mul_count != 0)
    {
        main2(2*mul_count-1);
        mul_count /= 2;
    }
    
    printf("H[%u]:\n", mul_count); matrix_show2((Double*)shm.H[mul_count], IN_ROWS, IN_ROWS); //here: mul_count=0
    
    // Stage 3: Calculation of R matrice
    matrix_mul((Double*)shm.H[0], (Double*)in, (Double*)shm.R, IN_ROWS, IN_ROWS, IN_ROWS, IN_COLS);
    puts("R:\n"); matrix_show2((Double*)shm.R, IN_ROWS, IN_COLS); 
    
    // Stage 4: Calculation of Q matrice
    matrix_transfer2((Double*) shm.H[0], (Double*)shm.Q, IN_ROWS, IN_ROWS);
    matrix_transpose((Double*)shm.Q, IN_ROWS, IN_ROWS);
    puts("Q:\n"); matrix_show2((Double*)shm.Q, IN_ROWS, IN_ROWS); 
    
    // Stage 5: Calculation of A matrice i.e. Q*R
    matrix_mul((Double*)shm.Q, (Double*)shm.R, (Double*)shm.A, IN_ROWS, IN_ROWS, IN_ROWS, IN_COLS);
    puts("A:\n"); matrix_show2((Double*)shm.A, IN_ROWS, IN_COLS); 

    return 0;
}