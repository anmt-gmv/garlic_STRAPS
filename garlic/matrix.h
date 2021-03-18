#ifndef _MATRIX_H
#define _MATRIX_H

#include <math.h>
#include "GICSRx_defines.h"

#define MATEPSILON 	1.e-5

/* A = B*C, where A, B, C are matrices */
#define matMul_f(A,B,C,n,m,l)  				\
{											\
	int i,j,k; 								\
	float aux;   							\
	for(i = 0; i < n; i++){        			\
		for(j = 0; j < m; j++){      		\
			aux = 0.0;              		\
			for(k = 0; k < l; k++)         	\
				aux += B[i][k] * C[k][j]; 	\
			A[i][j] = aux;					\
		}           						\
	}                               		\
}

/* y = A*x, where A matrix, x,y vectors */
#define matVecMul_f(Y, A, X, n, m)			\
{											\
	int i, j;								\
	float aux;								\
	for(i = 0; i < n; i++){        			\
		aux = 0;                    		\
		for(j = 0; j < m; j++)				\
			aux += A[i][j] * X[j];			\
		Y[i] = aux;               			\
	}      									\
}

/* Calculate the scalar product of A,B vectors */
#define DotProduct(A, B, n, Dprod)			\
{                                  			\
	int i;                     				\
	Dprod = 0;								\
   	for(i = 0; i < n; i++)             		\
	  	Dprod += A[i]*B[i]; 				\
}

/* Transpose a matrix A */
#define Transpose(AT, A, n, m)				\
{											\
   	int i,j;                            	\
   	for(i = 0; i < n; i++)                  \
		for(j = 0; j < m; j++)            	\
			AT[j][i] = A[i][j];    			\
}

/* A = B + C, where A,B,C are vectors */
#define matVecSum(A, B, C, n)				\
{                              				\
   	int i;                     				\
   	for(i = 0; i < n; i++)           		\
		A[i] = B[i] + C[i];  				\
}

/* Returns the norm of A 3-component vector */
#define Vec3Norm(A)	(sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]))

/* This structure stores the required data to define a sparse
 * matrix, it is, a matrix where many of its components are 0.
 * It allows for a faster way of matrix multiplication. */
typedef struct
{
	float values           [KNSTATES][KNSTATES];
	float diagonal         [KNSTATES];
	int   col_number       [KNSTATES][KNSTATES];
	int   num_elem_per_row [KNSTATES];

} sparse_matrix_t;

/* Functions related to matrices mapped in vectors */
double get_sym_element(const float *matrix, int row, int col);
void   set_sym_element(float *matrix, int row, int col, float value);
void   update_sym_element(float *matrix, int row, int col, float value);
void   update_matrix(int dim, const float *K, const float *PH, float *P);

/* Sparse matrices functions */
void reset_sparse_matrix(unsigned int size, sparse_matrix_t *matrix);
void add_element_to_sparse_matrix(int row, int col, float value, sparse_matrix_t *matrix);
void update_matrix_sparse(unsigned int size, const sparse_matrix_t *matrix, float *P_mat);

#if USE_FAST_MECHANIZATION == 1
void update_matrix_sparse_fast(const sparse_matrix_t *matrix, float *P_mat);
#endif

/* This function performs the multiplication between a symmetric matrix A and a vector X: Y = AX */
void matSymVecMul_f(float *Y, const float *A, const float *X, int n);

/* Functions used for matrix inversion */
int  choleskyDecomp  (float *pdMatrix, int siSize);
void triangularSolve (const float *pdMatrix, const float * pdInputVector, int siSize, float *pdSolutionVector);
void triangularSolveT(const float *pdMatrix, const float * pdInputVector, int siSize, float *pdSolutionVector);
void choleskyInv     (float *pdMatrixL, int siSize);

#endif //_MATRIX_H
