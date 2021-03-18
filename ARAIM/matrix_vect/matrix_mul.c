/*
 * matrix_mul.c
 *
 *  Created on: 6 ago. 2019
 *      Author: anmt
 */

#include "ARAIM_RINEX.h"

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            matrix_mul.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to compute the difference between two matrixes:
 *
 * INPUT PARAMETERS:
 *
 * 		- P:				number of rows of first matrix
 * 		- Q:				number of columns of first matrix (same as number of rows of second matrix)
 * 		- R:				number of columns of second matrix
 * 		- A[P][Q]:			first matrix
 * 		- B[Q][R]:			second matrix
 *
 *
 * OUTPUT PARAMETERS:
 *
 * 		- C1:				output matrix C = A * B
 */

void matrix_mul( int P, int Q, int R, double A[P][Q], double B[Q][R], double* C1){

	int k=0, p, q, r;

	for (p = 0; p < P; p++ ){

		for (r = 0; r < R; r++){

			C1[k]=0;

			for (q = 0; q < Q; q++){

				C1[k] = C1[k] + A[p][q]*B[q][r];
			}

			k++;

		}

	}

}





