/*
 * matrix_diff.c
 *
 *  Created on: 7 ago. 2019
 *      Author: anmt
 */

#include "ARAIM_RINEX.h"

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            matrix_diff.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to compute the difference between two matrixes:
 *
 * INPUT PARAMETERS:
 *
 * 		- n_rows:				number of matrix rows
 * 		- n_cols:				number of matrix columns
 * 		- A[n_rows][n_cols]:	first matrix
 * 		- B[n_rows][n_cols]:	second matrix
 *
 *
 * OUTPUT PARAMETERS:
 *
 * 		- C:					output matrix C = A - B
 */

void matrix_diff( int n_rows, int n_cols, double A[n_rows][n_cols], double B[n_rows][n_cols], double* C ){

	int i,j,idx=0;

	for( i = 0; i < n_rows; i++ ){

		for( j = 0; j < n_cols; j++ ){

			C[idx] = A[i][j]-B[i][j];
			idx++;

		}
	}

}

