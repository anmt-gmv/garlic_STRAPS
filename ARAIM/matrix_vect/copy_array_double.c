/*
 * copy_array_double.c
 *
 *  Created on: 7 ago. 2019
 *      Author: anmt
 */

#include "ARAIM_RINEX.h"

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            copy_array_double.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to make a copy of an array of double variables:
 *
 * INPUT PARAMETERS:
 *
 * 		- n_rows:				number of matrix rows
 * 		- n_cols:				number of matrix columns
 * 		- A[n_rows][n_cols]:	original matrix
 *
 *
 * OUTPUT PARAMETERS:
 *
 * 		- copy:					copy of the original matrix
 */

void copy_array_double( int n_rows, int n_cols, double A[n_rows][n_cols], double* copy ){

	int i,j,idx=0;

	for( i = 0; i < n_rows; i++ ){

		for( j = 0; j < n_cols; j++ ){

			copy[idx] = A[i][j];
			idx++;

		}

	}

}

