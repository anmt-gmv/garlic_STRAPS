/*
 * transpose.c
 *
 *  Created on: 6 ago. 2019
 *      Author: anmt
 */

#include "ARAIM_RINEX.h"

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            transpose.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to transpose a given input matrix
 *
 * INPUT PARAMETERS:
 *
 * 		- n_rows:					number of matrix rows
 * 		- n_cols:					number of matrix columns
 * 		- A[n_rows][n_cols]:		input matrix
 *
 * OUTPUT PARAMETERS:
 *
 * 		- At:						transpose input matrix At = A';
 */

void transpose(int n_rows, int n_cols, double A[n_rows][n_cols], double* At){

	int i, j, idx = 0;
	for( i = 0; i < n_cols; i++){

		for( j = 0; j < n_rows; j++ ){

			At[ idx ] = A[ j ][ i ];
			idx++;

		}

	}

}
