/*
 * float_to_double.c
 *
 *  Created on: 3 sept. 2019
 *      Author: anmt
 */

#include "ARAIM_RINEX.h"

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            float_to_double.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to converter the content of a 2D array of floats to double (the initial number of columns of the input matrix is init_cols):
 *
 * INPUT PARAMETERS:
 *
 * 		- init_cols:			initial number of columns of input matrix (in the case the matrix is pre-initialize to a number of cols > than the actual)
 * 		- init_rows:			initial number of rows of input matrix (in the case the matrix is pre-initialize to a number of rows > than the actual)
 * 		- n_rows:				number of rows (output)
 * 		- n_cols:				number of columns (output)
 * 		- A[n_rows][n_cols]:	matrix to be converted (float)
 *
 *
 * OUTPUT PARAMETERS:
 *
 * 		- A_db:					output matrix (double)
 */

void float_to_double( int init_rows, int init_cols,  int n_rows, int n_cols, float A[init_rows][init_cols], double* A_db ){

	int i , j_in, j_out, idx = 0;

	for( i = 0; i < n_rows; i++ ){

		j_out = 0;

		for( j_in = 0; j_in < init_cols; j_in++ ){

			if( j_in < n_cols ){

				A_db[idx] = (double)A[i][j_out];
				idx++;
				j_out++;

			}

		}

	}

}




