/*
 * minor_matrix.c
 *
 *  Created on: 1 ago. 2019
 *      Author: anmt
 */

#include "ARAIM_RINEX.h"

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            minor_matrix.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to compute the minor of matrix A given row idx_r and column idx_c (same matrix A but w/o row i-th and column j-th):
 *
 * INPUT PARAMETERS:
 *
 * 		- Nrows:				number of matrix rows
 * 		- Ncols:				number of matrix columns
 * 		- A[Nrows][Ncols]:		input matrix
 * 		- idx_r:				index of row to be dropped
 * 		- idx_c:				index of column to be dropped
 *
 * OUTPUT PARAMETERS:
 *
 * 		- p:					minor of matrix A given i,j
 */

void minor_matrix( int Nrows, int Ncols, double A[Nrows][Ncols], int idx_r, int idx_c, double *p ){

	int i,j, idx = 0;

	for( i = 0; i < Nrows; i++ ){

		for( j = 0; j < Ncols; j++ ){

			if( i != idx_r && j != idx_c ){

				p[idx] = A[i][j];
				idx++;

			}

		}

	}

}
