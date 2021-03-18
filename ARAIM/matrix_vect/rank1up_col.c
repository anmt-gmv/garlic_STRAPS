/*
 * rank1up_col.c
 *
 *  Created on: 9 oct. 2019
 *      Author: anmt
 */
#include "ARAIM_RINEX.h"

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            rank1up_col.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to update the covariance matrix C = inv( G' * W * G ) when a column (i.e., a whole constellation) is dropped from the geometry matrix G (saving the computational cost of matrix inversion)
 *
 * INPUT PARAMETERS:
 *
 * 		- Nx:					number of columns of G
 * 		- C[Nx][Nx]:			covariance matrix
 * 		- idx_col				index of the column to be dropped
 *
 * OUTPUT PARAMETERS:
 *
 * 		- C_up:					updated covariance matrix
 */

void rank1up_col( int Nx, double C[Nx][Nx], int idx_col, double* C_up ){

	int i,j,j_idx;
	double C_sorted_col[Nx][Nx], C_sorted[Nx][Nx], F11inv[Nx-1][Nx-1], F22inv, u2[Nx-1], u2u2T[Nx-1][Nx-1];

	// move column idx_col at the last column position of matrix
	for( i = 0; i < Nx; i++ ){
		j_idx = 0;
		for( j = 0; j < Nx; j++ ){
			if( j != idx_col ){
				C_sorted_col[i][j_idx++] = C[i][j];
			}
			else{
				C_sorted_col[i][Nx-1] = C[i][j];
			}
		}
	}

	// move row idx_col at the last row position of matrix
	for( i = 0; i < Nx; i++ ){
		j_idx = 0;
		for( j = 0; j < Nx; j++ ){
			if( j != idx_col ){
				C_sorted[j_idx++][i] = C_sorted_col[j][i];
			}
			else{
				C_sorted[Nx-1][i] = C_sorted_col[j][i];
			}
		}
	}

	// build the inv of F11 matrix
	for( i = 0; i < Nx - 1; i++ ){
		for( j = 0; j < Nx - 1; j++ ){
				F11inv[i][j] = C_sorted[i][j];
		}
	}

	// build the inv of F22 element
	F22inv = C_sorted[Nx-1][Nx-1];

	// build vector u2
	for( i = 0; i < Nx - 1; i ++ ){
		u2[i] = - C_sorted[i][Nx-1] / F22inv;
	}

	// build matrix u2u2T
	for( i = 0; i < Nx; i++ ){
		for( j = 0; j < Nx; j++ ){
			u2u2T[i][j] = u2[i] * u2[j] * F22inv;
		}
	}

	// update
	int idx = 0;

	for( i = 0; i < Nx -1; i++ ){
		for( j = 0; j < Nx - 1; j++ ){

			C_up[idx] = F11inv[i][j] - u2u2T[i][j];
			idx++;

		}
		idx++;
	}

}
