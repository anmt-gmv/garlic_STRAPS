/*
 * disp_matrix_int.c
 *
 *  Created on: 7 ago. 2019
 *      Author: anmt
 */

#include "ARAIM_RINEX.h"

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            disp_matrix_int.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to display the content of a matrix of integer variables:
 *
 * INPUT PARAMETERS:
 *
 * 		- n_rows:				number of matrix rows
 * 		- n_cols:				number of matrix columns
 * 		- A[n_rows][n_cols]:	matrix to be displayed
 *
 *
 * OUTPUT PARAMETERS:
 *
 * 		- none
 */

void disp_matrix_int( int n_rows, int n_cols, int A[n_rows][n_cols] ){

	int i,j;

	for( i = 0; i < n_rows; i++ ){

		printf("\n");

		for( j = 0; j < n_cols; j++ ){

			printf("%d ", A[i][j]);

		}
	}

}


