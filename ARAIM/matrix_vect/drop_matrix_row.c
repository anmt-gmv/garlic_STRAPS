/*
 * drop_matrix_row.c
 *
 *  Created on: 7 ago. 2019
 *      Author: anmt
 */

#include "ARAIM_RINEX.h"

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            drop_matrix_row.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to drop a given row of a matrix:
 *
 * INPUT PARAMETERS:
 *
 * 		- n_rows:					number of rows
 * 		- n_cols:					number of columns
 * 		- A[n_rows][n_cols]:		input matrix
 * 		- row_to_drop:				row to be dropped
 *
 *
 * OUTPUT PARAMETERS:
 *
 * 		- A_drop:		            output matrix without the chosen row
 */

void drop_matrix_row( int n_rows, int n_cols, double A[n_rows][n_cols], int row_to_drop, double* A_drop ){

	int i,j,idx=0;

	for( i = 0; i < n_rows; i++ ){

		for( j = 0; j < n_cols; j++){

			/* save all the matrix elements but those belonging to the "col_to_drop" column */
			if( i != row_to_drop ){
				A_drop[idx] = A[i][j];
				idx++;
			}

		}

	}

}
