/*
 * inv.c
 *
 *  Created on: 1 ago. 2019
 *      Author: anmt
 */

#include "ARAIM_RINEX.h"

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            inv.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to perform the inversion of a square matrix by Cramer's rule (ATTENTION!! Highly inefficient!):
 *
 * INPUT PARAMETERS:
 *
 * 		- Nrows:					number of rows
 * 		- Ncols:					number of columns
 * 		- A[Nrows][Ncols]:			input matrix
 *
 * OUTPUT PARAMETERS:
 *
 * 		- p:		            	inverted matrix, inv( A )
 */


int inv( int Nrows, int Ncols, double A[Nrows][Ncols], double *p ){

	if(Nrows != Ncols){
		printf("\nERROR: matrix needs to be a square matrix to be inverted!\n");
		return 0;
	}

	else{

		int i, j;
		double A_cof[Nrows][Ncols]; cofactor_matrix( Nrows, Ncols, A, &A_cof[0][0] );
		double det_A = det( Nrows, Ncols, A );

		for( i = 0; i < Nrows; i++ ){

			for( j = 0; j < Ncols; j++ ){

				*p = A_cof[j][i] / det_A;
				p = p + 1;
			}
		}

	}

	return  0;

}
