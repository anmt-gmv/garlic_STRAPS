/*
 * cofactor_matrix.c
 *
 *  Created on: 1 ago. 2019
 *      Author: anmt
 */

#include "ARAIM_RINEX.h"

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            cofactor_matrix.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to build the matrix of cofactors of a matrix:
 *
 * INPUT PARAMETERS:
 *
 * 		- Nrows:				number of matrix rows
 * 		- Ncols:				number of matrix columns
 * 		- A[Nrows][Ncols]:		input matrix
 *
 * OUTPUT PARAMETERS:
 *
 * 		- p:					matrix of cofactors
 */

void cofactor_matrix( int Nrows, int Ncols, double A[Nrows][Ncols], double *p ){

	/* if it is a 2x2 matrix then compute cofactor directly (inverting the diagonal elements and changing the sign of non-diagonal ones) */
	if( Nrows == 2 && Ncols == 2 ){

		p[0] = A[1][1];
		p[1] = -A[1][0];
		p[2] = -A[0][1];
		p[3] = A[0][0];

	}

	/* otherwise compute the minors and the corresponding determinants */
	else{

		int i, j, idx = 0;

		for( i = 0; i < Nrows; i++ ){

			for( j = 0; j < Ncols; j++ ){

				double A_min[ Nrows - 1 ][ Ncols - 1 ]; minor_matrix( Nrows, Ncols, A, i, j, &A_min[0][0] );

				/* alternate sign */
				p[idx] = pow( -1, i + j ) * det( Nrows-1, Ncols-1, A_min );
				idx++;

			}

		}

	}

}
