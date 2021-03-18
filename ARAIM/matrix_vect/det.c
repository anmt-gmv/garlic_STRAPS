/*
 * det.c
 *
 *  Created on: 1 ago. 2019
 *      Author: anmt
 */

#include "ARAIM_RINEX.h"

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            det.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function (recursive) to compute the determinant of a matrix A:
 *
 * INPUT PARAMETERS:
 *
 * 		- Nrows:				number of matrix rows
 * 		- Ncols:				number of matrix columns
 * 		- A[Nrows][Ncols]:	    input matrix
 *
 *
 * OUTPUT PARAMETERS:
 *
 * 		- res:					determinant of matrix A
 */

double det( int Nrows, int Ncols, double A[Nrows][Ncols] ){

	if(Nrows != Ncols){
		printf("\nERROR: matrix needs to be a square matrix to compute the determinant!\n");
		return 0;
	}

	else{

		int j;
		double res = 0;

		/* if the matrix is 2x2, then directly copmute the determinant */
		if ( Nrows == 2 ){

			res = A[0][0] * A[1][1] - A[0][1] * A[1][0];

		}

		/* otherwise start recursion */
		else{

			for( j = 0; j < Nrows; j++ ){

				double A_min[Nrows-1][Ncols-1]; minor_matrix( Nrows, Ncols, A, 0, j, &A_min[0][0] );
				res = res + pow( -1, j ) * A[0][j] * det( Nrows-1, Ncols-1 , A_min );

			}

		}

	return res;

	}

}
