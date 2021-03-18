/*
 * copy_vect_double.c
 *
 *  Created on: 17 sept. 2019
 *      Author: anmt
 */

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            copy_vect_double.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to make a copy of a vector of double variables:
 *
 * INPUT PARAMETERS:
 *
 * 		- n:			size of the vector
 * 		- A[n]:		    original vector
 *
 *
 * OUTPUT PARAMETERS:
 *
 * 		- copy:			copy of the original vector
 */

void copy_vect_double( int n, double A[n], double* copy ){

	int i;

	for( i = 0; i < n; i++ ){

		copy[i] = A[i];

	}

}
