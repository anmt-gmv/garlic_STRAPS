/*
 * copy_vect_int.c
 *
 *  Created on: 17 sept. 2019
 *      Author: anmt
 */


/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            copy_vect_int.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to make a copy of a vector of integer variables:
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


void copy_vect_int( int n, int A[n], int* copy ){

	int i;

	for( i = 0; i < n; i++ ){

		copy[i] = A[i];

	}

}
