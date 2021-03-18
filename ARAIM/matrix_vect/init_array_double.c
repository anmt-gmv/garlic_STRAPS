/*
 * init_array_double.c
 *
 *  Created on: 21 ago. 2019
 *      Author: anmt
 */

#include "ARAIM_RINEX.h"

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            init_array_double.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to initialize to 0 an array of double variables:
 *
 * INPUT PARAMETERS:
 *
 * 		- nrows:				number of matrix rows
 * 		- ncols:				number of matrix columns
 *
 *
 * OUTPUT PARAMETERS:
 *
 * 		- A_init:				array initialized to zero
 */

void init_array_double( int nrows, int ncols, double* A_init ){

	int i, j, idx = 0;

	for( i = 0; i < nrows; i++ ){

		for( j = 0; j < ncols; j++ ){

			A_init[idx] = 0;
			idx++;

		}
	}

}
