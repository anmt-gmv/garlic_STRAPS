/*
 * sum_elements_vect.c
 *
 *  Created on: 8 ago. 2019
 *      Author: anmt
 */

#include "ARAIM_RINEX.h"

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            sum_elements_vect.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to compute the sum of the elements of a given vector
 *
 * INPUT PARAMETERS:
 *
 * 		- size:					size of the vector
 * 		- a[size]:				input vector
 *
 * OUTPUT PARAMETERS:
 *
 * 		- sum:					sum of the elements of vector a
 */

double sum_elements_vect( int size, double a[size] ){

	int i;
	double sum = 0;
	for( i = 0; i < size; i++ ){

		sum += a[i];

	}

	return sum;

}
