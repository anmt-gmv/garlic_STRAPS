/*
 * disp_vect.c
 *
 *  Created on: 7 ago. 2019
 *      Author: anmt
 */

#include "ARAIM_RINEX.h"

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            disp_vect.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to display the content of a vector of double variables:
 *
 * INPUT PARAMETERS:
 *
 * 		- size:				size of the vector
 * 		- a[size]:			vector to be displayed
 *
 *
 * OUTPUT PARAMETERS:
 *
 * 		- none
 */

void disp_vect(int size, double a[size]){

	int i;

	for( i = 0; i < size; i++ ){

		printf("%.22f\n", a[i] );

	}

}
