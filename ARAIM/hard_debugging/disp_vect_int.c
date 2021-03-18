/*
 * disp_vect_int.c
 *
 *  Created on: 7 ago. 2019
 *      Author: anmt
 */

#include "ARAIM_RINEX.h"

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            disp_vect_int.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to display the content of a vector of integer variables:
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

void disp_vect_int(int size, int a[size]){

	int i;

	for( i = 0; i < size; i++ ){

		printf("%d ", a[i] );

	}

	printf("\n");



}

