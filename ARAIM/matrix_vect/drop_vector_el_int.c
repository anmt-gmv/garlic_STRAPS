/*
 * dopr_vector_el_int.c
 *
 *  Created on: 17 sept. 2019
 *      Author: anmt
 */

#include "ARAIM_RINEX.h"

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            dopr_vector_el_int.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to drop one element from a vector of integers:
 *
 * INPUT PARAMETERS:
 *
 * 		- Nsat:				   vector size
 * 		- pr_res[Nsat]:		   original vector
 * 		- idx_sat:			   index of the element to drop
 *
 * OUTPUT PARAMETERS:
 *
 * 		- pr_res_fde:		   output vector not containing the element idx_sat
 */

void drop_vector_el_int( int Nsat , int pr_res[Nsat], int idx_sat, int* pr_res_fde ){

	int i, idx = 0;

	for( i = 0; i < Nsat; i++ ){

		if( i != idx_sat ){

			pr_res_fde[idx] = pr_res[i];
			idx++;

		}

	}


}



