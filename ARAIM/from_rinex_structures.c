/*
 * from_rinex_structures.c
 *
 *  Created on: 21 oct. 2019
 *      Author: anmt
 */

#include "dbtypedef.h"
#include "ARAIM_RINEX.h"
#include <algebra.h>

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            from_rinex_structures.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to extract data from the RINEX structures and to store them into vectors/matrices suitable as input for the ARAIM algorithm
 *
 * INPUT PARAMETERS:
 *
 * 		- lsq_state:					structure containing state information after the LS algorithm
 * 		- lsq_mat:		 				structure containing the matrices/vectors output of the LS algorithm
 * 		- gnssdata:						structure containing the data extractd from the RINEX file
 * 		- ARAIM:					    struct containing all the ARAIM configuration parameters (see function "load_sim_parameters.c" )
 * 		- Nsat: 						number of satellites in view
 * 		- Nx: 							number of parameters (unknowns)
 *
 * OUTPUT PARAMETERS:
 *
 * 		- lsq_num_LOS[NUM_SYS]:   	    define the number of staellites in view per each constellation (e.g., num_LOS[0] = 2 means that there are 2 GPS satellites in view)
 * 		- PRN_used:						list of the PRNs of the satellites in view
 * 		- G:							observation matrix in ENU-local frame coordinates
 * 		- w_ura:  					    weighting matrix used for integrity purposes
 * 		- w_ure:   						weighting matrix used for continuity/accuracy purposes
 * 		- c_ure:  					    inverse of the weighting matrix used for continuity/accuracy purposes, e.g., w_ure = inv( c_ure )
 * 		- pr_res:		 			    pseudorange residuals
 */

void from_rinex_structures( lsq_state_t *lsq_state, lsq_matrix_t *lsq_mat, gnssdata_t *gnssdata, configuration *ARAIM,
		int Nsat, int Nx, int* lsq_num_LOS, int* PRN_used, double* G, double* w_ura, double* w_ure, double* c_ure, double* pr_res ){

	int i;

	// number of constellations
	for( i = 0; i < NUM_SYS; i++ ) lsq_num_LOS[i] = lsq_state->num_LOS[i];

	// extract converged position in Lat Long Heigh
	double geod_pos[3]; ECEFtoNAV_pos(lsq_state->pos, geod_pos);

	for( i = 0; i < lsq_state->total_obs; i++ ) PRN_used[i] = lsq_state->prn_used[i];

	// build obs matrix in ENU
	double G_double[Nsat][Nx]; float_to_double( MAX_CHANNELS_G , NUM_PARS, Nsat, Nx , lsq_mat->G , &G_double[0][0] ); /* convert obs matrix to double */

	ECEF_to_ENU( Nsat,  Nx, PRN_used, G_double , geod_pos, G );

	init_array_double( Nsat , Nsat , w_ura );
	init_array_double( Nsat , Nsat , w_ure );
	init_array_double( Nsat , Nsat , c_ure );
	init_array_double( Nsat , 1 , pr_res );

}
