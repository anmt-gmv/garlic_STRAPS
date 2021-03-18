/*
 * ECEF_to_ENU.c
 *
 *  Created on: 11 oct. 2019
 *      Author: anmt
 */

#include "ARAIM_RINEX.h"

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            ECEF_to_ENU.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to convert the geometry matrix G from ECEF to ENU local frame
 *
 * INPUT PARAMETERS:
 *
 * 		- Nsat:							total number of satellites
 * 		- Nx:		 				    number of parameters (number of unknowns)
 * 		- PRN[Nsat]:  					list of PRNs of the satellites in view
 * 		- G_ECEF[Nsat][Nx]: 			geometry matrix in ECEF coordinates
 * 		- usrpos[3]:    				user positions for the conversion in ENU coordinates
 *
 * OUTPUT PARAMETERS:
 *
 * 		- G_ENU:			  			geometry matrix in ENU-local frame
 */

void ECEF_to_ENU( int Nsat, int Nx, int PRN[Nsat], double G_ECEF[Nsat][Nx], double usrpos[3], double* G_ENU ){

	int i, idx = 0;
	int Nconst = Nx - 3;

	// latitude, longitud en height of the user position at the output os LS algorithm
	double la = usrpos[0];
	double lo = usrpos[1];

	// transformation matrix
	double R[3][3];

	R[0][0] = -sin( lo ); R[0][1] = -cos( lo ) * sin( la ); R[0][2] = cos( lo ) * cos( la );
	R[1][0] = cos( lo ); R[1][1] = - sin( lo ) * sin( la ); R[1][2] = sin( lo ) * cos( la );
	R[2][0] = 0; R[2][1] = cos( la ); R[2][2] = sin( la );

	for( i = 0; i < Nsat; i++ ){

		G_ENU[idx++] = G_ECEF[i][0] * R[0][0] + G_ECEF[i][1] * R[1][0] + G_ECEF[i][2] * R[2][0];
		G_ENU[idx++] = G_ECEF[i][0] * R[0][1] + G_ECEF[i][1] * R[1][1] + G_ECEF[i][2] * R[2][1];
		G_ENU[idx++] = G_ECEF[i][0] * R[0][2] + G_ECEF[i][1] * R[1][2] + G_ECEF[i][2] * R[2][2];

		if( Nconst == 2 ){

			/* for GPS add 1 to the 4th column and zero to the 5th*/
			if( PRN[i] <= MAXNUMGPS_G ){

				G_ENU[idx++] = 1;
				G_ENU[idx++] = 0;

			}
			else{

				G_ENU[idx++] = 1;
				G_ENU[idx++] = 1;

			}

		}
		else if( Nconst == 3 ){

			/* for GPS add 1 to the 4th column and zero to the 5th*/
			if( PRN[i] <= MAXNUMGPS_G ){

				G_ENU[idx++] = 1;
				G_ENU[idx++] = 0;
				G_ENU[idx++] = 0;

			}
			else if( PRN[i] > MAXNUMGPS_G && PRN[i] <= MAXNUMGPS_G + MAXNUMGAL_G ){ // In case GALILEO add 1s accordingly

				G_ENU[idx++] = 1;
				G_ENU[idx++] = 1;
				G_ENU[idx++] = 0;

			}
			else{

				G_ENU[idx++] = 1;
				G_ENU[idx++] = 0;
				G_ENU[idx++] = 1;

			}
		}
		else if( Nconst == 4 ){

			/* for GPS add 1 to the 4th column and zero to the 5th*/
			if( PRN[i] <= MAXNUMGPS_G ){

				G_ENU[idx++] = 1;
				G_ENU[idx++] = 0;
				G_ENU[idx++] = 0;
				G_ENU[idx++] = 0;

			}
			else if( PRN[i] > MAXNUMGPS_G && PRN[i] <= MAXNUMGPS_G + MAXNUMGAL_G ){ // In case GALILEO add 1s accordingly

				G_ENU[idx++] = 1;
				G_ENU[idx++] = 1;
				G_ENU[idx++] = 0;
				G_ENU[idx++] = 0;

			}
			else if( PRN[i] > MAXNUMGPS_G + MAXNUMGAL_G && PRN[i] <= MAXNUMGPS_G + MAXNUMGAL_G + MAXNUMGLO_G ){

				G_ENU[idx++] = 1;
				G_ENU[idx++] = 0;
				G_ENU[idx++] = 1;
				G_ENU[idx++] = 0;

			}
			else{

				G_ENU[idx++] = 1;
				G_ENU[idx++] = 0;
				G_ENU[idx++] = 0;
				G_ENU[idx++] = 1;

			}
		}


	}




}
