/*
 * rank1up.c
 *
 *  Created on: 8 oct. 2019
 *      Author: anmt
 */
#include "ARAIM_RINEX.h"

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            rank1up.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to update the covariance matrix C = inv( G' * W * G ) when a row (i.e., a satellites) is dropped from the geometry matrix G (saving the computational cost of matrix inversion)
 *
 * INPUT PARAMETERS:
 *
 * 		- Nx:					number of columns of G
 * 		- Nsat:					number of satellites
 * 		- C[Nx][Nx]:			covariance matrix
 * 		- G[Nsat][Nx]:			geometry matrix
 * 		- wi:					element of weighting matrix corresponding to the satellite to be dropped
 * 		- idx_sat				index of the row (satellite) to be dropped
 *
 * OUTPUT PARAMETERS:
 *
 * 		- C_red:				updated covariance matrix
 */


 void rank1up( int Nx, int Nsat, double C[Nx][Nx], double G[Nsat][Nx], double wi, int idx_sat, double* C_red ){

	// init variables
	int i,j;
	double tmp_num[Nx][Nx], tmp_den = 0, num[Nx][Nx], den, gi[Nx];

	// extract gi
	for( i = 0; i < Nx; i++ ){
		gi[i] = G[idx_sat][i];
	}

	// compute denominator
	for( i = 0; i < Nx; i++ ){
		for( j = 0; j < Nx; j++ ){
			tmp_den += C[i][j] * gi[i] * gi[j] * wi;
		}
	}

	den = 1 - tmp_den;

	// compute numerator
	for( i = 0; i < Nx; i++ ){
		for( j = 0; j < Nx; j++ ){
			tmp_num[i][j] = wi * gi[i] * gi[j] / den;
		}
	}

	double C_tmp[Nx][Nx]; matrix_mul( Nx , Nx, Nx , C, tmp_num , &C_tmp[0][0] ); // C_tmp = C * (gi*wi*gi');
	matrix_mul( Nx , Nx, Nx , C_tmp, C , &num[0][0] ); // C_num = C * (gi*wi*gi') * C;

	// update
	matrix_sum( Nx, Nx, C, num, C_red );

}
