/*
 * check_new_sat_in.c
 *
 *  Created on: 16 mar. 2021
 *      Author: anmt
 */
#include <string.h>
#include <dbdefine.h>
#include <dbtypedef.h>

int check_new_sat_in( gnssdata_t *gnssdata, int nsat, int *nsat_previous,
		int sat_iw[MAX_CHANNELS_G], int sat_iw_previous[MAX_CHANNELS_G], int sat_new[MAX_CHANNELS_G] ){

	memcpy( sat_iw , gnssdata->SAT_ID , MAX_CHANNELS_G*4 );
	nsat = gnssdata->noOfChannelsAv;
	int sat_new_in = 0, sat_in_set;
	// check if new sat came in
	for( int i_iw = 0; i_iw < nsat; i_iw++ ){
		sat_in_set = 0;
		for( int j_iw = 0; j_iw < *nsat_previous; j_iw++ ){
			if( sat_iw[i_iw] == sat_iw_previous[j_iw] ){
				sat_in_set = 1;
			}
		}
		if( !sat_in_set ){
			sat_new[sat_new_in++] = sat_iw[i_iw];
		}
	}
	memcpy( sat_iw_previous , sat_iw , MAX_CHANNELS_G*4 );
	*nsat_previous = nsat;

	return sat_new_in;

}
