/*
 * GICSRxRINEX.h
 *
 *  Created on: 19/01/2015
 *      Author: cmvv
 */

#ifndef GICSRXRINEX_H_
#define GICSRXRINEX_H_

#include <stdio.h>

#include "dbtypedef.h"

typedef enum
{
	RINEX_v2 = 0,
	RINEX_v3 = 1

} rinex_version_t;

typedef struct
{
	rinex_version_t obs_version;
	rinex_version_t navgps_version;
	rinex_version_t navglo_version;
	rinex_version_t navgal_version;
	rinex_version_t navbei_version;

	char   mixed_navigation;

	long   first_week;
	double first_tow;
	int    leap_seconds;

	double epochs_interval;

	double XYZ_antenna  [3];
	double delta_antenna[3];

	int  number_of_observations[NUM_SYS];
	char observation_list[NUM_SYS][40][3];

} rinex_header_t;

#define MAXLEAPS 18
static const double BDT0 []={2006,1,1,0,0,0}; /* beidou time reference */
static const double LEAPS[MAXLEAPS][7]={ /* leap seconds (y,m,d,h,m,s,gpst-utc) */
    {2017,1,1,0,0,0,18},
    {2015,7,1,0,0,0,17},
    {2012,7,1,0,0,0,16},
    {2009,1,1,0,0,0,15},
    {2006,1,1,0,0,0,14},
    {1999,1,1,0,0,0,13},
    {1997,7,1,0,0,0,12},
    {1996,1,1,0,0,0,11},
    {1994,7,1,0,0,0,10},
    {1993,7,1,0,0,0, 9},
    {1992,7,1,0,0,0, 8},
    {1991,1,1,0,0,0, 7},
    {1990,1,1,0,0,0, 6},
    {1988,1,1,0,0,0, 5},
    {1985,7,1,0,0,0, 4},
    {1983,7,1,0,0,0, 3},
    {1982,7,1,0,0,0, 2},
    {1981,7,1,0,0,0, 1},
};

void get_rinex_info(rinex_header_t *info);

char read_nav_header_RINEX(FILE *fp_nav, gnssdata_t *gnssdata);
void read_obs_header_RINEX(FILE *fp_obs, gnssdata_t *gnssdata, lsq_state_t *lsq_state);

void initialize_ephemeris  (FILE *fp_ngps, FILE *fp_nglo, FILE *fp_ngal, FILE *fp_nbei, gnssdata_t *gnssdata);
//void check_ephemeris_status(FILE *fp_ngps, FILE *fp_nglo, FILE *fp_ngal, FILE *fp_nbei, gnssdata_t *gnssdata);
void check_ephemeris_status(FILE *fp_ngps, FILE *fp_nglo, FILE *fp_ngal, FILE *fp_nbei, gnssdata_t *gnssdata,
		int sat_new_in, int epoch_update_eph, int new_sat[MAX_CHANNELS_G]);

void read_next_obs_RINEX_v2(FILE *fp_obs, gnssdata_t *gnssdata);
void read_next_obs_RINEX_v3(FILE *fp_obs, gnssdata_t *gnssdata);

#endif /* GICSRXRINEX_H_ */
