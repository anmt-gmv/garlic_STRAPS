/////////////////////////////////////////////////////////////////////////////////
// Project:                     GICSRX
// Purpose:                     Functions for LSQ Algorithm and IBPL computation.
// File:                        GICSRxWNavSol.c
// Version:                     1.0
// Author:                      GMV, S.A.
// Language:                    c
// Compiler:                    ANSI
// Creation date:               -
// Last change date:            -
// Changes:
// +--------------------------------------------------------------------------+
// | Version |   Date   | Author       | Change reason                        |
// |---------+----------+--------------+--------------------------------------|
// |   1.00  | 13/12/12 | cmvv         | Initial Release                      |
// +--------------------------------------------------------------------------+
//
// DESCRIPTION: this file includes the definitions of the least-squares function,
//              atmospheric correction functions and IBPL integrity functions.
//
//
/////////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"
#include "algebra.h"
#if CALC_ATMSPHCORR == 1
#include "calendar.h"
#endif
#include "GICSRxObs.h"
#include "GICSRxPosition.h"
#include "GICSRxWNavSol.h"
#include "GICSRxIBPL.h"
#include "garlic_functions.h"

#define MAX_ITER_LSQ (first_position ? 5 : 15)
#define DELTA_X_LSQ  10

static char first_position = 0;

/* Extract state information from the navigation data structures */
static void extract_nav_solution(int week, double tow, const lsq_state_t *lsq_state, nav_sol_t *nav_sol, quality_t quality)
{
	int i;
	if(quality >= nav_sol->quality)
	{
		nav_sol->week = week;
		nav_sol->tow  = tow;

		nav_sol->kalman_valid   = 0;
		nav_sol->position[0]    = lsq_state->pos[0];
		nav_sol->position[1]    = lsq_state->pos[1];
		nav_sol->position[2]    = lsq_state->pos[2];
		nav_sol->velocity[0]    = lsq_state->vel[0];
		nav_sol->velocity[1]    = lsq_state->vel[1];
		nav_sol->velocity[2]    = lsq_state->vel[2];

		// In case that there is not clock bias for GPS, it is taken from the the first next
		// system clock valid and the isb for that other constellation is 0.
		int ref = 0, first = 1;
		memset(nav_sol->isb, 0, sizeof(nav_sol->isb));
		for (i=0; i<NUM_SYS; i++)
		{
			if (lsq_state->num_LOS[i] != 0 && first)
			{
				nav_sol->clock_bias = lsq_state->pos[3+i];
				ref = i;
				first = 0;
			}
			if (i>0 && lsq_state->num_LOS[i] != 0)
			{
				nav_sol->isb[i-1] = lsq_state->pos[3+i] - lsq_state->pos[3+ref];
			}
		}


		for (i=0; i<NUM_FREQ-1; i++) { nav_sol->ifreqBias[i] = lsq_state->pos[3+NUM_SYS+i]; }

		nav_sol->clock_drift[0] = lsq_state->vel[3];
		if (TWO_RX == 1) { nav_sol->clock_drift[1] = lsq_state->vel[4]; }

		nav_sol->prev_used_sats = lsq_state->total_obs_prev;
		nav_sol->curr_used_sats = lsq_state->total_obs;

		memcpy(nav_sol->num_LOS, lsq_state->num_LOS, sizeof(nav_sol->num_LOS));

		// Copy mask of used satellites in solution.
		memcpy(nav_sol->prn_used, lsq_state->prn_used, sizeof(lsq_state->prn_used));

		nav_sol->prop_time = 1.0;
		nav_sol->quality   = quality;

		nav_sol->hdop = lsq_state->hdop;
		nav_sol->vdop = lsq_state->vdop;
	}
}

/* This function reconstructs the measurements and compute the associated partial derivatives and residuals.
 * With this purpose, calculate first the satellite state (position and velocity), and updates the ephemeris
 * if necessary. The information calculated for the satellite is used after in the computation of the Kalman Filter.
 */
static void reconstruct_obs(gnssdata_t *gnssdata, SP3STRUCT_D *pSp3, lsq_state_t *lsq_state, lsq_matrix_t *lsq_mat, char prop_mode, propephem_t *propephem, kconf_t *kconf, int count, char *update_eph)
{
	double GeoPos[3];
	float  LatLonTrig[4];

	// Initialization of the geometry matrix.
	memset(lsq_mat->G, 0, sizeof(lsq_mat->G));

	// Some trigonometric operations regarding longitude and latitude
	// are made previously to compute faster azimuth and elevation.
	ECEFtoNAV_pos(lsq_state->pos, GeoPos);

	LatLonTrig[0] = (float)cos(GeoPos[0]);
	LatLonTrig[1] = (float)sin(GeoPos[0]);
	LatLonTrig[2] = (float)cos(GeoPos[1]);
	LatLonTrig[3] = (float)sin(GeoPos[1]);

	int imonth = 0;
	double day = 0;
#if CALC_ATMSPHCORR == 1
	int iyear, iday, ihour, iminute;
	double second;

	GPSToCal_G(gnssdata->week, gnssdata->tow, &iyear, &imonth, &iday, &ihour, &iminute, &second);
	day = DateToDay(iyear, imonth, iday);
#endif

	// Iterate for each valid GPS/GLONASS satellite in prn_visib_sats.
	int i;
	for(i = 0; i < gnssdata->noOfChannelsAv; i++)
	{
		// Obtain for the current iteration the value of the innovations / residuals for Least-Squares algorithm.
		reconstruct_observation(imonth, day, gnssdata->tow, kconf, count, gnssdata->EPH, pSp3, prop_mode, propephem, LatLonTrig, lsq_state,
								&gnssdata->OBS[i], lsq_mat->G[i], &lsq_state->scaled_bias[i], update_eph, &gnssdata->iono_utc);
	}

} // END of function reconstruct_obs

/* This function calculates the geometry matrix for each available satellite (GPS and GLONASS)
 * in the channel mask. It builds also the weighting matrix.
 */
static int build_obs_matrices(gnssdata_t *gnssdata, lsq_state_t *lsq_state, lsq_matrix_t *lsq_mat, kconf_t *kconf)
{
	memset(lsq_state->num_LOS, 0, sizeof(lsq_state->num_LOS));
	lsq_state->total_LOS = 0;
	memset(lsq_state->used_freq,0, sizeof(lsq_state->used_freq));
	memset(lsq_state->prn_used, 0, sizeof(lsq_state->prn_used));
	lsq_state->total_obs = 0;

	int n = 0, gal_ref = 0, gps_ref = 0;

	// Iterate for each valid GPS/GLONASS satellite in prn_visib_sats.
	int i;
	for(i = 0; i < gnssdata->noOfChannelsAv; i++)
	{

		if (gnssdata->OBS[i].prange_status == OBS_VALID)
		{
			// Set weighting matrix.
#if USE_WEIGHTING_MODEL == 1
			prange_weighting_model (&gnssdata->OBS[i], 10.0, kconf->kvelStopThr, kconf->kvelStopFactor);
			doppler_weighting_model(&gnssdata->OBS[i]);
#else
			gnssdata->OBS[i].prange_sigma  = kconf->sigma_prange;
			gnssdata->OBS[i].doppler_sigma = kconf->sigma_doppler;
#endif
			lsq_mat->Wpos[n] = 1.0 / (gnssdata->OBS[i].prange_sigma  * gnssdata->OBS[i].prange_sigma);
			lsq_mat->Wvel[n] = 1.0 / (gnssdata->OBS[i].doppler_sigma * gnssdata->OBS[i].doppler_sigma);

			// Residual range to be minimized is calculated.
			lsq_mat->zpos[n] = gnssdata->OBS[i].prange_residual;
			lsq_mat->zvel[n] = gnssdata->OBS[i].doppler_residual;

			int k;
			for (k = 0; k < NUM_PARS; k++)
			{
				lsq_mat->G[n][k] = lsq_mat->G[i][k];
			}

			// Increase the number of valid satellites (PRNs)
			lsq_state->prn_used[n++] = gnssdata->OBS[i].PRN;

			char repeated = 0;
			for (k = n-2; k >= 0; k--)
			{
				if (lsq_state->prn_used[n-1] == lsq_state->prn_used[k])	{
					repeated = 1;
					break;
				}
			}

			int sys = signal2system(gnssdata->OBS[i].sigFlag);
			int freq = signal2freq(gnssdata->OBS[i].sigFlag);

			// Increase the number of lines-of-sight per constellation only once
			if (!repeated)
			{
				if (sys >= 0) { lsq_state->num_LOS[sys]++; }
				lsq_state->total_LOS++;
			}

			// Increase the number of measurement used per frequency band
			if (freq >= 0) { lsq_state->used_freq[freq]++; }

			// Check reference signal (E1/L1) for GPS and GAL are present
			if (sys == SYS_GAL && freq == FREQ_L1E1) { gal_ref = 1; }
			if (sys == SYS_GPS && freq == FREQ_L1E1) { gps_ref = 1; }
		}
	} // END for loop

	// Number of valid measurements without considering dummy observations
	lsq_state->total_obs = n;

	/* For multi-frequency constellation (GPS and GAL) it is checked if there are
	 * 	measurements for the reference frequency, which is L1C/A for GPS and E1 for GAL
	 * 	In case of missing the reference frequency, the next one with measurements
	 * 	available is considered as reference. */
	if (gal_ref == 0)
	{
		if (lsq_state->used_freq[FREQ_E5a] > 0)
		{
			for (i = 0; i < n; i++)
			{
				lsq_mat->G[i][3+NUM_SYS+FREQ_E5a-1] = 0;
				lsq_state->used_freq[FREQ_E5a] = 0;
			}
		}
		else if (lsq_state->used_freq[FREQ_E5b] > 0)
		{
			for (i = 0; i < n; i++)
			{
				lsq_mat->G[i][3+NUM_SYS+FREQ_E5b-1] = 0;
				lsq_state->used_freq[FREQ_E5b] = 0;
			}
		}
	}

	if (gps_ref == 0)
	{
		if (lsq_state->used_freq[FREQ_L5] > 0)
		{
			for (i = 0; i < n; i++)
			{
				lsq_mat->G[i][3+NUM_SYS+FREQ_L5-1] = 0;
				lsq_state->used_freq[FREQ_L5] = 0;
			}
		}
	}


	/* Since matrixes are dimensioned to support the maximum number of constellations (NUM_SYS)
	 * and frequencies (NUM_FREQ) in case of missing measurements for any of them, to avoid
	 * singularity of the matrix, it is added dummy measurements */
	int num_dummy_obs = 0;
	// Add dummy observation if case of missing constellations to avoid singular geometry matrix.
	for (i=0; i<NUM_SYS; i++)
	{
		if (lsq_state->num_LOS[i] == 0)
		{
			memset(lsq_mat->G[n],0,sizeof(lsq_mat->G[n]));
			lsq_mat->G[n][3+i] 		 = 1;
			lsq_mat->Wpos[n]   		 = 1;
			lsq_mat->zpos[n]   		 = 0;
			lsq_state->prn_used[n++] = -1;
			num_dummy_obs++;
		}
	}

	// Add dummy observation if case of missing frequency to avoid singular geometry matrix.
	for (i=1; i<NUM_FREQ; i++)
	{
		if (lsq_state->used_freq[i] == 0)
		{
			memset(lsq_mat->G[n],0,sizeof(lsq_mat->G[n]));
			lsq_mat->G[n][3+NUM_SYS+i-1] = 1;
			lsq_mat->Wpos[n]   		     = 1;
			lsq_mat->zpos[n]   		     = 0;
			lsq_state->prn_used[n++]     = -1;
			num_dummy_obs++;
		}
	}

	return num_dummy_obs;

} // END of function build_obs_matrices

#if CHECK_OUTLIERS == 1
/*
 * This function performs simple checks over the pseudo-range and doppler residuals, and rejects the measurements
 * if the deviation with respect to the typical noise is higher than a certain threshold.
 */
static void check_outliers(gnssdata_t *gnssdata, int count, double *init_pos)
{
	int i;
	int num_prange  = 0;
	int num_doppler = 0;

	float median_prange = 0, median_doppler = 0, mean_prange = 0;
	float list_range[MAX_CHANNELS_G], list_doppler[MAX_CHANNELS_G];

	const float threshold_prange  = 50;
	const float threshold_doppler = 10;

	if(first_position == 0)
	{
		return;
	}

	float factor = ((count == 0) ? 10.0 : 1.0);

	for(i = 0; i < gnssdata->noOfChannelsAv; i++)
	{
		// Pseudo-range check.
		if (gnssdata->OBS[i].prange_status == OBS_VALID)
		{
			list_range[num_prange++] = gnssdata->OBS[i].prange_residual;
			mean_prange             += gnssdata->OBS[i].prange_residual;
		}
		// Doppler check.
		if (gnssdata->OBS[i].doppler_status == OBS_VALID)
		{
			list_doppler[num_doppler++] = gnssdata->OBS[i].doppler_residual;
		}
	}

	compute_median(list_range,   num_prange,  &median_prange);
	compute_median(list_doppler, num_doppler, &median_doppler);

	// Iterate for each available channel.
	for(i = 0; i < gnssdata->noOfChannelsAv; i++)
	{
		// Check, when the pseudo-range observation is a valid measurement, that the residual is not higher than a threshold.
		// If it is higher, discard current measurement. Notice that it is required that at least two iterations of the
		// Least-Squares algorithm have been performed. In this loop the residual rms value for each type of observation
		// is computed.
		if (gnssdata->OBS[i].prange_status == OBS_VALID)
		{
			if (fabs(gnssdata->OBS[i].prange_residual - median_prange) > factor*threshold_prange)
			{
				gnssdata->OBS[i].prange_status = REJ_CHISQUARE_TEST;
			}
		}

		// Check, when the doppler observation is a valid measurement, that the residual is not higher than a threshold.
		// If it is higher, discard current measurement. Notice that it is required that at least two iterations of the
		// Least-Squares algorithm have been performed.
		if (gnssdata->OBS[i].doppler_status == OBS_VALID)
		{
			if (fabs(gnssdata->OBS[i].doppler_residual - median_doppler) > factor*threshold_doppler)
			{
				gnssdata->OBS[i].doppler_status = REJ_CHISQUARE_TEST;
			}
		}
	}
} // END of function check_outliers

#endif

/* This function calculates the module of the pseudo-range and doppler residuals vectors, properly weighted.
 * These values are required to calculate the Least-Squares IBPL integrity solution.
 */
static void compute_lsq_indicators(lsq_state_t *lsq_state, lsq_matrix_t *lsq_mat, char num_obs, char result, int num_dummy_obs_ibpl )
{
	int i;

	// Store current value to be considered in the next epoch.
	lsq_state->total_obs_prev = lsq_state->total_obs;

	// Initialization of statistics for pseudo-range.
	lsq_state->num_pos_dof = num_obs - lsq_state->num_pos_pars + num_dummy_obs_ibpl;
	lsq_state->prange_residuals = 0.0;

	// Calculate the weighted squared module of the pseudo-range residuals vector.
	for(i = 0; i < num_obs; i++)
	{
		lsq_state->prange_residuals += lsq_mat->zpos[i]*lsq_mat->Wpos[i]*lsq_mat->zpos[i];
	}

	// Initialization of statistics for doppler.
	lsq_state->num_vel_dof = num_obs - lsq_state->num_vel_pars + num_dummy_obs_ibpl;
	lsq_state->doppler_residuals = 0.0;

	// Calculate the weighted squared module of the doppler residuals vector.
	for(i = 0; i < num_obs; i++)
	{
		lsq_state->doppler_residuals += lsq_mat->zvel[i]*lsq_mat->Wvel[i]*lsq_mat->zvel[i];
	}

	// Normalize the module of pseudo-range residuals vector with the number of DoF.
	if (lsq_state->num_pos_dof > 0)
	{
		lsq_state->prange_residuals = sqrt(lsq_state->prange_residuals / lsq_state->num_pos_dof);
	}

	// Normalize the module of doppler residuals vector with the number of DoF.
	if (lsq_state->num_vel_dof > 0)
	{
		lsq_state->doppler_residuals = sqrt(lsq_state->doppler_residuals / lsq_state->num_vel_dof);
	}

#if CALCULATE_IBPL == 1

	// Mark IBPL flag as unavailable.
	lsq_state->ibpl_valid = 0;

	compute_dop(lsq_state->pos, lsq_state->num_pos_pars, lsq_state->DOPmat_pos, &lsq_state->hdop, &lsq_state->vdop);
	compute_dop(lsq_state->pos, lsq_state->num_vel_pars, lsq_state->DOPmat_vel, &lsq_state->hdop_vel, &lsq_state->vdop_vel);

	// Calculate IBPL.
	if(result == 1 && lsq_state->num_pos_dof > 0)
	{
		lsq_state->ibpl_valid = GICSRx_LSQ_IBPL(lsq_state);
	}

#endif

	float varinv_pos = 0;
	float varinv_vel = 0;

	for(i = 0; i < num_obs; i++)
	{
		varinv_pos += lsq_mat->Wpos[i];
		varinv_vel += lsq_mat->Wvel[i];
	}

	if(num_obs >= 4)
	{
		lsq_state->hdop *= sqrt(varinv_pos/num_obs);
		lsq_state->vdop *= sqrt(varinv_pos/num_obs);

		lsq_state->hdop_vel *= sqrt(varinv_vel/num_obs);
		lsq_state->vdop_vel *= sqrt(varinv_vel/num_obs);
	}
	else
	{
		lsq_state->hdop = 1e6;
		lsq_state->vdop = 1e6;
	}
}

/* This function calculates the DOP matrix, considering for that purpose the geometry and the weighting matrices,
 * since DOPmat = inv(H'*W*H). Cholesky method is used to perform the inverse */
static void compute_DOP_matrix(const float HTWH[NORMAL_MATRIX_SIZE], float DOPmat[NORMAL_MATRIX_SIZE], int dim)
{
	// Initialization for choleskyInv: copy HTWH in the DOP output matrix.
	int i;
	for (i = 0; i < NORMAL_MATRIX_SIZE; i++)
	{
		DOPmat[i] = HTWH[i];
	}
	// Calculate the inverse matrix of (H'*W*H) using the Cholesky decomposition.
	choleskyInv(DOPmat, dim);
}

/* This function calculates the DOP matrix as the inverse of (H'*W*H). Previous processing is required to avoid
 * singularities in the geometry matrix due to the possibility of working with a single constellation.
 */
static char least_squares_dop(double state[NUM_PARS], float H[MAX_CHANNELS_G][NUM_PARS], float W[MAX_CHANNELS_G],
		float HTWH[NORMAL_MATRIX_SIZE], float HTW[NUM_PARS][MAX_CHANNELS_G], float DOP_mat[NORMAL_MATRIX_SIZE],
		int num_obs, int num_pars, char DOPmat_io)
{
	// Declaration of variables.
	int i, j, k;

	char result = 0;

	// Check that there are enough measurements to calculate a valid solution.
	if(num_obs >= num_pars)
	{
		// Check whether the DOP matrix, calculated as inv(H'*W*H), has been previously computed.
		// If not, calculate it.
		if(DOPmat_io == 0)
		{
			for (i = 0; i < num_pars; i++)
			{
				for (j = 0; j < num_obs; j++)
				{
					HTW[i][j] = H[j][i] * W[j];
				}
			}

			// Compute normal matrix HTWH (packed form)
			int index = 0;
			for (i = 0; i < num_pars; i++)
			{
				for (j = 0; j <= i; j++)
				{
					HTWH[index] = 0.0;
					for (k = 0; k < num_obs; k++)
					{
						HTWH[index] += HTW[i][k] * H[k][j];
					}
					index++;
				}
			}

			// Decompose HTWH matrix into a lower triangular matrix L.
			if (choleskyDecomp(HTWH, num_pars) == 1)
			{
				// Calculate INV = pinv(HTWH) using L.
				compute_DOP_matrix(HTWH, DOP_mat, num_pars);

				// Result is valid, set flag to 0.
				result = 1;
			}
		}
		else
		{
			// DOP matrix is trusted to have been calculated properly.
			result = 1;
		}
	}
	return result;

} // END of function least_squares_dop

/*
 * This function solves the least-square approach for the estimation problem. Provided with the geometry
 * matrix, weighting matrix, and the innovation / residuals vector, calculates the new value for the state
 * vector. It also calculates the DOP matrix whose statistics can be used by other external functions.
 */
static char GICSRxLeastSquares(float H[MAX_CHANNELS_G][NUM_PARS], float W[MAX_CHANNELS_G], float z[MAX_CHANNELS_G],
		double state[NUM_PARS], float DOP_mat[NORMAL_MATRIX_SIZE], double *eps, int num_obs, int num_pars,
		float HTWH[NORMAL_MATRIX_SIZE], float HTW[NUM_PARS][MAX_CHANNELS_G], char DOPmat_io)
{
	// Declaration of variables.
	int i, j, k;

	char result = 0;

	float dx[NUM_PARS];
	float HTWz[NUM_PARS];

	char invertible = least_squares_dop(state, H, W, HTWH, HTW, DOP_mat, num_obs, num_pars, DOPmat_io);
	if(invertible == 1)
	{
		// Compute right-hand side vector HT*W *z
		for (i = 0; i < num_pars; i++)
		{
			HTWz[i] = 0.0;
			for (k = 0; k < num_obs; k++)
			{
				HTWz[i] += HTW[i][k] * z[k];
			}
		}

		triangularSolve (HTWH, HTWz, num_pars, HTWz);
		triangularSolveT(HTWH, HTWz, num_pars, dx);

		// Add innovation.
		for(i = 0; i < num_pars; i++)
		{
			state[i] += dx[i];
		}

		// Update convergence parameters.
		for(*eps = 0, i = 0; i < num_pars; i++)
		{
			*eps += dx[i]*dx[i];
		}

		// update residuals
		for( i = 0; i < num_obs; i++ ){
			for( j = 0; j < NUM_PARS; j++  ){
				z[i] -= H[i][j] * dx[j];
			}
		}

		// Result is valid, set flag to 1.
		result = 1;
	}
	return result;

} // END of function GICSRxLeastSquares

/* This function calculates the GARLIC least-squares navigation solution and the IBPL integrity solution.
 * Two operations mode are available: calculating each epoch the satellite positions and velocities with the
 * ephemeris functions, which is time-consuming, or using ephemeris propagation, which provides a very accurate
 * approach for the satellite state vector but optimizing the CPU resources.
 */
#if USE_EPHEMERIS_PROP == 1
char GICSRxWNavSol(gnssdata_t *gnssdata, SP3STRUCT_D *pSp3, lsq_state_t *lsq_state, lsq_matrix_t *lsq_mat, kconf_t *kconf, propephem_t *propephem, nav_sol_t *nav_sol, double *init_pos, double *init_vel)
#else
char GICSRxWNavSol(gnssdata_t *gnssdata, SP3STRUCT_D *pSp3, lsq_state_t *lsq_state, lsq_matrix_t *lsq_mat, kconf_t *kconf, nav_sol_t *nav_sol, double *init_pos, double *init_vel)
#endif
{
	int count = 0;
	double eps1 = DELTA_X_LSQ + 1.;
	double eps2 = DELTA_X_LSQ + 1.;

	float HTWH_pos [NORMAL_MATRIX_SIZE];
	float HTWH_vel [NORMAL_MATRIX_SIZE];
	float HTW_pos  [NUM_PARS][MAX_CHANNELS_G];
	float HTW_vel  [NUM_PARS][MAX_CHANNELS_G];

	// Solution flag set to false.
	char result = 0;
	char update_eph = 0;
	int num_dummy_obs_ibpl = 0;

	// In case the receiver has been set in time mark re-synchronization state, the first_fix flag in gnssdata
	// has been set to 0, while first_position may be 1. Set LSQ position flag also to 0 to start again.
	/*if(first_position == 1 && gnssdata->first_fix == 0)
	{
		first_position = 0;
	}*/

	// An estimation of the position is computed.
	do
	{
		// Initialization of indexer.
#if USE_EPHEMERIS_PROP == 1
		reconstruct_obs(gnssdata, pSp3, lsq_state, lsq_mat, (RECEIVER_MODEL == 0 ? result : 1), propephem, kconf, count, &update_eph);
#else
		reconstruct_obs(gnssdata, pSp3, lsq_state, lsq_mat, (RECEIVER_MODEL == 0 ? result : 1), NULL, kconf, count, &update_eph);
#endif

#if CHECK_OUTLIERS == 1
		check_outliers(gnssdata, count, init_pos);
#endif
		int num_dummy_obs = build_obs_matrices(gnssdata, lsq_state, lsq_mat, kconf);
		num_dummy_obs_ibpl = num_dummy_obs;

		// Calculate user position with a weighted least squares method (WLSQ).
		lsq_state->num_pos_pars = NUM_PARS;
		if (lsq_state->total_LOS + num_dummy_obs >= 3 + NUM_SYS + NUM_FREQ - 1)
		{
			result = GICSRxLeastSquares(lsq_mat->G, lsq_mat->Wpos, lsq_mat->zpos, lsq_state->pos, lsq_state->DOPmat_pos, &eps1, lsq_state->total_obs+num_dummy_obs, lsq_state->num_pos_pars, HTWH_pos, HTW_pos, (result && count > 1));
		} else {
			result = 0;
		}

		// Calculate user velocity with a weighted least squares method (WLSQ).
		/* Re-design of the matrix to assign all the measurements (from any
		 * constellation) to the same clock drift.
		 * In case that there are 2 measurements sources (2 receivers combined)
		 * the measurements are split to the corresponding receiver source
		 */
		if(count > 0 && result)
		{
			// Reset number of dummy measurements
			num_dummy_obs = 0;
			int k, i;

			// Set all the measurements for the clock drift estimation (one receiver option)
			if (TWO_RX == 0)
			{
				lsq_state->num_vel_pars = 4;
				for (k=0; k < lsq_state->total_obs; k++)
				{
					lsq_mat->G[k][3] = 1;
				}
			}
			// Split the measurements per receiver source for each clock drift estimation
			else
			{
				lsq_state->num_vel_pars = 5;
				int rx1_flag = 0, rx2_flag = 0, num = 0;
				for(i=0; i < gnssdata->noOfChannelsAv; i++)
				{
					if (gnssdata->OBS[i].prange_status == OBS_VALID)
					{
						// Assumed that GPS and GAL measurements are from the first receiver
						int sys = signal2system(gnssdata->OBS[i].sigFlag);
						if (sys == SYS_GPS || sys == SYS_GAL)
						{
							lsq_mat->G[num][3] = 1;
							lsq_mat->G[num][4] = 0;
							rx1_flag = 1;
						}
						// Assumed that GLO and BEI measurements are from the second receiver
						else if (sys == SYS_GLO || sys == SYS_BEI)
						{
							lsq_mat->G[num][3] = 0;
							lsq_mat->G[num][4] = 1;
							rx2_flag = 1;
						}
						num++;
					}
				}
				// Add dummy measurements in case of missing measurements for any of the receiver sources
				if (rx1_flag == 0)
				{
					memset(lsq_mat->G[num],0,sizeof(lsq_mat->G[num]));
					lsq_mat->G[num][3]		 = 1;
					lsq_mat->Wvel[num] 		 = 1;
					lsq_mat->zvel[num] 		 = 0;
					num_dummy_obs++;
				}
				if (rx2_flag == 0)
				{
					memset(lsq_mat->G[num],0,sizeof(lsq_mat->G[num]));
					lsq_mat->G[num][4]		 = 1;
					lsq_mat->Wvel[num] 		 = 1;
					lsq_mat->zvel[num] 		 = 0;
					num_dummy_obs++;
				}
			}
			result = GICSRxLeastSquares(lsq_mat->G, lsq_mat->Wvel, lsq_mat->zvel, lsq_state->vel, lsq_state->DOPmat_vel, &eps2, lsq_state->total_obs+num_dummy_obs, lsq_state->num_vel_pars, HTWH_vel, HTW_vel, (result && count > 1));
		}

		// Increase iteration counter.
		count++;

	} // END do-while: continue if convergence has not been achieved and maximum iterations have not been reached.
	while((eps1 > DELTA_X_LSQ || eps2 > DELTA_X_LSQ) && count < MAX_ITER_LSQ && result == 1);

	// If the maximum number of iterations has been reached, proceed to mark the solution as wrong.
	if(count == MAX_ITER_LSQ)
	{
		result = 0;
	}

	// Store the number of iterations required until convergence.
	lsq_state->num_iter = count;

	// Set position validity flag.
	lsq_state->pos_valid = result;



	// Compute IBPL.
	compute_lsq_indicators(lsq_state, lsq_mat, lsq_state->total_LOS, result , num_dummy_obs_ibpl);

	// Condition for the first position fix.
	if(first_position == 0 && result == 1)
	{
		first_position = result;
	}

	// Check if a solution is not available for the current epoch.
	if((result &= first_position) == 0)
	{
		// Initialize.
		memcpy(lsq_state->pos, init_pos, sizeof(lsq_state->pos));
		memcpy(lsq_state->vel, init_vel, sizeof(lsq_state->vel));
	}

	// Store the value of the navigation solution in the output structure.
	extract_nav_solution(gnssdata->week, gnssdata->tow, lsq_state, nav_sol, (result == 0 ? NOSOL : WLSQ_SOL));

#if USE_EPHEMERIS_PROP == 1
	// If a valid clock bias estimation is available, and no ephemeris has been already updated in the current epoch, refresh one channel.
	if(result == 1 && update_eph == 0)
	{
		update_ephemeris(gnssdata, lsq_state, propephem, pSp3);
	}
#endif

	// Return the status of the solution calculation.
	return result;

} // END of function GICSRxWNavSol

