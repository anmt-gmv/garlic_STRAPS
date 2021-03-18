#include <stdlib.h>
#include <string.h>

#include "GICSRxKalman.h"

#include "matrix.h"
#include "algebra.h"
#include "GICSRxObs.h"
#include "GICSRxIBPL.h"

#if CALC_GNSS_ONLY == 1

// PR validity window.
static int using_PR[MAXWINOBS] = {1};

// GPS time when last Kalman update was produced.
static double last_update = -1.0;

// GPS time at which last reset was produced.
static double last_reset = -1.0;

/* Returns the current module of the user velocity, including the innovations performed up to
 * this point in the Kalman Filter through the update function.
 */
static float get_velocity(const exchange_t *exchange)
{
	float tmpvel_vec[3];

	tmpvel_vec[0] =  exchange->gnss_state[3] + exchange->gnss_delta[3];
	tmpvel_vec[1] =  exchange->gnss_state[4] + exchange->gnss_delta[4];
	tmpvel_vec[2] =  exchange->gnss_state[5] + exchange->gnss_delta[5];

	return Vec3Norm(tmpvel_vec);
}

/* Initialize the GNSS-only Kalman Filter, using the Least-Squares information */
static void initialize (int week, double tow, gicsrx_t *gicsrx, const lsq_state_t *lsq_state)
{
	static int firstInit = 1;

	// Initialize state vector and covariance state matrix.
	gicsrx->exchange.week      = week;
	gicsrx->exchange.tow       = tow;
	gicsrx->exchange.prop_time = 0.0;

	// GNSS initialization.
	initialize_gnss_state(GNSS_KSTATES, lsq_state, gicsrx->exchange.gnss_cov_mat, gicsrx->exchange.gnss_state);

	// Initialize Kalman-based IBPL elements.
#if CALCULATE_IBPL == 1
	initialize_kibpl_indicators(&gicsrx->kconf, lsq_state, gicsrx->exchange.tow, &gicsrx->exchange.kibpl_gnss_only);
#endif

	// Store last update epoch and exit from reset state.
	last_update = last_reset = gicsrx->exchange.tow;
	gicsrx->exchange.kalman2_on = 1;

	if(gicsrx->exchange.kalman_on == 1)
	{
		firstInit = 0;
	}
	if(firstInit == 1)
	{
		gicsrx->exchange.timeOfLastReset = gicsrx->exchange.tow;
	}

	// Debug message.
	fprintf(stderr, "RESET B at %.3f\n", gicsrx->exchange.tow);

} // END of function initialize

/* Extract state information from the navigation data structures */
static void extract_nav_solution(exchange_t *exchange, nav_sol_t *nav_sol, quality_t quality)
{
	int k;
	for(k = 0; k < GNSS_KSTATES; k++)
	{
		if(isnan(exchange->gnss_state[k]))
		{
			exchange->kalman2_on = 0;
			return;
		}
	}

	if(quality >= nav_sol->quality)
	{
		nav_sol->week 			   = exchange->week;
		nav_sol->tow			   = exchange->tow;

		nav_sol->kalman_valid      = exchange->kalman2_on;
		nav_sol->position[0]       = exchange->gnss_state[0];
		nav_sol->position[1]       = exchange->gnss_state[1];
		nav_sol->position[2]       = exchange->gnss_state[2];
		nav_sol->velocity[0]       = exchange->gnss_state[3];
		nav_sol->velocity[1]       = exchange->gnss_state[4];
		nav_sol->velocity[2]       = exchange->gnss_state[5];
		nav_sol->clock_bias        = exchange->gnss_state[6];
		nav_sol->isb_g1       = exchange->gnss_state[7];
#if KAC > 0
		nav_sol->isb_g2       = exchange->gnss_state[8];
#endif
		nav_sol->clock_drift  = exchange->gnss_state[8+KAC];

		if(quality == KALMAN_ONLY_0)
		{
			nav_sol->curr_used_sats = nav_sol->prev_used_sats = 0;
			nav_sol->used_gps_sats  = nav_sol->used_glo_sats  = nav_sol->used_gal_sats  = 0;
		}
		else
		{
			nav_sol->prev_used_sats    = exchange->gnss_prev_sats;
			nav_sol->curr_used_sats    = exchange->gnss_used_sats;

			nav_sol->used_gps_sats  = exchange->gnss_gps_sats;
			nav_sol->used_glo_sats  = exchange->gnss_glo_sats;
			nav_sol->used_gal_sats  = exchange->gnss_gal_sats;
		}
		nav_sol->prop_time         = exchange->prop_time;
		nav_sol->quality           = quality;
	}
}

/* This function propagates the GNSS-only Kalman Filter state linearly and its covariance matrix */
static void propagate_gnss(gicsrx_t *gicsrx, gnssdata_t *gnssdata)
{
	double delta_t    = gicsrx->exchange.prop_time;
	double *pIState_0 = gicsrx->exchange.gnss_state;

	// State vector propagation. State(k) = M*State(k-1)
	//
	// INS absolute state.
	pIState_0[0] = pIState_0[0] + delta_t*pIState_0[3];
	pIState_0[1] = pIState_0[1] + delta_t*pIState_0[4];
	pIState_0[2] = pIState_0[2] + delta_t*pIState_0[5];

#if CLK_STEERING == 0
	pIState_0[6] = pIState_0[6] + delta_t*pIState_0[8+KAC];
#endif

	memset(gicsrx->exchange.gnss_delta, 0, GNSS_KSTATES*sizeof(double));

	// State noise covariance propagation
	// P(k) = M*P(k-1)*M' + Q(k);
	prop_cov_mat(GNSS_KSTATES, gicsrx->exchange.gnss_cov_mat, delta_t, &gicsrx->kconf);

} // END of function propagate_gnss

/* This function takes a pseudo-range measurement and uses it to update the Kalman Filter */
static char update_prange(exchange_t *exchange, obsdata_t *obs, double tow, float chi2_thr)
{
	// Ancillary line of sight vector (9 elements).
	float h_vec[GNSS_KSTATES];

	// Compute the innovation for the input measurement and the current state vector.
	compute_prange_residual (exchange->gnss_state, exchange->gnss_delta, obs, tow, exchange->tow, &obs->prange_residual, h_vec, GNSS_KSTATES);

	// Use the computed innovation to update the state and covariance matrix in KF.
	char passed = update_kalman(GNSS_KSTATES, h_vec, obs->prange_sigma, obs->prange_residual, exchange->gnss_cov_mat, exchange->gnss_delta, chi2_thr);

	if(passed == 1)
	{
		switch(obs->sysFlag)
		{
		case SYS_GPS:
			exchange->gnss_gps_sats++; break;
		case SYS_GLO:
			exchange->gnss_glo_sats++; break;
		case SYS_GALOS:
		case SYS_GALPRS:
			exchange->gnss_gal_sats++; break;
		default:;
		}
	}

	return passed;

} // END of function update_prange

/* This function takes a doppler measurement and uses it to update the Kalman Filter */
static char update_doppler(exchange_t *exchange, obsdata_t *obs, float chi2_thr)
{
	// Ancillary line of sight vector (9 elements).
	float h_vec[GNSS_KSTATES];

	// Compute the innovation for the input measurement and the current state vector.
	compute_doppler_residual(exchange->gnss_state, exchange->gnss_delta, obs, &obs->doppler_residual, h_vec, GNSS_KSTATES);

	// Use the computed innovation to update the state and covariance matrix in KF.
	char passed = update_kalman(GNSS_KSTATES, h_vec, obs->doppler_sigma, obs->doppler_residual, exchange->gnss_cov_mat, exchange->gnss_delta, chi2_thr);

	return passed;

} // END of function update_Doppler

#if USE_CONSTANT_HEIGHT == 1
/* This function adds a restriction to the filter considering that the height of the vehicle is constant. */
static void use_constant_height(exchange_t *exchange)
{
	float h_vec[GNSS_KSTATES], Cen[3][3], ecefvel[3], residual;
	float sigma_height = 1.00, chi2_thr = 1e6;

	double ecefpos[3], geodpos[3];

	// Obtain user position (in ECEF).
	ecefpos[0] = exchange->gnss_state[0] + exchange->gnss_delta[0];
	ecefpos[1] = exchange->gnss_state[1] + exchange->gnss_delta[1];
	ecefpos[2] = exchange->gnss_state[2] + exchange->gnss_delta[2];

	// Calculate user position (in geodetics).
	ECEFtoNAV_pos(ecefpos, geodpos);

	// Calculate rotation matrix from ECEF to NAV coordinate frame.
	ECEFtoNAV_mat(Cen, geodpos[0], geodpos[1]);

	// Initialize partials vector
	memset(h_vec, 0, sizeof(h_vec));

	// Obtain user velocity (in ECEF).
	ecefvel[0] = exchange->gnss_state[3] + exchange->gnss_delta[3];
	ecefvel[1] = exchange->gnss_state[4] + exchange->gnss_delta[4];
	ecefvel[2] = exchange->gnss_state[5] + exchange->gnss_delta[5];

	// Set observation vector.
	h_vec[3]  = Cen[2][0];
	h_vec[4]  = Cen[2][1];
	h_vec[5]  = Cen[2][2];

	// Calculate residual for the down-component of the velocity.
	residual = -(Cen[2][0]*ecefvel[0] + Cen[2][1]*ecefvel[1] + Cen[2][2]*ecefvel[2]);

	// Update KF status.
	update_kalman(GNSS_KSTATES, h_vec, sigma_height, residual, exchange->gnss_cov_mat, exchange->gnss_delta, chi2_thr);

}
#endif

/* This function takes as input the GNSS information (satellite states, pseudo-ranges and doppler measurements)
 * and uses them to update the GNSS-only Kalman Filter state. This function is very similar to the function used
 * to update the Hybrid Kalman Filter, and it has been separated for the sake of clarity in this .c file. Besides,
 * it shares several functionalities, as the capability of sorting the PRN mask in terms of CN0 or elevation, or
 * flagging measurements as unavailable. This is taken into account in order not to perform twice the same operations.
 */
static void update_gnss(gicsrx_t *gicsrx, gnssdata_t *gnssdata)
{
	// Declaration of algorithmic variables.
	int j, k, count_PR = 0, used_PR = 0, count_PR_glo = 0;
	float acc_PR, normvel;

	t_sort_index array_sort[MAX_CHANNELS_G];

	// Sort in terms of elevation (from maximum to minimum) satellites with doppler available.
	sort_gnss_observations(gicsrx, gnssdata, 0, array_sort);

#if USE_CONSTANT_HEIGHT == 1
	use_constant_height(&gicsrx->exchange);
#endif
	// Perform Doppler observation updates. Doppler is used in the first place since it is supposed
	// to be a more reliable measure.
	for(k = 0; k < gnssdata->noOfChannelsAv; k++)
	{
		// Check whether it is a valid satellite (elevation higher than 0 degrees).
		if(array_sort[k].CN0 != 0)
		{
			j = array_sort[k].index;  // Get index in OBS structure.

			if (gnssdata->OBS[j].doppler_status == OBS_VALID)
			{
#if USE_WEIGHTING_MODEL == 1
				doppler_weighting_model(&gnssdata->OBS[j]);
#else
				gnssdata->OBS[j].doppler_sigma = gicsrx->kconf.sigma_doppler;
#endif
				update_doppler(&gicsrx->exchange, &gnssdata->OBS[j], gicsrx->kconf.chi_test2);
			}
		}
	}

	normvel = get_velocity(&gicsrx->exchange);

	char filterConverged = (gnssdata->tow - last_reset) > gicsrx->kconf.ktimeNoSmooth;

	// Apply PR measurement updates.
	for(k = 0; k < gnssdata->noOfChannelsAv; k++)
	{
		if(array_sort[k].CN0 != 0)
		{
			j = array_sort[k].index;	// Get index in OBS structure.

			if (gnssdata->OBS[j].prange_status == OBS_VALID)
			{
#if USE_WEIGHTING_MODEL == 1
				prange_weighting_model(&gnssdata->OBS[j], normvel, gicsrx->kconf.kvelStopThr, filterConverged ? gicsrx->kconf.kvelStopFactor : 1);
#else
				gnssdata->OBS[j].prange_sigma = gicsrx->kconf.sigma_prange;
#endif
				if(update_prange(&gicsrx->exchange, &gnssdata->OBS[j], gnssdata->tow, gicsrx->kconf.chi_test2))
				{
					used_PR++;			// Increase the number of used PR.
				}
				if(gnssdata->OBS[j].sysFlag == SYS_GLO)
				{
					count_PR_glo++;
				}
			}
			count_PR++;				// Increase the number of available PR.
		}
	}

#if CALCULATE_IBPL == 1
	gicsrx->exchange.kibpl_gnss_only.common.ibpl_epoch_flag = 1;
	compute_kibpl_indicators(
			gicsrx->exchange.tow,
			GNSS_KSTATES,
			&gicsrx->kconf,
			gicsrx->kconf.acc_indicator_factor_gnss,
			gicsrx->exchange.gnss_state,
			gicsrx->exchange.gnss_delta,
			gicsrx->exchange.gnss_cov_mat,
			&gicsrx->exchange.kibpl_gnss_only.common,
			gnssdata);
#endif
	// Update state vector.
	for(k = 0; k < GNSS_KSTATES; k++)
	{
		gicsrx->exchange.gnss_state[k] += gicsrx->exchange.gnss_delta[k];
	}

	// Checking PR window rejection algorithm applied for reset criterion.
	if(count_PR > 0)
	{
		// Used PR percentage.
		acc_PR = ((float)(used_PR)/count_PR)*100;

		// Shift values in PR-observation window.
		for (k = gicsrx->kconf.win_PR-1; k > 0; k--)
		{
			using_PR[k] = using_PR[k-1];
		}
		// Check if PR percentage is above the defined threshold per_PR.
		using_PR[0] = (acc_PR > gicsrx->kconf.per_PR) && 
					 !(count_PR_glo >= 3 && gicsrx->exchange.gnss_glo_sats == 0);
	}

	if (count_PR > 0)
	{
		last_update = gicsrx->exchange.tow;
	}
	gicsrx->exchange.gnss_used_sats = used_PR;

} // END of function update

/* This function implements a GNSS-only Kalman Filter, which shall be invoked to every epoch in order to perform properly
 * the propagation and update (if new measurements for the current time are provided) operations in the navigation algorithm.
 */
void kalman_gnss_only(gicsrx_t *gicsrx, gnssdata_t *gnssdata, char result, lsq_state_t *lsq_state, nav_sol_t *nav_sol)
{
	// The frequency counter is used to guarantee that update in Kalman Filter is only performed each
	// N seconds, where N is equal to the value of seconds determined in GNSS_rate.
	gicsrx->exchange.gnss_frequency_counter++;

	// Check that in this epoch update of KF shall be calculated.
	if(gicsrx->exchange.gnss_frequency_counter >= gicsrx->kconf.GNSS_rate)
	{
		// Declaration of variables.
		int k, prsum;

		// Check the number of rejected epochs in window.
		// using_PR positions store 1 if at least certain percentage of PR was used (determined
		// by configuration), or 0 if it was not, during an observation window of length win_PR.
		for(k = 0, prsum = 0; k < gicsrx->kconf.win_PR; k++) { prsum += using_PR[k]; }

		// Reset the filter if there have been rejected more epochs than allowed
		// or elapsed time is higher than threshold.
		if(prsum < gicsrx->kconf.tol_PR || (last_update > 0 && ((gicsrx->exchange.tow - last_update) > gicsrx->kconf.maxKalx)))
		{
			// RESET Kalman filter.
			gicsrx->exchange.kalman2_on = 0;
		}

		// Check whether the filter is initialized or not.
		if(gicsrx->exchange.kalman2_on == 0)
		{
			// Initialize observation window for PR rejection.
			for(k = 0; k < gicsrx->kconf.win_PR; k++) { using_PR[k] = 1; }

			// If it has been possible to obtain a LS solution, re-initialize to these values.
			// Include some additional conditions to initialize the filter.
			if(result == 1 && lsq_state->num_used >= gicsrx->kconf.noSats_min && lsq_state->hdop < 2.5)
			{
				initialize(gnssdata->week, gnssdata->tow, gicsrx, lsq_state);
				return;
			}
		}
		else
		{
			// Propagate for current epoch.
			propagate_gnss(gicsrx, gnssdata);

			// Update status of the navigation solution.
			extract_nav_solution(&gicsrx->exchange, nav_sol, KALMAN_ONLY_0);
		}

		// Update the number of satellites used in the previous epoch.
		gicsrx->exchange.gnss_prev_sats = gicsrx->exchange.gnss_used_sats;
		gicsrx->exchange.gnss_used_sats = 0;

		// Set to zero the number of satellites.
		gicsrx->exchange.gnss_gps_sats  = 0;
		gicsrx->exchange.gnss_glo_sats  = 0;
		gicsrx->exchange.gnss_gal_sats  = 0;

		if (gnssdata->noOfChannelsAv > 0)
		{
			if (gicsrx->exchange.kalman2_on == 1)
			{
#if RECEIVER_MODEL == 0
				// Correct clock bias taking into account the step value in clock drift.
				propagate_clock_bias(gicsrx->exchange.gnss_state, gnssdata, lsq_state, CLK_MAXITER);
#endif
				// Update error vector with GNSS measurements.
				update_gnss(gicsrx, gnssdata);

#if CALCULATE_IBPL == 1
				kibpl_code_t kibpl_return_code =
						GICSRx_KF_IBPL(0, &gicsrx->kconf, gicsrx->exchange.tow, gicsrx->exchange.gnss_used_sats, &gicsrx->exchange.kibpl_gnss_only, nav_sol);

				if (kibpl_return_code == KIBPL_OK)
				{
					gicsrx->exchange.kibpl_gnss_only.common.last_ibpl_tow = gicsrx->exchange.tow;
				}
				nav_sol->kalman_ibpl_valid = (kibpl_return_code == KIBPL_OK);
#endif
			}
			// Check that at least one pseudo-range measurement has been successfully updated the filter state.
			if(gicsrx->exchange.gnss_used_sats > 0)
			{
				// Update status of the navigation solution.
				extract_nav_solution(&gicsrx->exchange, nav_sol, KALMAN_ONLY_1);
			}
		}
	}
	else
	{
		if (gicsrx->exchange.kalman2_on == 1)
		{
			// Propagate GNSS-only KF, if it is initialized.
			propagate_gnss(gicsrx, gnssdata);

			if(gnssdata->noOfChannelsAv > 0)
			{
#if RECEIVER_MODEL == 0
				// Correct clock bias taking into account the step value in clock drift.
				propagate_clock_bias(gicsrx->exchange.gnss_state, gnssdata, lsq_state, CLK_MAXITER);
#endif
			}
			// Update status of the navigation solution.
			extract_nav_solution(&gicsrx->exchange, nav_sol, KALMAN_ONLY_0);
		}
	}

} // END of function GICSRxKalman
#endif
