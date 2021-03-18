#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "GICSRxKalman.h"

#include "matrix.h"
#include "algebra.h"
#include "GICSRxObs.h"
#include "GICSRxIBPL.h"
#include "garlic_functions.h"

// PR validity window.
static int using_PR[MAXWINOBS] = {1};

// GPS time when last Kalman update was produced.
static double last_update = -1.0;

// Sort in terms of elevation (from minimum to maximum) satellites with doppler available.
t_sort_index array_sort [MAX_CHANNELS_G];

#if USE_WEIGHTING_MODEL == 1
/* Returns the current module of the user velocity, including the innovations performed up to
 * this point in the Kalman Filter through the update function.
 */
static float get_velocity(const exchange_t *exchange)
{
	float tmpvel_vec[3];

	tmpvel_vec[0] =  exchange->ins_state_ecef[3] + exchange->kal_state_ecef[3];
	tmpvel_vec[1] =  exchange->ins_state_ecef[4] + exchange->kal_state_ecef[4];
	tmpvel_vec[2] =  exchange->ins_state_ecef[5] + exchange->kal_state_ecef[5];

	return Vec3Norm(tmpvel_vec);
}
#endif

#if USE_IMU == 1

/* This function stores the KF state vector (in ECEF coordinates) into another vector in
 * navigation frame coordinates (position in geodetics and velocity in NED)
 */
static void get_nav_state_from_ecef(double *ins_state_ecef, double *ins_state_nav)
{
	// Convert position and velocity into the proper system.
	ECEFtoNAV_pos(&ins_state_ecef[0], &ins_state_nav[0]);
	ECEFtoNAV_vel(&ins_state_ecef[3], &ins_state_nav[3], ins_state_nav[0], ins_state_nav[1]);

	// The rest of parameters remain equal.
	int k;
	for(k = 6; k < KNSTATES; k++)
	{
		ins_state_nav[k] = ins_state_ecef[k];
	}
}

/* This function stores the KF state vector (in navigation frame coordinates, i.e. position in geodetics
 * and velocity in NED) into another vector in ECEF coordinates.
 */
static void get_ecef_state_from_nav(double *ins_state_nav, double *ins_state_ecef)
{
	// Convert position and velocity into the proper system.
	NAVtoECEF_pos(&ins_state_ecef[0], &ins_state_nav[0]);
	NAVtoECEF_vel(&ins_state_ecef[3], &ins_state_nav[3], ins_state_nav [0], ins_state_nav[1]);

	// The rest of parameters remain equal.
	int k;
	for(k = 6; k < KNSTATES; k++)
	{
		ins_state_ecef[k] = ins_state_nav[k];
	}
}

/* Update the value of a state variable using a Gauss-Markov model (noise sample is provided) */
static double gauss_markov(double state, double noise_sample, float delta_t, float tau)
{
	double state_out;

	state_out = state*exp(-delta_t/tau) + noise_sample;

	return state_out;
}
#endif

/* Auxiliary function used to apply quicksort. */
static int compare_structs(const void *a, const void *b)
{
	int result = 0;

	t_sort_index *struct_a = (t_sort_index *) a;
	t_sort_index *struct_b = (t_sort_index *) b;

	if(struct_a->CN0 < struct_b->CN0)
	{
		result = +1;
	}
	else if(struct_a->CN0 > struct_b->CN0)
	{
		result = -1;
	}

	return result;

} // END of function compare_structs

/* Sort the GNSS observations list */
void sort_gnss_observations(gicsrx_t *gicsrx, gnssdata_t *gnssdata, char ignore_call, t_sort_index array_copy[MAX_CHANNELS_G])
{
	int k;

	// The call to this function must be ignored if it is called from the hybrid kalman filter
	// and the GNSS-only configuration is also activated, since the sorting has already been done.
	if(ignore_call == 0)
	{
		// Sort measurements in terms of elevation. Also check that CN0 and elevation masks are
		// satisfied.
		for(k = 0; k < gnssdata->noOfChannelsAv; k++)
		{
			if (check_valid_obs(gicsrx, &gnssdata->OBS[k]))
			{
				array_sort[k].CN0  = gnssdata->OBS[k].S1;
			}
			else
			{
				array_sort[k].CN0 = 0;
			}
			array_sort[k].index = k;
		}

		// Apply sorting algorithm to sort the satellites in terms of their elevation.
		// Higher elevations shall be used first since they are supposed to be the better ones.
		qsort(array_sort, gnssdata->noOfChannelsAv, sizeof(t_sort_index), compare_structs);
	}

	if(array_copy != NULL)
	{
		memcpy(array_copy, array_sort, sizeof(array_sort));
	}
}

/* This function updates the Kalman Filter state and covariance matrix, taking as input the innovation / residual and the noise measurement.
 * It can be used both for updating the Hybrid KF (18-states) or the GNSS-only KF (9-states), by indicating properly the size of the state.
 * Another required inputs are, obviously, the KF state and covariance, and the line of sight for which the innovation has been calculated.
 * Notice that the measurement might be rejected if a chi-square test, which compares the magnitude of the innovation with the expected variance,
 * fails.
 */
char update_kalman(unsigned int size, const float *h_vec, float sigma_meas, float residual, float *cov_mat, double *state_vec, float chi2_thr)
{
	int i;

	// Declaration of variables
	float res_var, chi_square;

	// Declaration of matrix variables.
	float PH   [KNSTATES];
	float K_mat[KNSTATES];

	// Calculate PH = P*h
	matSymVecMul_f(PH, cov_mat, h_vec, size);

	// Calculate h*P*h, that determines the variance of the residual
	// for the current measurement.
	DotProduct(h_vec, PH, size, res_var);

	// Add inherent variance of the measurement.
	res_var = res_var + sigma_meas*sigma_meas;

	// Obtain chi-square value.
	chi_square = residual*residual/res_var;

#if USE_CHI2TEST == 1
	// Check whether the square of the modulus of residual vector
	// normalized by the variance of the residual is bounded by the
	// configurable parameter chi_test, in order to pass the test.
	if (chi_square > chi2_thr)
	{
		return 0;
	}
#endif

	// Update state vector.
	float res_var_inv = 1.0/res_var;
	for (i = 0; i < size; i++)
	{
		K_mat[i]       = PH[i] * res_var_inv;
		state_vec[i]  += K_mat[i] * residual;
	}

	// Update the noise state covariance matrix.
	// P+ = (I - K*H)*P-
	update_matrix(size, K_mat, PH, cov_mat);

	// Return test result.
	return 1;

} // END of function update_kalman

/* This function is used to initialize the elements of the state and covariance matrix of the Kalman Filter related
 * to GNSS information, i.e. position and velocity, and clock information, using the navigation solution calculated
 * by the least-squares algorithm. It can be used both for the GNSS-only KF and the Hybrid KF, indicating properly
 * the size of the state vector and covariance matrix. In the case of Hybrid KF, the part of the state and covariance
 * related to attitude and IMU elements (accelerometer and gyro biases) still need to be initialized.
 */
void initialize_gnss_state (unsigned int size, const lsq_state_t *lsq_state, float *kal_cov_mat, double *kal_state)
{
	unsigned int cov_matrix_size = (((size + 1) * size) / 2);

	int m, n;

	for (m = 0; m < cov_matrix_size; m++)
	{
		kal_cov_mat[m] = 0.0;
	}

	// Store covariance values in P matrix.
	// Position and velocity variables covariances
	for (m = 0; m < 3; m++)
	{
		for(n = 0 ; n <= m; n++)
		{
			set_sym_element(kal_cov_mat, m, n, get_sym_element(lsq_state->DOPmat_pos, m, n));
			set_sym_element(kal_cov_mat, m+3, n+3, get_sym_element(lsq_state->DOPmat_vel, m, n));
		}
	}

	// Clocks with position covariances
	for (m = 0; m < 3; m++)
	{
		for (n = 0; n < NUM_SYS; n++)
		{
			if (lsq_state->num_LOS[n] > 0)
			{
				set_sym_element(kal_cov_mat, 6+n, m, get_sym_element(lsq_state->DOPmat_pos, 3+n, m));
			}
		}
		// Drift with velocity covariances
		set_sym_element(kal_cov_mat, 6+NUM_SYS+NUM_FREQ-1, 3+m, get_sym_element(lsq_state->DOPmat_vel, 3, m));
		if (TWO_RX == 1) { set_sym_element(kal_cov_mat, 6+NUM_SYS+NUM_FREQ-1+TWO_RX, 3+m, get_sym_element(lsq_state->DOPmat_vel, 3+TWO_RX, m)); }

	}

	// Clock and isb covariances
	double ref_clk_cov = get_sym_element(lsq_state->DOPmat_pos, 3, 3);
	set_sym_element(kal_cov_mat, 6, 6, ref_clk_cov);

	for (m = 0; m < NUM_SYS-1; m++)
	{
		if (lsq_state->num_LOS[m+1] > 0)
		{
			double isb_cov = get_sym_element(lsq_state->DOPmat_pos, 4+m, 4+m);
			set_sym_element(kal_cov_mat, 7+m, 6, sqrt(ref_clk_cov*ref_clk_cov + isb_cov*isb_cov));
			set_sym_element(kal_cov_mat, 7+m, 7+m, get_sym_element(lsq_state->DOPmat_pos, 4+m, 4+m));
		}
		else
		{
			set_sym_element(kal_cov_mat, 7+m, 7+m, 10000);
		}
	}

	// Drift covariance
	set_sym_element(kal_cov_mat, 6+NUM_SYS+NUM_FREQ-1, 6+NUM_SYS+NUM_FREQ-1, get_sym_element(lsq_state->DOPmat_vel, 3, 3));
	if (TWO_RX == 1) { set_sym_element(kal_cov_mat, 6+NUM_SYS+NUM_FREQ-1+TWO_RX, 6+NUM_SYS+NUM_FREQ-1+TWO_RX, get_sym_element(lsq_state->DOPmat_vel, 3+TWO_RX, 3+TWO_RX)); }

	// Inter-freq covariances
	for (m = 0; m < NUM_FREQ-1; m++)
	{
		if (lsq_state->used_freq[m+1] > 0)
		{
			for (n = 0; n < 3; n++)
			{
				set_sym_element(kal_cov_mat, 6+NUM_SYS+m, n, get_sym_element(lsq_state->DOPmat_pos, 3+NUM_SYS+m, n));
			}
			set_sym_element(kal_cov_mat, 6+NUM_SYS+m, 6+NUM_SYS+m, get_sym_element(lsq_state->DOPmat_pos, 3+NUM_SYS+m, 3+NUM_SYS+m));
		}
		else
		{
			set_sym_element(kal_cov_mat, 6+NUM_SYS+m, 6+NUM_SYS+m, 10000);
		}
	}

	// Store initial state vector with the results of the Least-Squares algorithm.
	// Position
	kal_state[0] = lsq_state->pos[0];
	kal_state[1] = lsq_state->pos[1];
	kal_state[2] = lsq_state->pos[2];
	// Velocity
	kal_state[3] = lsq_state->vel[0];
	kal_state[4] = lsq_state->vel[1];
	kal_state[5] = lsq_state->vel[2];
	// Clock
	kal_state[6] = lsq_state->pos[3];
	// ISB: init to 0 if there is no measurements for that constellation
	for (m = 1; m < NUM_SYS; m++)
	{
		if (lsq_state->num_LOS[m] > 0)
		{
			if (TWO_RX == 1)
			{
				if      (m == SYS_GLO) { kal_state[6+m] = lsq_state->pos[3+m]; }
				else if (m == SYS_BEI) { kal_state[6+m] = lsq_state->pos[3+m] - kal_state[6+SYS_GLO];}
				else                   { kal_state[6+m] = lsq_state->pos[3+m] - lsq_state->pos[3];}
			}
			else
			{
				kal_state[6+m] = lsq_state->pos[3+m] - lsq_state->pos[3];
			}
		}
		else
		{
			kal_state[6+m] = 0.0;
		}
	}
	// Inter-frequency bias: init to 0 if there is no measurements for that frequency
	for (m = 0; m < NUM_FREQ-1; m++)
	{
		if (lsq_state->used_freq[m+1] > 0)
		{
			kal_state[6+NUM_SYS+m] = lsq_state->pos[3+NUM_SYS+m];
		}
		else
		{
			kal_state[6+NUM_SYS+m] = 0.0;
		}
	}
	// Drift
	kal_state[6+NUM_SYS+NUM_FREQ-1] = lsq_state->vel[3];
	if (TWO_RX == 1) { kal_state[6+NUM_SYS+NUM_FREQ-1+TWO_RX] = lsq_state->vel[3+TWO_RX]; }

} // END of function initialize_gnss_state

#if CALC_GNSS_ONLY == 1
static void copy_gnss_state(exchange_t *exchange)
{
	int m, n;
	for (m = 0; m < COV_MATRIX_SIZE; m++) { exchange->kal_cov_mat[m] = 0.0; }

	// Store covariance values in P matrix.
	for (m = 0; m < GNSS_KSTATES; m++)
	{
		// Store initial state vector with the results of the GNSS-only KF algorithm.
		exchange->ins_state_ecef[m] = exchange->gnss_state[m];

		for(n = 0 ; n <= m; n++)
		{
			set_sym_element(exchange->kal_cov_mat, m, n, get_sym_element(exchange->gnss_cov_mat, m, n));
		}
	}
}
#endif

/* This function is used by the Hybrid Kalman Filter to initialize the part of the state and covariance matrix
 * related to attitude and IMU elements (accelerometer and gyro biases).
 */
#if USE_IMU == 1
static void initialize_nav_state (gicsrx_t *gicsrx)
{
	double *gnss_state = gicsrx->exchange.ins_state_ecef;

	float Cbn[3][3];

	ECEFtoNAV_pos(&gnss_state[0], &gicsrx->exchange.ins_state_nav[0]);
	ECEFtoNAV_vel(&gnss_state[3], &gicsrx->exchange.ins_state_nav[3], gicsrx->exchange.ins_state_nav[0], gicsrx->exchange.ins_state_nav[1]);

	int m, n;
	for(m = 6; m < 9+KAC; m++)
	{
		gicsrx->exchange.ins_state_nav[m] = gnss_state[m];
	}

	gicsrx->exchange.ins_state_nav[ 9+KAC] = atan2(gicsrx->exchange.ins_state_nav[4], gicsrx->exchange.ins_state_nav[3]);
	gicsrx->exchange.ins_state_nav[10+KAC] = 0.00;
	gicsrx->exchange.ins_state_nav[11+KAC] = 0.00;

	memset(Cbn, 0, sizeof(Cbn));

	Cbn[0][0] = +cos(gicsrx->exchange.ins_state_nav[9+KAC]);
	Cbn[0][1] = -sin(gicsrx->exchange.ins_state_nav[9+KAC]);
	Cbn[1][0] = +sin(gicsrx->exchange.ins_state_nav[9+KAC]);
	Cbn[1][1] = +cos(gicsrx->exchange.ins_state_nav[9+KAC]);
	Cbn[2][2] = 1;

	mat2quaternion(Cbn, gicsrx->exchange.quaternion);

	set_sym_element(gicsrx->exchange.kal_cov_mat,  9+KAC,  9+KAC, 1.0);
	set_sym_element(gicsrx->exchange.kal_cov_mat, 10+KAC, 10+KAC, 1.0);
	set_sym_element(gicsrx->exchange.kal_cov_mat, 11+KAC, 11+KAC, 1.0);

	float cov_0 = 1e-2;
	for(m = 12+KAC; m < KNSTATES; m++)
	{
		gicsrx->exchange.ins_state_nav[m] = 0.00;
		for(n = 12+KAC; n < KNSTATES; n++)
		{
			float cov_value = (m == n) ? cov_0 : 0.00;
			set_sym_element(gicsrx->exchange.kal_cov_mat, m, n, cov_value);
		}
	}

} // END of function initialize_nav_state
#endif

/* Initializes the state and covariance matrix of the Hybrid Kalman Filter */
static void initialize (int week, double tow, gicsrx_t *gicsrx, lsq_state_t *lsq_state)
{
	// Initialize state vector and covariance state matrix.
	gicsrx->exchange.week      = week;
	gicsrx->exchange.tow       = tow;
	gicsrx->exchange.prop_time = 0.0;

	// Initialize position, velocity and clock bias and drift.
#if CALC_GNSS_ONLY == 0
	initialize_gnss_state(KNSTATES, lsq_state, gicsrx->exchange.kal_cov_mat, gicsrx->exchange.ins_state_ecef);
#else
	copy_gnss_state(&gicsrx->exchange);
#endif

	// Initialize attitude and accelerometer and gyro biases.
#if USE_IMU == 1
	initialize_nav_state (gicsrx);
#endif

	// Initialize Kalman-based IBPL elements.
#if CALCULATE_IBPL == 1
	initialize_kibpl_indicators(&gicsrx->kconf, lsq_state, gicsrx->exchange.tow, &gicsrx->exchange.kibpl);
#endif
	// Store last update epoch and exit from reset state.
	gicsrx->exchange.kalman_on = 1;

	last_update = gicsrx->exchange.tow;
	gicsrx->exchange.timeOfLastReset = gicsrx->exchange.tow;

	// Debug message.
//	fprintf(stderr, "\nRESET A at %.3f\n", gicsrx->exchange.tow);

} // END of function initialize

/* Extract state information from the navigation data structures */
static void extract_nav_solution(exchange_t *exchange, nav_sol_t *nav_sol, quality_t quality)
{
	int k;
	for(k = 0; k < KNSTATES; k++)
	{
		if(isnan(exchange->ins_state_ecef[k]))
		{
			exchange->kalman_on = 0;
			return;
		}
	}

	char distance_condition = 1;
#if USE_IMU == 1
	distance_condition = (quality == KALMAN_HIMU_0 || (exchange->kal_prev_sats != 0 && exchange->kal_used_sats != 0));
#if CALC_GNSS_ONLY == 1
	distance_condition |= (exchange->kalman2_on == 0);
#endif
#endif

	if(quality >= nav_sol->quality && distance_condition == 1)
	{
		nav_sol->week = exchange->week;
		nav_sol->tow  = exchange->tow;

		nav_sol->kalman_valid = exchange->kalman_on;
		nav_sol->position[0]  = exchange->ins_state_ecef[0];
		nav_sol->position[1]  = exchange->ins_state_ecef[1];
		nav_sol->position[2]  = exchange->ins_state_ecef[2];
		nav_sol->velocity[0]  = exchange->ins_state_ecef[3];
		nav_sol->velocity[1]  = exchange->ins_state_ecef[4];
		nav_sol->velocity[2]  = exchange->ins_state_ecef[5];
		nav_sol->clock_bias   = exchange->ins_state_ecef[6];
		for (k = 0; k < NUM_SYS-1; k++) { nav_sol->isb[k] = exchange->ins_state_ecef[6+k+1]; }
		for (k = 0; k < NUM_FREQ-1; k++){ nav_sol->ifreqBias[k] = exchange->ins_state_ecef[6+NUM_SYS+k]; }
		nav_sol->clock_drift[0]  = exchange->ins_state_ecef[6+NUM_SYS+NUM_FREQ-1];
		if (TWO_RX == 1) {nav_sol->clock_drift[1]  = exchange->ins_state_ecef[6+NUM_SYS+NUM_FREQ-1+TWO_RX];}

		if(quality == KALMAN_ONLY_0 || quality == KALMAN_HIMU_0)
		{
			nav_sol->curr_used_sats = nav_sol->prev_used_sats = 0;
			for (k = 0; k < NUM_SYS; k++) { nav_sol->num_LOS[k] = 0; }
		}
		else
		{
			nav_sol->prev_used_sats = exchange->kal_prev_sats;
			nav_sol->curr_used_sats = exchange->kal_used_sats;
			for (k = 0; k < NUM_SYS; k++) { nav_sol->num_LOS[k] = exchange->num_LOS[k]; }
		}

		memcpy(nav_sol->prn_used, exchange->kal_prn_used, sizeof(exchange->kal_prn_used));

		nav_sol->prop_time = exchange->prop_time;
		nav_sol->quality   = quality;
	}

#if CALC_GNSS_ONLY == 1
	const int gnss_sats_thr = 4;
	if(quality == KALMAN_HIMU_0 && exchange->gnss_used_sats < gnss_sats_thr)
	{
		nav_sol->kalman_valid = exchange->kalman_on;
		nav_sol->position[0]  = exchange->ins_state_ecef[0];
		nav_sol->position[1]  = exchange->ins_state_ecef[1];
		nav_sol->position[2]  = exchange->ins_state_ecef[2];
		nav_sol->clock_bias   = exchange->ins_state_ecef[6];
		for (k = 0; k < NUM_SYS-1; k++) { nav_sol->isb[k] = exchange->ins_state_ecef[7+k]; }
		nav_sol->quality      = quality;
	}
#endif
}

#if (CALCULATE_IBPL == 1 || CALC_GNSS_ONLY == 1 || USE_IMU == 0)
static void prop_cov_mat_deterministic(unsigned int size, float *cov_mat, float delta_t)
{
	sparse_matrix_t F_mat;
	reset_sparse_matrix(size, &F_mat);

	// State noise covariance propagation
	// P(k) = M*P(k-1)*M';
	add_element_to_sparse_matrix(0, 3, delta_t, &F_mat);
	add_element_to_sparse_matrix(1, 4, delta_t, &F_mat);
	add_element_to_sparse_matrix(2, 5, delta_t, &F_mat);

	// Clock drift term.
	add_element_to_sparse_matrix(6, 6+NUM_SYS+NUM_FREQ-1, delta_t, &F_mat);
	if (TWO_RX == 1) { add_element_to_sparse_matrix(6+SYS_GLO, 6+NUM_SYS+NUM_FREQ-1+TWO_RX, delta_t, &F_mat); }

	update_matrix_sparse(size, &F_mat, cov_mat);
}

static void prop_cov_mat_stochastic(float *cov_mat, float delta_t, float sigma2_v, float sigma2_clk, float *sigma2_isb, float *sigma2_ifreqb)
{
	// Get P(k) = M*P(k-1)*M' + Q(k)
	float delta_t2 = delta_t  * delta_t;
	float delta_t3 = delta_t2 * delta_t;

	int i;

	// Position and velocity
	for(i = 0; i < 3; i++)
	{
		update_sym_element(cov_mat, i,   i,   sigma2_v * delta_t3 / 3.0);
		update_sym_element(cov_mat, i+3, i,   sigma2_v * delta_t2 / 2.0);
		update_sym_element(cov_mat, i+3, i+3, sigma2_v * delta_t);
	}
	// Clock
	update_sym_element(cov_mat, 6, 6, sigma2_clk  * delta_t3 / 3.0);
	if (TWO_RX == 1) {update_sym_element(cov_mat, 6+SYS_GLO, 6+SYS_GLO, sigma2_clk  * delta_t3 / 3.0);}
	// ISB
	for(i = 0; i < NUM_SYS-1; i++)
	{
		update_sym_element(cov_mat, 7+i, 7+i, sigma2_isb[i] * delta_t);
	}
	// Inter-frequency bias
	for(i = 0; i < NUM_FREQ-1; i++)
	{
		update_sym_element(cov_mat, 6+NUM_SYS+i, 6+NUM_SYS+i, sigma2_ifreqb[i] * delta_t);
	}
	// Drift
	update_sym_element(cov_mat, 6+NUM_SYS+NUM_FREQ-1, 6, sigma2_clk  * delta_t2 / 2.0);
	update_sym_element(cov_mat, 6+NUM_SYS+NUM_FREQ-1, 6+NUM_SYS+NUM_FREQ-1, sigma2_clk  * delta_t);
	if (TWO_RX == 1)
	{
		update_sym_element(cov_mat, 6+NUM_SYS+NUM_FREQ-1+TWO_RX, 6+SYS_GLO, sigma2_clk  * delta_t2 / 2.0);
		update_sym_element(cov_mat, 6+NUM_SYS+NUM_FREQ-1+TWO_RX, 6+NUM_SYS+NUM_FREQ-1+TWO_RX, sigma2_clk  * delta_t);
	}
}

/* This function can be used to propagate the covariance matrix of the GNSS-only KF. It is declared as a global function
 * since it is invoked externally when the Kalman-based IBPL is computed.
 */
void prop_cov_mat(unsigned int size, float *cov_mat, float delta_t, kconf_t *kconf)
{
	prop_cov_mat_deterministic(size, cov_mat, delta_t);
	prop_cov_mat_stochastic(cov_mat, delta_t, kconf->sigma2_vel, kconf->sigma2_clk, kconf->sigma2_isb, kconf->sigma2_ifreqb);
}
#endif

#if USE_IMU == 0
static void propagate_gnss(gicsrx_t *gicsrx, gnssdata_t *gnssdata)
{
	// Calculate propagation time in the Kalman filter.
	gicsrx->exchange.prop_time = gnssdata->tow - gicsrx->exchange.tow;

	// Store timestamp.
	gicsrx->exchange.week = gnssdata->week;
	gicsrx->exchange.tow  = gnssdata->tow;

	double delta_t    = gicsrx->exchange.prop_time;
	double *pIState_0 = gicsrx->exchange.ins_state_ecef;

	// State vector propagation. State(k) = M*State(k-1)
	//
	// INS absolute state.
	pIState_0[0] = pIState_0[0] + delta_t*pIState_0[3];
	pIState_0[1] = pIState_0[1] + delta_t*pIState_0[4];
	pIState_0[2] = pIState_0[2] + delta_t*pIState_0[5];

#if CLK_STEERING == 0
	pIState_0[6] = pIState_0[6] + delta_t*pIState_0[6+NUM_SYS+NUM_FREQ-1];
	if (TWO_RX == 1) {pIState_0[6+SYS_GLO] = pIState_0[6+SYS_GLO] + delta_t*pIState_0[6+NUM_SYS+NUM_FREQ-1+TWO_RX];}
#endif

	memset(gicsrx->exchange.kal_state_ecef, 0, KNSTATES*sizeof(double));

	// State noise covariance propagation
	// P(k) = M*P(k-1)*M' + Q(k);
	prop_cov_mat(KNSTATES, gicsrx->exchange.kal_cov_mat, delta_t, &gicsrx->kconf);

} // END of function propagate
#endif

#if USE_IMU == 1
/* This function is invoked internally by the function update_filter_state() to update the attitude matrix (the quaternion vector)
 * using the innovations calculated by the update function in the Kalman Filter. Two different methods for the attitude rotation
 * are defined: the first one corresponds to the implementation performed in ATLANTIDA Simulink model, while the second is
 * considered to be a more accurate approach, implemented specifically for the GARLIC software.
 */
static void update_attitude_matrix(float drot_x, float drot_y, float drot_z, double latitude, double longitude, float quaternion[4], double *att_angles)
{
	float Cbn [3][3], Cbe[3][3], Cne[3][3], Cen[3][3], Cbn_u[3][3], Cbe_u[3][3];
	float Rmat[3][3];

#if UPDATE_ATT_MODE == 0
	// Update method used in ATLANTIDA.
	Rmat[0][0] = +cos(drot_z)*cos(drot_y);
	Rmat[0][1] = +sin(drot_z)*cos(drot_x) + cos(drot_z)*sin(drot_y)*sin(drot_x);
	Rmat[0][2] = +sin(drot_z)*sin(drot_x) - cos(drot_z)*sin(drot_y)*cos(drot_x);
	Rmat[1][0] = -sin(drot_z)*cos(drot_y);
	Rmat[1][1] = +cos(drot_z)*cos(drot_x) - sin(drot_z)*sin(drot_y)*sin(drot_x);
	Rmat[1][2] = +cos(drot_z)*sin(drot_x) + sin(drot_z)*sin(drot_y)*cos(drot_x);
	Rmat[2][0] = +sin(drot_y);
	Rmat[2][1] = -cos(drot_y)*sin(drot_x);
	Rmat[2][2] = +cos(drot_y)*cos(drot_x);

#else
	// Update method used in GARLIC.
	float normrot = sqrt(drot_x*drot_x + drot_y*drot_y + drot_z*drot_z);

	// Check whether the rotation is high enough to avoid numerical divergences.
	if(normrot < 1e-5)
	{
		Rmat[0][0] = 1.00;
		Rmat[0][1] = +drot_z;
		Rmat[0][2] = -drot_y;
		Rmat[1][0] = -drot_z;
		Rmat[1][1] = 1.00;
		Rmat[1][2] = +drot_x;
		Rmat[2][0] = +drot_y;
		Rmat[2][1] = -drot_x;
		Rmat[2][2] = 1.00;
	}
	else
	{
		// The rotation matrix calculated in this section applies to the attitude
		// matrix expressed as the required change of coordinates to convert from
		// body-frame to ECEF coordinate frame. It is calculated as:
		//
		// Rmat = I + sin(rot)/rot*R + (1-cos(rot))/(rot^2)*R^2
		//
		// where rot is the module of the rotation vector, and R the skew matrix
		// associated to the rotation vector (the attitude innovation in the KF).

		float Rmat1[3][3], Rmat2[3][3];

		float sinrot = sin(normrot);
		float cosrot = cos(normrot);

		Rmat1[0][0] = 0.00;
		Rmat1[0][1] = +drot_z/normrot;
		Rmat1[0][2] = -drot_y/normrot;
		Rmat1[1][0] = -drot_z/normrot;
		Rmat1[1][1] = 0.00;
		Rmat1[1][2] = +drot_x/normrot;
		Rmat1[2][0] = +drot_y/normrot;
		Rmat1[2][1] = -drot_x/normrot;
		Rmat1[2][2] = 0.00;

		matMul_f(Rmat2, Rmat1, Rmat1, 3, 3, 3);

		Rmat[0][0] = 1 + sinrot*Rmat1[0][0] + (1-cosrot)*Rmat2[0][0];
		Rmat[0][1] = 0 + sinrot*Rmat1[0][1] + (1-cosrot)*Rmat2[0][1];
		Rmat[0][2] = 0 + sinrot*Rmat1[0][2] + (1-cosrot)*Rmat2[0][2];
		Rmat[1][0] = 0 + sinrot*Rmat1[1][0] + (1-cosrot)*Rmat2[1][0];
		Rmat[1][1] = 1 + sinrot*Rmat1[1][1] + (1-cosrot)*Rmat2[1][1];
		Rmat[1][2] = 0 + sinrot*Rmat1[1][2] + (1-cosrot)*Rmat2[1][2];
		Rmat[2][0] = 0 + sinrot*Rmat1[2][0] + (1-cosrot)*Rmat2[2][0];
		Rmat[2][1] = 0 + sinrot*Rmat1[2][1] + (1-cosrot)*Rmat2[2][1];
		Rmat[2][2] = 1 + sinrot*Rmat1[2][2] + (1-cosrot)*Rmat2[2][2];
	}
#endif
	/////////////////////////////////////////////////////////////////////////////////////////////////

	// Obtain rotation matrix from ECEF to navigation frame.
	NAVtoECEF_mat(Cne, latitude, longitude);
	Transpose(Cen,Cne,3,3);

	// Obtain rotation matrix from body to navigation frame.
	quaternion2mat(quaternion, Cbn);

	// Calculate the rotation matrix from body to ECEF and update it with the update of the rotation
	// matrix calculated in the previous step.
	matMul_f(Cbe,   Cne,  Cbn,   3,3,3);
	matMul_f(Cbe_u, Rmat, Cbe,   3,3,3);
	matMul_f(Cbn_u, Cen,  Cbe_u, 3,3,3);

	// Store the new attitude matrix (body to navigation).
	mat2quaternion(Cbn_u, quaternion);

	// Calculate attitude angles from quaternions.
	get_attitude_angles_from_matrix(Cbn_u, &att_angles[0], &att_angles[1], &att_angles[2]);
}
#endif

/* This function takes a doppler measurement and uses it to update the Kalman Filter */
static char update_prange(exchange_t *exchange, obsdata_t *obs, double tow, float chi2_thr)
{
	// Auxiliar line of sight vector (9 or 18 elements).
	float h_vec[KNSTATES];

	// Compute the innovation for the input measurement and the current state vector.
	compute_prange_residual(exchange->ins_state_ecef, exchange->kal_state_ecef, obs, tow, exchange->tow, &obs->prange_residual, h_vec, KNSTATES);

	// Use the computed innovation to update the state and covariance matrix in KF.
	char passed = update_kalman(KNSTATES, h_vec, obs->prange_sigma, obs->prange_residual, exchange->kal_cov_mat, exchange->kal_state_ecef, chi2_thr);

	// Flag pseudo-range observation as rejected.
	if (passed == 1)
	{
		// FIXME: find better way to count LOS when is not from L1 band
		if(signal2freq(obs->sigFlag) == FREQ_L1E1) {
			exchange->num_LOS[signal2system(obs->sigFlag)]++;
		}
	}
	else
	{
		obs->prange_status = REJ_CHISQUARE_TEST;
	}
	return passed;

} // END of function update_prange

/* This function takes a doppler measurement and uses it to update the Kalman Filter */
static char update_doppler(exchange_t *exchange, obsdata_t * obs, float chi2_thr)
{
	// Auxiliar line of sight vector (9 or 18 elements).
	float h_vec[KNSTATES];

	// Compute the innovation for the input measurement and the current state vector.
	compute_doppler_residual(exchange->ins_state_ecef, exchange->kal_state_ecef, obs, &obs->doppler_residual, h_vec, KNSTATES);

	// Check whether the residual is not above the threshold.
	if (fabs(obs->doppler_residual) > 100.0)
	{
		obs->doppler_status = REJ_CHISQUARE_TEST;
		return 0;
	}

	// Use the computed innovation to update the state and covariance matrix in KF.
	char passed = update_kalman(KNSTATES, h_vec, obs->doppler_sigma, obs->doppler_residual, exchange->kal_cov_mat, exchange->kal_state_ecef, chi2_thr);

	// Flag doppler observation as rejected.
	if (passed == 0)
	{
		obs->doppler_status = REJ_CHISQUARE_TEST;
	}
	return passed;

} // END of function update_Doppler

#if USE_IMU == 1

static void compute_ceb_matrix(float quaternion[4], double latitude, double longitude, float Ceb[3][3])
{
	float Cbn[3][3], Cnb[3][3], Cen[3][3];

	quaternion2mat(quaternion, Cbn);
	ECEFtoNAV_mat(Cen, latitude, longitude);

	Transpose(Cnb, Cbn, 3,3);
	matMul_f (Ceb, Cnb, Cen, 3,3,3);
}

/* This function implements the non-holonomic restriction in the Hybrid Kalman Filter.
 * The non-holonomic concept is related with the assumption that a ground vehicle does not skid. This
 * implies that the velocity vector of the vehicle is directly related with the attitude, which gives us
 * constraints between the state variables and reduce the number of degrees of freedom of the navigation
 * system. This assumption is introduced in the Kalman Filter as an artificial observation, including
 * a noise measurement, since this assumption is not absolutely true (ground vehicles do skid).
 */
static char update_non_holonomic(exchange_t *exchange)
{
	float h_vec[KNSTATES], user_vel_e[3];
	float Ceb[3][3];

	float residual_1, residual_2;

	// Temporary Kalman state and covariance matrix
	float kal_cov_mat     [COV_MATRIX_SIZE];
	double kal_state_ecef [KNSTATES];

	// Different noise measurements for each constraint.
	float sigma_nh_hor = 0.05;
	float sigma_nh_ver = 0.30;
	float chi2_thr = 1e6;
	float residuals_threshold = 0.5;

	// Copy current Kalman state in temporary variables
	int k;
	for (k = 0; k < COV_MATRIX_SIZE; k++)
	{
		kal_cov_mat[k] = exchange->kal_cov_mat[k];
	}

	for (k = 0; k < KNSTATES; k++)
	{
		kal_state_ecef[k] = exchange->kal_state_ecef[k];
	}

	// Compute transformation matrix from body axis to ECEF
	compute_ceb_matrix( exchange->quaternion, exchange->ins_state_nav[0], exchange->ins_state_nav[1], Ceb);

	// Initialize partials vector
	memset(h_vec, 0, sizeof(h_vec));

	// Obtain first user velocity (in ECEF).
	user_vel_e[0] = exchange->ins_state_ecef[3] + exchange->kal_state_ecef[3];
	user_vel_e[1] = exchange->ins_state_ecef[4] + exchange->kal_state_ecef[4];
	user_vel_e[2] = exchange->ins_state_ecef[5] + exchange->kal_state_ecef[5];

	// Build observation vector with the respective partial derivatives.
	h_vec[3]  = Ceb[1][0];
	h_vec[4]  = Ceb[1][1];
	h_vec[5]  = Ceb[1][2];

	h_vec[ 9+KAC] = - user_vel_e[2] * Ceb[1][1] + user_vel_e[1] * Ceb[1][2];
	h_vec[10+KAC] = + user_vel_e[2] * Ceb[1][0] - user_vel_e[0] * Ceb[1][2];
	h_vec[11+KAC] = - user_vel_e[1] * Ceb[1][0] + user_vel_e[0] * Ceb[1][1];

	// Calculate residual in the horizontal plane for the non-holonomic restriction.
	residual_1 = -(Ceb[1][0]*user_vel_e[0] + Ceb[1][1]*user_vel_e[1] + Ceb[1][2]*user_vel_e[2]);

	// Add constraint if residual does not exceed threshold
	if (fabs(residual_1) < residuals_threshold)
	{
		update_kalman(KNSTATES, h_vec, sigma_nh_hor, residual_1, kal_cov_mat, kal_state_ecef, chi2_thr);
	}

	// Obtain second user velocity (in ECEF).
	user_vel_e[0] = exchange->ins_state_ecef[3] + exchange->kal_state_ecef[3];
	user_vel_e[1] = exchange->ins_state_ecef[4] + exchange->kal_state_ecef[4];
	user_vel_e[2] = exchange->ins_state_ecef[5] + exchange->kal_state_ecef[5];

	// Build observation vector with the respective partial derivatives.
	h_vec[3]  = Ceb[2][0];
	h_vec[4]  = Ceb[2][1];
	h_vec[5]  = Ceb[2][2];

	h_vec[ 9+KAC] = - user_vel_e[2] * Ceb[2][1] + user_vel_e[1] * Ceb[2][2];
	h_vec[10+KAC] = + user_vel_e[2] * Ceb[2][0] - user_vel_e[0] * Ceb[2][2];
	h_vec[11+KAC] = - user_vel_e[1] * Ceb[2][0] + user_vel_e[0] * Ceb[2][1];

	// Calculate residual in the vertical plane for the non-holonomic restriction.
	residual_2 = -(Ceb[2][0]*user_vel_e[0] + Ceb[2][1]*user_vel_e[1] + Ceb[2][2]*user_vel_e[2]);

	// Add constraint if residual does not exceed threshold
	if (fabs(residual_2) < residuals_threshold)
	{
		update_kalman(KNSTATES, h_vec, sigma_nh_ver, residual_2, kal_cov_mat, kal_state_ecef, chi2_thr);
	}

	// Copy updated state into final state
	if (fabs(kal_state_ecef[9+KAC]) < 0.005 && fabs(kal_state_ecef[10+KAC]) < 0.005 && fabs(kal_state_ecef[11+KAC]) < 0.005)
	{
		for (k = 0; k < COV_MATRIX_SIZE; k++)
		{
			exchange->kal_cov_mat[k] = kal_cov_mat[k];
		}

		for (k = 0; k < KNSTATES; k++)
		{
			exchange->kal_state_ecef[k] = kal_state_ecef[k];
		}
	}

	return 0;
}
#endif

#if USE_IMU == 0 && USE_CONSTANT_HEIGHT == 1
/* This function adds a restriction to the filter considering that the height of the vehicle is constant. */
static void use_constant_height(exchange_t *exchange)
{
	float h_vec[KNSTATES], Cen[3][3], ecefvel[3], residual;
	float sigma_height = 1.00, chi2_thr = 1e6;

	double ecefpos[3], geodpos[3];

	// Obtain user position (in ECEF).
	ecefpos[0] = exchange->ins_state_ecef[0] + exchange->kal_state_ecef[0];
	ecefpos[1] = exchange->ins_state_ecef[1] + exchange->kal_state_ecef[1];
	ecefpos[2] = exchange->ins_state_ecef[2] + exchange->kal_state_ecef[2];

	// Calculate user position (in geodetics).
	ECEFtoNAV_pos(ecefpos, geodpos);

	// Calculate rotation matrix from ECEF to NAV coordinate frame.
	ECEFtoNAV_mat(Cen, geodpos[0], geodpos[1]);

	// Initialize partials vector
	memset(h_vec, 0, sizeof(h_vec));

	// Obtain user velocity (in ECEF).
	ecefvel[0] = exchange->ins_state_ecef[3] + exchange->kal_state_ecef[3];
	ecefvel[1] = exchange->ins_state_ecef[4] + exchange->kal_state_ecef[4];
	ecefvel[2] = exchange->ins_state_ecef[5] + exchange->kal_state_ecef[5];

	// Set observation vector.
	h_vec[3]  = Cen[2][0];
	h_vec[4]  = Cen[2][1];
	h_vec[5]  = Cen[2][2];

	// Calculate residual for the down-component of the velocity.
	residual = -(Cen[2][0]*ecefvel[0] + Cen[2][1]*ecefvel[1] + Cen[2][2]*ecefvel[2]);

	// Update KF status.
	update_kalman(KNSTATES, h_vec, sigma_height, residual, exchange->kal_cov_mat, exchange->kal_state_ecef, chi2_thr);
}
#endif

/* This function is used to update the state vector in the Kalman Filter (Hybrid or GNSS-only).
 * Since the update function is invoked sequentially, and the innovations are accumulated in an
 * independent vector, when the update process has finished the new information must be moved to
 * the state vector. These innovations are added linearly in the case of position and velocity state
 * variables. Nevertheless, in the case of attitude and sensor biases, these innovations are applied
 * differently: the attitude innovations correspond to a rotation in the attitude matrix, and the
 * bias innovations, noise that shall be added using a gauss-markov process.
 */
#if USE_IMU == 1
static void update_filter_state(exchange_t *exchange, float tau_acc, float tau_gyr)
#else
static void update_filter_state(exchange_t *exchange)
#endif
{
	// Update state vector (rx, ry, rz, vx, vy, vz, clk, [isb1, isb2, ...], dri).
	int k;
	for(k = 0; k < KNSTATES; k++)
	{
		exchange->ins_state_ecef[k]  += exchange->kal_state_ecef[k];
	}

#if USE_IMU == 1

	// Update attitude (heading, roll, pitch).
	float drot_x = exchange->kal_state_ecef[ 9+KAC];
	float drot_y = exchange->kal_state_ecef[10+KAC];
	float drot_z = exchange->kal_state_ecef[11+KAC];

	update_attitude_matrix(drot_x, drot_y, drot_z, exchange->ins_state_nav[0], exchange->ins_state_nav[1], exchange->quaternion, &exchange->ins_state_ecef[9+KAC]);

	// Update accelerometer bias using Gauss-Markov model. Take into account that what we obtain in the
	// innovation vector for the accelerometer and gyro biases have opposite sign.
	exchange->ins_state_ecef[12+KAC] = gauss_markov(exchange->ins_state_ecef[12+KAC], -exchange->kal_state_ecef[12+KAC], exchange->prop_time, tau_acc);
	exchange->ins_state_ecef[13+KAC] = gauss_markov(exchange->ins_state_ecef[13+KAC], -exchange->kal_state_ecef[13+KAC], exchange->prop_time, tau_acc);
	exchange->ins_state_ecef[14+KAC] = gauss_markov(exchange->ins_state_ecef[14+KAC], -exchange->kal_state_ecef[14+KAC], exchange->prop_time, tau_acc);

	// Update gyro bias using Gauss-Markov model.
	exchange->ins_state_ecef[15+KAC] = gauss_markov(exchange->ins_state_ecef[15+KAC], -exchange->kal_state_ecef[15+KAC], exchange->prop_time, tau_gyr);
	exchange->ins_state_ecef[16+KAC] = gauss_markov(exchange->ins_state_ecef[16+KAC], -exchange->kal_state_ecef[16+KAC], exchange->prop_time, tau_gyr);
	exchange->ins_state_ecef[17+KAC] = gauss_markov(exchange->ins_state_ecef[17+KAC], -exchange->kal_state_ecef[17+KAC], exchange->prop_time, tau_gyr);

	get_nav_state_from_ecef(exchange->ins_state_ecef, exchange->ins_state_nav);
#endif

	memset(exchange->kal_state_ecef, 0, sizeof(exchange->kal_state_ecef));
}

/* This function takes as input the GNSS information (satellite states, pseudo-ranges and doppler measurements)
 * and uses them to update the Hybrid Kalman Filter state. This function is very similar to the function used
 * to update the GNSS-only Kalman Filter, and it has been separated for the sake of clarity in this .c file. Besides,
 * it shares several functionalities, as the capability of sorting the PRN mask in terms of CN0 or elevation, or
 * flagging measurements as unavailable. This is taken into account in order not to perform twice the same operations.
 * In particular, sorting of measurements is not performed when the GNSS-only KF is turned on (flag CALC_GNSS_ONLY).
 */
static void update_gnss(gicsrx_t *gicsrx, gnssdata_t *gnssdata)
{
	// Declaration of algorithmic variables.
	int j, k, count_PR = 0, used_PR = 0;
	float acc_PR;

	// Sort in terms of elevation (from maximum to minimum) satellites with doppler available.
	sort_gnss_observations(gicsrx, gnssdata, CALC_GNSS_ONLY, NULL);

	// Apply non-holonomic restriction (velocity and attitude are directly related, i.e. the vehicle does not skid).
	// If GNSS-only KF is configured, apply constant height restriction if it is defined.
#if USE_IMU == 1
	update_non_holonomic(&gicsrx->exchange);
#elif USE_CONSTANT_HEIGHT == 1
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
				update_doppler(&gicsrx->exchange, &gnssdata->OBS[j], gicsrx->kconf.chi_test1);
			}
		}
	}

#if USE_WEIGHTING_MODEL == 1
	float normvel = get_velocity(&gicsrx->exchange);
#endif

	char filterConverged = (gnssdata->tow - gicsrx->exchange.timeOfLastReset) > gicsrx->kconf.ktimeNoSmooth;

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
				if(update_prange(&gicsrx->exchange, &gnssdata->OBS[j], gnssdata->tow, gicsrx->kconf.chi_test1))
				{
					// Increase the number of used PR.
					gicsrx->exchange.kal_prn_used[used_PR++] = gnssdata->OBS[j].PRN;
				}
			}
			if(gnssdata->OBS[j].prange_status != OBS_OK_BUT_UNUSED)
			{
				// Increase the number of available PR.
				count_PR++;
			}
		}
	}

#if CALCULATE_IBPL == 1
	gicsrx->exchange.kibpl.common.ibpl_epoch_flag = 1;
#if USE_IMU == 1
	float acc_indicator_factor = gicsrx->kconf.acc_indicator_factor_imu;
#else
	float acc_indicator_factor = gicsrx->kconf.acc_indicator_factor_gnss;
#endif
	compute_kibpl_indicators(
			gicsrx->exchange.tow,
			KNSTATES,
			&gicsrx->kconf,
			acc_indicator_factor,
			gicsrx->exchange.ins_state_ecef,
			gicsrx->exchange.kal_state_ecef,
			gicsrx->exchange.kal_cov_mat,
			&gicsrx->exchange.kibpl.common,
			gnssdata);
#endif

#if USE_IMU == 1
	// Update INS state vector.
	update_filter_state(&gicsrx->exchange, gicsrx->kconf.tau_acc, gicsrx->kconf.tau_gyr);
#else
	// Update GNSS state vector.
	update_filter_state(&gicsrx->exchange);
#endif

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
		using_PR[0] = (acc_PR > gicsrx->kconf.per_PR);
	}

	if (count_PR > 0)
	{
		last_update = gicsrx->exchange.tow;
	}
	gicsrx->exchange.kal_used_sats = used_PR;

} // END of function update

static void update_residuals(exchange_t *exchange, gnssdata_t *gnssdata)
{
	float h_vec[KNSTATES];

	int k;
	for(k = 0; k < gnssdata->noOfChannelsAv; k++)
	{
		// Compute the residual for the input measurement and the current state vector.
		compute_doppler_residual(exchange->ins_state_ecef, exchange->kal_state_ecef, &gnssdata->OBS[k], &gnssdata->OBS[k].doppler_residual, h_vec, KNSTATES);

		// Compute the residual for the input measurement and the current state vector.
		compute_prange_residual(exchange->ins_state_ecef, exchange->kal_state_ecef, &gnssdata->OBS[k], exchange->tow, exchange->tow, &gnssdata->OBS[k].prange_residual, h_vec, KNSTATES);
	}
}

/* This function implements the initialization and GNSS update of the Hybrid Kalman Filter. Prior to invoke this function
 * it is assumed that the programmer has propagated the state vector and covariance matrix through a mechanization process,
 * using the information from accelerometers and gyros. This capability is already implemented in GICSRxMechanization.c, in
 * the function kalman_hybrid.
 * In case that the programmer has set the flag USE_IMU to 0, this function implements the initialization, propagation and
 * update functions of a GNSS-only Kalman Filter.
 */
void kalman_gnss(gicsrx_t *gicsrx, gnssdata_t *gnssdata, char result, lsq_state_t *lsq_state, nav_sol_t *nav_sol)
{
#if CALC_INSTALL_MAT == 1
	if(gicsrx->kconf.Imat_ready == 0)
	{
		return;
	}
#endif
	// Temporary gicsrx_t data structure used to optimize IBPL computation.
	gicsrx_t gicsrx_prov;
	memset(&gicsrx_prov, 0, sizeof(gicsrx_t));

#if CALC_GNSS_ONLY == 0
	// The frequency counter is used to guarantee that update in Kalman Filter is only performed each
	// N seconds, where N is equal to the value of seconds determined in GNSS_rate.
	gicsrx->exchange.gnss_frequency_counter++;
#endif

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
			// Reset Kalman filter.
			gicsrx->exchange.kalman_on = 0;

			// If KF reset is due to more than maxKalx epochs without measurements, indicate that
			// something is wrong with the measurements/timestamp in the receiver.
			if(last_update > 0 && ((gicsrx->exchange.tow - last_update) > gicsrx->kconf.maxKalx))
			{
				gnssdata->first_fix = 0;
			}
		}

		// Check whether the filter is initialized or not.
		if(gicsrx->exchange.kalman_on == 0)
		{
			// Initialize observation window for PR rejection.
			for(k = 0; k < gicsrx->kconf.win_PR; k++) { using_PR[k] = 1; }

			// If it has been possible to obtain a LS solution, re-initialize to these values.
			// Include some additional conditions to initialize the filter.
			if(result && lsq_state->total_LOS >= gicsrx->kconf.noSats_min && lsq_state->hdop < 2.5)
			{
#if CALC_GNSS_ONLY == 1
				if(nav_sol->quality == KALMAN_ONLY_1 && Vec3Norm(lsq_state->vel) > gicsrx->kconf.hkvel0)
				{
#elif USE_IMU == 1
				if(Vec3Norm(lsq_state->vel) > gicsrx->kconf.hkvel0)
				{
#endif
				// Initialize state vector and covariance state matrix.
				initialize(gnssdata->week, gnssdata->tow, gicsrx, lsq_state);
				return;
#if USE_IMU == 1
			}
#endif
			}
		}
#if USE_IMU == 0
		else
		{
			char filterConverged = (gnssdata->tow - gicsrx->exchange.timeOfLastReset) > gicsrx->kconf.ktimeNoSmooth;

			// Propagate for current epoch.
			propagate_gnss(gicsrx, gnssdata);
			extract_nav_solution(&gicsrx->exchange, nav_sol, filterConverged ? KALMAN_ONLY_0 : NOSOL);
		}
#endif

#if USE_IMU == 1
		// Prepare Kalman filter. Convert from navigation vector to ECEF.
		get_ecef_state_from_nav(gicsrx->exchange.ins_state_nav, gicsrx->exchange.ins_state_ecef);

		if(gicsrx->exchange.kalman_on == 1)
		{
			// Update status of the navigation solution.
			extract_nav_solution(&gicsrx->exchange, nav_sol, KALMAN_HIMU_0);
		}
#endif
		// Update the number of satellites used in the previous epoch.
		gicsrx->exchange.kal_prev_sats = gicsrx->exchange.kal_used_sats;
		gicsrx->exchange.kal_used_sats = 0;

		// Initialize variables to 0.
		for(k = 0; k < NUM_SYS; k++) { gicsrx->exchange.num_LOS[k] = 0; }

		if (gnssdata->noOfChannelsAv > 0)
		{
			if (gicsrx->exchange.kalman_on == 1)
			{
#if RECEIVER_MODEL == 0
				// Correct clock bias taking into account the step value in clock drift.
				propagate_clock_bias(gicsrx->exchange.ins_state_ecef, gnssdata, lsq_state, CLK_MAXITER);
#endif

#if USE_IMU == 1
				// Update clock bias also in the propagation state vector.
				gicsrx->exchange.ins_state_nav[6] = gicsrx->exchange.ins_state_ecef[6];
#endif

				// Copy current state of gicsrx data structure to avoid wrong propagations in IBPL computation.
				gicsrx_prov = *gicsrx;

				// Update error vector with GNSS measurements.
				update_gnss(&gicsrx_prov, gnssdata);

#if CALCULATE_IBPL == 1
				kibpl_code_t kibpl_return_code =
						GICSRx_KF_IBPL(1, &gicsrx_prov.kconf, gicsrx->exchange.tow, gicsrx_prov.exchange.kal_used_sats, &gicsrx_prov.exchange.kibpl, nav_sol);

				if (kibpl_return_code == KIBPL_OK)
				{
					gicsrx->exchange = gicsrx_prov.exchange;
					if ((gicsrx->exchange.tow - gicsrx->exchange.kibpl.common.last_ibpl_tow) > gicsrx->kconf.maxKalx*2)
					{
						gicsrx->exchange.kalman_on = 0;
					}
					gicsrx->exchange.kibpl.common.last_ibpl_tow = gicsrx->exchange.tow;
				}

				if (kibpl_return_code == KIBPL_RESET_FILTER)
				{
					gicsrx->exchange.kalman_on = 0;
				}

				nav_sol->kalman_ibpl_valid = (kibpl_return_code == KIBPL_OK);
#else
				gicsrx->exchange = gicsrx_prov.exchange;
#endif
				update_residuals(&gicsrx_prov.exchange, gnssdata);
			}
			// Check that at least one pseudo-range measurement has been successfully updated the filter state.
			if(gicsrx_prov.exchange.kal_used_sats > 0)
			{
				// Update status of the navigation solution.
				extract_nav_solution(&gicsrx_prov.exchange, nav_sol, (USE_IMU == 0 ? KALMAN_ONLY_1 : KALMAN_HIMU_1));
			}
#if USE_IMU == 0
			if(nav_sol->quality >= KALMAN_ONLY_0)
			{
				nav_sol->prev_used_sats = gicsrx->exchange.kal_prev_sats;
				nav_sol->curr_used_sats = gicsrx->exchange.kal_used_sats;
			}
#endif
		}
		// Set GNSS update counter to 0.
		gicsrx->exchange.gnss_frequency_counter = 0;
	}
	else
	{
		if (gicsrx->exchange.kalman_on == 1)
		{
#if USE_IMU == 1
			// Prepare Kalman filter. Convert from navigation vector to ECEF.
			get_ecef_state_from_nav(gicsrx->exchange.ins_state_nav, gicsrx->exchange.ins_state_ecef);
#else
			// Propagate for current epoch.
			propagate_gnss(gicsrx, gnssdata);
#endif
			if(gnssdata->noOfChannelsAv > 0)
			{
#if RECEIVER_MODEL == 0
				// Correct clock bias taking into account the step value in clock drift.
				propagate_clock_bias(gicsrx->exchange.ins_state_ecef, gnssdata, lsq_state, CLK_MAXITER);
#endif
#if USE_IMU == 1
				// Update clock bias also in the propagation state vector.
				gicsrx->exchange.ins_state_nav [6] = gicsrx->exchange.ins_state_ecef[6];
#endif
			}
			// Update status of the navigation solution.
			extract_nav_solution(&gicsrx->exchange, nav_sol, (USE_IMU == 0 ? KALMAN_ONLY_0 : KALMAN_HIMU_0));
		}
	}

#if USE_IMU == 0
	// Compute traveled distance. Only if hybrid mode (IMU) is disabled
	calculate_distance_metrics(gicsrx, nav_sol, &gicsrx->exchange.travelled_distance);
#endif

} // END of function kalman_gnss

