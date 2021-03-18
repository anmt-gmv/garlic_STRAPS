#ifndef GICSRXKALMAN_H_
#define GICSRXKALMAN_H_

#include "dbtypedef.h"
#include "GICSRx_defines.h"
#include "Summary_Recorder.h"
#include "Indicators_Recorder.h"

#define MAXWINOBS 15		// Maximum size of observation window.

/* Data structure used to apply quicksort algorithm in terms of CN0 */
typedef struct
{
	double CN0;
	int    index;

} t_sort_index;

/* Sort the GNSS observations list */
void sort_gnss_observations(gicsrx_t *gicsrx, gnssdata_t *gnssdata, char ignore_call, t_sort_index array_sort[MAX_CHANNELS_G]);

/* This function is used to initialize the elements of the state and covariance matrix of the Kalman Filter related
 * to GNSS information, i.e. position and velocity, and clock information, using the navigation solution calculated
 * by the least-squares algorithm. It can be used both for the GNSS-only KF and the Hybrid KF, indicating properly
 * the size of the state vector and covariance matrix. In the case of Hybrid KF, the part of the state and covariance
 * related to IMU elements (attitude, biases) still need to be initialized.
 */
void initialize_gnss_state (unsigned int size, const lsq_state_t *lsq_state, float *kal_cov_mat, double *kal_state);

/* This function updates the Kalman Filter state and covariance matrix, taking as input the innovation / residual and the noise measurement.
 * It can be used both for updating the Hybrid KF (18-states) or the GNSS-only KF (9-states), by indicating properly the size of the state.
 * Another required inputs are, obviously, the KF state and covariance, and the line of sight for which the innovation has been calculated.
 * Notice that the measurement might be rejected if a chi-square test, which compares the magnitude of the innovation with the expected variance,
 * fails.
 */
char update_kalman(unsigned int size, const float *h_vec, float sigma_meas, float residual, float *cov_mat, double *state_vec, float chi2_thr);

#if (CALCULATE_IBPL == 1 || CALC_GNSS_ONLY == 1 || USE_IMU == 0)
/* This function can be used to propagate the covariance matrix of the GNSS-only KF. It is declared as a global function
 * since it is invoked externally when the Kalman-based IBPL is computed.
 */
void prop_cov_mat(unsigned int size, float *cov_mat, float delta_t, kconf_t *kconf);
#endif

/* This function implements the initialization and GNSS update of the Hybrid Kalman Filter. Prior to invoke this function
 * it is assumed that the programmer has propagated the state vector and covariance matrix through a mechanization process,
 * using the information from accelerometers and gyros. This capability is already implemented in GICSRxMechanization.c, in
 * the function kalman_hybrid.
 * In case that the programmer has set the flag USE_IMU to 0, this function implements the initialization, propagation and
 * update functions of a GNSS-only Kalman Filter.
 */
void kalman_gnss(gicsrx_t *gicsrx, gnssdata_t *gnssdata, char result, lsq_state_t *lsq_state, nav_sol_t *nav_sol);

#endif
