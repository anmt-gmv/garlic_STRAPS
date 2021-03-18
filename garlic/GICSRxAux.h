#ifndef GICSRXAUX_H_
#define GICSRXAUX_H_

#include "GICSRx_defines.h"

#if USE_IMU == 1

#if USE_ANTISPOOFING == 1
/* This function implements an anti-spoofing algorithm, which checks the consistency between the IMU
 * acceleration (properly corrected with the corresponding accelerometer biases, attitude and gravity, and
 * expressed in navigation-frame), and GNSS-based acceleration, calculated as the difference between
 * consecutive velocity estimations in a time lapse of 1-second, since the current implementation of the
 * Kalman Filter does not include acceleration in the state vector.
 * The IMU acceleration can be obtained, for instance, as the averaged value over a second of 100 samples,
 * assuming that its value remain constant throughout this period.
 */
double antispoofing_algorithm(gicsrx_t *gicsrx, double acc_imu[3], double acc_kal[3], nav_sol_t *nav_sol);
#endif

#if CALC_GNSS_ONLY == 1
/* This function calculates, when both Hybrid and GNSS-only Kalman Filters are active, the discrepancy
 * between both solutions, in terms of position. When the distance between both is higher than a
 * pre-determined threshold, an alarm flag is returned, so the operator can take specific decisions.
 */
char check_filter_inconsistencies(double hybrid_state[KNSTATES], double gnss_state[GNSS_KSTATES], char noOfSats_GNSS);
#endif

/* This function calculates the initial distance, using for this purpose a buffer that stores the
 * raw imu accelerations (averaged 100 samples in one second) during a certain period (defined by
 * the length of the distance window). In REQ-GARLIC6-PRF-040, a 60 second period previous to the
 * first fix shall be integrated to obtain the initial distance.
 */
float calc_initial_distance(double *tow_imu, double raw_imu[3], exchange_t *exchange, char get_curr_dist);
#endif

/* This function updates the value of the distance travelled by the vehicle, considering for that purpose the
 * Kalman Filter velocity. To avoid divergences, two status flag, referred to the quality of the current and the
 * previous epoch are required. Distance is only computed when both flags are set to 1. If they are not, last
 * valid position is stored to update the distance count when possible.
 */
void  update_distance(const nav_sol_t *nav_sol, float delta_t, char prev_state, char curr_state, double *distance);

#endif /* GICSRXAUX_H_ */
