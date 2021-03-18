#ifndef GICSRXOBS_H
#define GICSRXOBS_H

#include "dbtypedef.h"
#include "GICSRx_defines.h"

#if USE_WEIGHTING_MODEL == 1

/* This function assigns the pseudo-range noise deviation in terms of CN0 and velocity. The
 * value of this noise deviation has been determined using a post-processing statistical model.
 */
void prange_weighting_model(obsdata_t *obs, float velocity, float kvelStopThr, float kvelStopFactor);

/* This function assigns the doppler noise deviation in terms of CN0. The value of this
 * noise deviation has been determined using a post-processing statistical model.
 */
void doppler_weighting_model(obsdata_t *obs);

#endif

/* This function calculates the satellite state vector (position and velocity) for a certain reception time, given the
 * satellite transmission time. The use of satellite ephemeris propagation can also be activated to not invoking every
 * call the ephemeris functions but propagating a node position with a second-order model. In this case, the node satellite
 * positions are also calculated in this function if the node time is out of bounds of the propagation window.
 */
char GICSRxSatPos(double rx_time, double tx_time, obsdata_t *obs, ephdata_t *eph, SP3STRUCT_D *pSp3, char prop_mode, propephem_t *propephem);

/* This function is invoked to calculate the pseudo-range and doppler information related to a specific PRN.
 * This information includes residuals / innovations, line of sights, scale factors for clock bias, satellite
 * state (position, clock and velocity) and noise variances.
 */
char reconstruct_observation(int month, double day, double tow, const kconf_t *kconf, int count, ephdata_t *eph_list, SP3STRUCT_D *pSp3, char prop_mode,
		propephem_t *propephem, float LatLonTrig[4], lsq_state_t *lsq_state, obsdata_t *obsdata, float *Gcol, double *scaled_bias, char *update_eph, ionoutc_t *iono_utc);

#if USE_EPHEMERIS_PROP == 1
/* This function is used to update the ephemeris of the GPS and GLONASS satellites included in the PRN mask.
 * It shall be invoked after the Least-Squares estimation process has been performed, if it is wanted that
 * at least one satellite updates its ephemeris each second, since the LSQ estimation includes implicitly a
 * refresh of ephemeris, if any measurement has become outdated. If not, the currChannel element is refreshed.
 */
void update_ephemeris(gnssdata_t *gnssdata, lsq_state_t *lsq_state, propephem_t *propephem, SP3STRUCT_D *pSp3);
#endif

#if RECEIVER_MODEL == 0	// 0: Teseo-II

/* Given a state vector of 9 components properly ordered, and which clock bias component has been linearly propagated,
 * updates the clock bias value to the current epoch, calculating the drift step, taken into account a pivot observation.
 */
void propagate_clock_bias(double *state_vec, gnssdata_t *gnssdata, lsq_state_t *lsq_state, int maxIter);

#endif

/* This function check that the current observation is valid in terms of being locked, CN0 and elevation.
 * If the criteria are not met, the observation is set as invalid, with the corresponding information.
 */
char check_valid_obs(gicsrx_t *gicsrx, obsdata_t *obs);

/* This function calculates the pseudo-range innovation / residual for a specific observation, considering the state given by the sum
 * of the linear_state and delta_state vectors. In normal operation, the first vector shall provide a fix value of the state, and
 * the second the accumulated innovations throughout the update function in the corresponding Kalman Filter. The function can be
 * used both for the 18-states (Hybrid) and 9-states (GNSS-only) Kalman Filter, properly identified with num_states input parameter.
 * The only consideration required is that the satellite state (position, clock and velocity) has been previously computed.
 * Due to the possibility of working with synchronized IMU information, it may be possible that the solution time does not correspond
 * exactly with the GNSS time. This difference is obtained by the difference of the input variables gnsstow and kaltow, and is taken
 * into account in order to calculate the navigation solution consistently.
 */
void compute_prange_residual(const double *linear_state, const double *delta_state, const obsdata_t *obs, double gnsstow, double kaltow, float *residual, float *h_vec, int numStates);

/* This function calculates the doppler innovation / residual for a specific observation, considering the state given by the sum
 * of the linear_state and delta_state vectors. In normal operation, the first vector shall provide a fix value of the state, and
 * the second the accumulated innovations throughout the update function in the corresponding Kalman Filter. The function can be
 * used both for the 18-states (Hybrid) and 9-states (GNSS-only) Kalman Filter, properly identified with num_states input parameter.
 * The only consideration required is that the satellite state (position, clock and velocity) has been previously computed.
 */
void compute_doppler_residual(const double *linear_state, const double *delta_state, const obsdata_t *obs, float *residual, float *h_vec, int numStates);

#endif
