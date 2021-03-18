/////////////////////////////////////////////////////////////////////////////////
// Project:                     GICSRX
// Purpose:                     Functions for observations reconstruction.
// File:                        GICSRxObs.c
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
// |   1.00  | 13/09/26 | cmvv         | Initial Release                      |
// +--------------------------------------------------------------------------+
//
// DESCRIPTION: this file includes the definitions of the least-squares function,
//              atmospheric correction functions and IBPL integrity functions.
//
//
/////////////////////////////////////////////////////////////////////////////////
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"
#include "algebra.h"
#include "atmsphcorr.h"
#include "GICSRxObs.h"
#include "GICSRxPosition.h"
#include "garlic_functions.h"

#if USE_WEIGHTING_MODEL == 1

/* This function assigns the pseudo-range noise deviation in terms of CN0 and velocity. The
 * value of this noise deviation has been determined using a post-processing statistical model.
 */
void prange_weighting_model(obsdata_t *obs, float velocity, float kvelStopThr, float kvelStopFactor)
{
	int index = (int)obs->S1 - 15;
	if (index < 0) {
		index = 0;
	} else if (index > 41) {
		index = 41;
	}
	obs->prange_sigma = ((velocity >= kvelStopThr) ? sigmas_PR[index] : kvelStopFactor*sigmas_PR[index]);

#if	RECEIVER_MODEL == 0
	// Experimental increase factor in GLONASS satellites for the pseudo-range measurement sigma.
	// Different values for the multiplication factor can be tested. Good accuracy results with
	// a factor of 3, although the performance of the distance computation algorithm is worse.
	if(signal2system(obs->sigFlag) == SYS_GLO) { obs->prange_sigma *= 3.00; }
#endif
}

/* This function assigns the doppler noise deviation in terms of CN0 and velocity. The
 * value of this noise deviation has been determined using a post-processing statistical model.
 */
void doppler_weighting_model(obsdata_t *obs)
{
	int index = (int)obs->S1 - 15;
	if (index < 0) {
		index = 0;
	} else if (index > 41) {
		index = 41;
	}
	obs->doppler_sigma = sigmas_DP[index];

#if	RECEIVER_MODEL == 0
	// Experimental increase factor in GLONASS satellites for the doppler measurement sigma.
	if(signal2system(obs->sigFlag) == SYS_GLO) { obs->doppler_sigma *= 3.00; }
#endif
}

#endif

/* This function calculates the satellite state vector (position and velocity) for a certain reception time, given the
 * satellite transmission time. The use of satellite ephemeris propagation can also be activated to not invoking every
 * call the ephemeris functions but propagating a node position with a second-order model. In this case, the node satellite
 * positions are also calculated in this function if the node time is out of bounds of the propagation window.
 */
char GICSRxSatPos(double rx_time, double tx_time, obsdata_t *obs, ephdata_t *eph, SP3STRUCT_D *pSp3, char prop_mode, propephem_t *propephem)
{
	char   update_eph = 0;

	int    prop_id = 0;
	double auxdata[3], acc[3], r, r3;
	double dt, cosrot, sinrot;

	// Pointer to the satellite state vector.
	ephemsat_t *pephsat;

	// If there is no ephemeris propagation data structure, turn off propagation mode.
	if(propephem == NULL)
	{
		prop_mode = 0;
	}
	else
	{
		// Calculate the time gap between the node and the current epoch.
		dt      = rx_time - propephem->tow_node[obs->PRN-1];
		prop_id = (int)(dt + (dt > 0 ? +1:-1) * 0.5);
	}

	// Check whether current epoch is out of bounds of the propagation window.
	// In this case, calculate again node (central propagation window) satellite position.
	// Also, if propephem is NULL, no propagation is used, thus satellite position must be
	// calculated.
	if(propephem == NULL || fabs(prop_id) >= PROP_WINDOW/2)
	{
		// Propagation mode de-activated. Pointer to current observation data.
		if(prop_mode == 0)
		{
			pephsat = &(obs->sat);
		}
		// Propagation mode activated. Shift current reception time to the center of the
		// propagation window.
		else
		{
			// Move current time to the next propagation window.
			rx_time  += PROP_WINDOW/2;
			tx_time  += PROP_WINDOW/2;
			// Get pointer to the satellite node state vector, used for propagation.
			pephsat   = &(propephem->sat[obs->PRN-1]);

			// Assign new time of week for the window node.
			propephem->tow_node[obs->PRN-1] = rx_time;
			update_eph = 1;
		}

		if(pSp3 != NULL)
		{
			GICSRxSatPos_SP3(pSp3, garlic2sysIndex(obs->PRN-1,signal2system(obs->sigFlag))+1, rx_time, tx_time, pephsat);
		}
		else
		{
			double bgd = 0;
			// Calculate satellite position and velocity in RxTime.
			switch(signal2system(obs->sigFlag))
			{
			// GPS satellite.
			case SYS_GPS:
				// Use GPS received ephemeris. Position, velocity and clock bias.
				GICSRxGPSPos(&(eph->GPS), rx_time, tx_time, pephsat, 1);
				break;

			// GLONASS satellite.
			case SYS_GLO:
				// Use GLONASS received ephemeris. Position, velocity and clock bias.
				GICSRxGLOPos(&(eph->GLO), rx_time, tx_time, pephsat, 1);
				break;

			// Galileo satellite.
			case SYS_GAL:

				// In case of iono-free combination there is not need to correct BGDs
				if (obs->iono_free == 0)
				{
					switch(signal2freq(obs->sigFlag))
					{
					case FREQ_L1E1:
						bgd = eph->GAL.bgd_E1E5b;
						break;
					case FREQ_E5a:
						bgd = eph->GAL.bgd_E1E5a * (BANDFREQL1*BANDFREQL1)/(BANDFREQE5a*BANDFREQE5a);
						break;
					case FREQ_E5b:
						bgd = eph->GAL.bgd_E1E5b * (BANDFREQL1*BANDFREQL1)/(BANDFREQE5b*BANDFREQE5b);
						break;
					case FREQ_E6:
						bgd = eph->GAL.bgd_E1E6A;
						break;
					}
				}
				// Use Galileo received ephemeris. Position, velocity and clock bias.
				GICSRxGalPos(&(eph->GAL), rx_time, tx_time, pephsat, bgd);
				break;

			// BeiDou satellite.
			case SYS_BEI:
				// Use GPS received ephemeris. Position, velocity and clock bias.
				GICSRxBEIPos(&(eph->BEI), rx_time, tx_time, pephsat, 1);
				break;

			// Unknown satellite type.
			default:
				break;
			}
		}

		// Go back to current week second if a future position has been calculated.
		if(prop_mode == 1)
		{
			rx_time -= PROP_WINDOW/2;
			tx_time -= PROP_WINDOW/2;
		}
	}

	// If propagation is activated (in mode 1), propagate current window position.
	// The function can also be invoked with a value not 0 nor 1, for instance, 2. This is
	// useful to recalculate the node timestamp without providing the current position.
	if(prop_mode == 1)
	{
		// Obtain window node for this prn and calculate the time difference between
		// current epoch and node time.
		pephsat = &(propephem->sat[obs->PRN-1]);
		dt      = rx_time - propephem->tow_node[obs->PRN-1];

		// Rotation matrix elements (rotate from time in node to current reception time).
		cosrot = cos(OMEGA_EARTH*dt);
		sinrot = sin(OMEGA_EARTH*dt);

		// Calculate the norm of position vector, and its 3-power, to calculate gravity vector.
		r  = Vec3Norm(pephsat->pos);
		r3 = r*r*r;

		acc[0] = -MU_EARTH/r3*pephsat->pos[0];
		acc[1] = -MU_EARTH/r3*pephsat->pos[1];
		acc[2] = -MU_EARTH/r3*pephsat->pos[2];

		// Calculate current position and velocity from the node, considering velocity and position
		// constants during the propagation period. Take also into account the appropriate rotations.
		auxdata[0] = pephsat->pos[0] + dt*pephsat->vel[0] + 0.5*dt*dt*acc[0];
		auxdata[1] = pephsat->pos[1] + dt*pephsat->vel[1] + 0.5*dt*dt*acc[1];
		auxdata[2] = pephsat->pos[2] + dt*pephsat->vel[2] + 0.5*dt*dt*acc[2];

		// Correct earth rotation (variation in ECEF frame from node time to reception time) in position.
		obs->sat.pos[0] = +cosrot*auxdata[0] + sinrot*auxdata[1];
		obs->sat.pos[1] = -sinrot*auxdata[0] + cosrot*auxdata[1];
		obs->sat.pos[2] = auxdata[2];

		auxdata[0] = pephsat->vel[0] + dt*acc[0];
		auxdata[1] = pephsat->vel[1] + dt*acc[1];
		auxdata[2] = pephsat->vel[2] + dt*acc[2];

		// Correct earth rotation (variation in ECEF frame from node time to reception time) in velocity.
		obs->sat.vel[0] = +cosrot*auxdata[0] + sinrot*auxdata[1];
		obs->sat.vel[1] = -sinrot*auxdata[0] + cosrot*auxdata[1];
		obs->sat.vel[2] = auxdata[2];

		// Propagate clock elements.
		obs->sat.clk_bias  = pephsat->clk_bias  + dt*pephsat->clk_drift + 0.5*dt*dt*pephsat->clk_driftrate;
		obs->sat.clk_drift = pephsat->clk_drift + dt*pephsat->clk_driftrate;

		// Store clock drift rate.
		obs->sat.clk_driftrate = pephsat->clk_driftrate;
	}
	return update_eph;
}

/* This function is invoked to calculate the pseudo-range and doppler information related to a specific PRN.
 * This information includes residuals / innovations, line of sights, scale factors for clock bias, satellite
 * state (position, clock and velocity) and noise variances.
 */
char reconstruct_observation(int month, double day, double tow, const kconf_t *kconf, int count, ephdata_t *eph_list, SP3STRUCT_D *pSp3, char prop_mode,
		propephem_t *propephem, float LatLonTrig[4], lsq_state_t *lsq_state, obsdata_t *obsdata, float *Gcol, double *scaled_bias, char *update_eph, ionoutc_t *iono_utc)
{
	int i;
	char sys_clk[NUM_SYS] = {0};
	double clk_bias;

	// Check that satellite is locked and ephemeris is available.
	if (obsdata->prange_status == REJ_SV_NOT_LOCKED && obsdata->doppler_status == REJ_SV_NOT_LOCKED) {
		return 0;
	}
	if (obsdata->lock == 0)
	{
		obsdata->prange_status  = REJ_EPHEM_NOT_LOCKED;
		obsdata->doppler_status = REJ_EPHEM_NOT_LOCKED;
		return 0;
	}

	// Check that prn is inside array bounds.
	if(garlic2sysIndex(obsdata->PRN-1, signal2system(obsdata->sigFlag)) < 0)
	{
		obsdata->prange_status  = REJ_INVALID_PRN;
		obsdata->doppler_status = REJ_INVALID_PRN;
		return 0;
	}

	// Check which system clock bias shall be considered and correct the tx_time accordingly
	sysFlag_t sys = signal2system(obsdata->sigFlag);
	sys_clk[sys] = 1;
	clk_bias = lsq_state->pos[3+sys] + FSTX_ISBIAS[sys];

	// Calculate transmission time.
//	double tx_time = tow - (obsdata->C1 + obsdata->sat.clk_bias)/SPEED_OF_LIGHT;
//	double tx_time = tow - (obsdata->C1 - (*scaled_bias + clk_bias) + obsdata->sat.clk_bias)/SPEED_OF_LIGHT;

	tow -= clk_bias/SPEED_OF_LIGHT;
	double tx_time = tow - (obsdata->C1 - (*scaled_bias + clk_bias) + obsdata->sat.clk_bias)/SPEED_OF_LIGHT;

#if USE_EPHEMERIS_PROP == 1
	// Calculate the position and velocity of the satellite at reception time. Use ephemeris propagation.
	*update_eph = GICSRxSatPos(tow, tx_time, obsdata, &eph_list[obsdata->PRN-1], pSp3, prop_mode, propephem);
#else
	// Calculate the position and velocity of the satellite at reception time. No propagation.
	*update_eph = GICSRxSatPos(tow, tx_time, obsdata, &eph_list[obsdata->PRN-1], pSp3, 0, NULL);
#endif

	// Computes PR and unit vector after corrections.
	float  rsu[3];
	double geom_range = GICSRxSatRange(lsq_state->pos, obsdata->sat.pos, rsu);

	// Computes elevation and azimuth angles.
#if CALC_ATMSPHCORR == 1
	GICSRxSatElevAzim(LatLonTrig, rsu, &obsdata->elev, &obsdata->azim);
#else
	GICSRxSatElev(LatLonTrig, rsu, &obsdata->elev);
	obsdata->azim = 0.0;
#endif

	// Store row in the geometry matrix.
	int j;
	for(j = 0; j < 3; j++) { Gcol[j] = -rsu[j]; }

	// Terms corresponding to system clock biases.
	for(i=0; i<NUM_SYS; i++) { Gcol[3+i] = sys_clk[i]; }

	// Terms corresponding to inter-frequency biases.
	double ifreq_bias = 0;
	freqFlag_t freq = signal2freq(obsdata->sigFlag);
	if (freq > 0)
	{
		Gcol[3+NUM_SYS+freq-1] = 1;
		ifreq_bias = lsq_state->pos[3+NUM_SYS+freq-1];
	}

#if CALC_ATMSPHCORR == 1
	if(count > 0)
	{
		double geod_pos[3];
		ECEFtoNAV_pos(lsq_state->pos, geod_pos);
		calculate_atmsphCorr(month, day, tow, geod_pos, iono_utc, obsdata, kconf->iono_model);
	}
#endif

	// Calculate the additional bias over the pseudo-range measurement (0.0*c, 0.2*c or 0.4*c in Teseo II).
	*scaled_bias = ((int)((obsdata->C1 - geom_range)/(0.2*SPEED_OF_LIGHT) + 0.5))*0.2*SPEED_OF_LIGHT;

	// Residual range to be minimized is calculated.
	obsdata->prange_residual = obsdata->C1 - geom_range - (*scaled_bias + clk_bias) - ifreq_bias + obsdata->sat.clk_bias - obsdata->atmsphCorr;

	if(count > 0)
	{
		// Residual velocity to be minimized is computed.
		double vel_doppler = GICSRxSatDoppler(lsq_state->vel, obsdata->sat.vel, rsu, lsq_state->pos);

		double clk_drift = lsq_state->vel[3];
		if (TWO_RX == 1 && (sys == SYS_GLO || sys == SYS_BEI))
		{
			clk_drift = lsq_state->vel[4];
		}
		// Calculate zvel measure for current satellite. Lambda value is defined for f1.
		obsdata->doppler_residual = obsdata->D1 - vel_doppler - clk_drift + obsdata->sat.clk_drift;
	}

	// Check whether the calculated results are numerically valid.
	if(obsdata->S1 < kconf->CN0_thr)
	{
		obsdata->prange_status  = REJ_LOW_CN0;
		obsdata->doppler_status = REJ_LOW_CN0;
		return 0;
	}

	// Apply elevation mask.
	if(count > 0 && obsdata->elev < kconf->mask_angle)
	{
		obsdata->prange_status  = REJ_LOW_ELEV;
		obsdata->doppler_status = REJ_LOW_ELEV;
		return 0;
	}

	if (isnan(obsdata->prange_residual))
	{
		obsdata->prange_status = REJ_OUT_OF_RANGE;
	}

	if (isnan(obsdata->doppler_residual))
	{
		obsdata->doppler_status = REJ_OUT_OF_RANGE;
	}

	if (obsdata->prange_status > OBS_OK_BUT_UNUSED)
	{
		return 0;
	}

	return 1;

} // END of function reconstruct_observation

#if USE_EPHEMERIS_PROP == 1
/* This function is used to update the ephemeris of the GPS and GLONASS satellites included in the PRN mask.
 * It shall be invoked after the Least-Squares estimation process has been performed, if it is wanted that
 * at least one satellite updates its ephemeris each second, since the LSQ estimation includes implicitly a
 * refresh of ephemeris, if any measurement has become outdated. If not, the currChannel element is refreshed.
 */
void update_ephemeris(gnssdata_t *gnssdata, lsq_state_t *lsq_state, propephem_t *propephem, SP3STRUCT_D *pSp3)
{

	// Check that pointer to the channel to be updated is not out of bounds.
	if(propephem->currChannel < 0 || propephem->currChannel >= gnssdata->noOfChannelsAv)
	{
		propephem->currChannel = 0;
	}

	// Refresh only an available channel not recently updated.
	//while(propephem->currChannel < gnssdata->noOfChannelsAv && gnssdata->OBS[propephem->currChannel].lock == 0)
	while(propephem->currChannel < gnssdata->noOfChannelsAv && gnssdata->OBS[propephem->currChannel].prange_status != OBS_VALID)
	{
		propephem->currChannel++;
	}

	// If there is an available channel, update the ephemeris matrix.
	if(propephem->currChannel < gnssdata->noOfChannelsAv)
	{
		// Get observation channel id and position to the ephemeris.
		int i   = propephem->currChannel;
		int prn = (gnssdata->OBS[i].PRN - 1);

		// Check which system clock bias shall be considered and correct the
		// tx_time accordingly
		sysFlag_t sys = signal2system(gnssdata->OBS[i].sigFlag);
		double clk_bias = lsq_state->pos[3+sys] + FSTX_ISBIAS[sys];

		// Calculate the transmission time for the selected satellite.
		double tx_time = gnssdata->tow - (gnssdata->OBS[i].C1 - (lsq_state->scaled_bias[i] + clk_bias) + gnssdata->OBS[i].sat.clk_bias)/SPEED_OF_LIGHT;

		// Set the ephemeris matrix of the selected satellite as unavailable.
		propephem->tow_node[prn] = -DAYSECS;

		// Call the refresh function.
		// FIXME: Revisit flag to decide prop mode
		GICSRxSatPos(gnssdata->tow, tx_time, &(gnssdata->OBS[i]), &(gnssdata->EPH[prn]), pSp3,
				2*(lsq_state->num_LOS[sys] > 4), propephem);
	}
	// Increase channel counter.
	(propephem->currChannel)++;
}
#endif

#if RECEIVER_MODEL == 0

static int calculate_pivot_satellite(const gnssdata_t *gnssdata)
{
	int k, max_cn0_ind = -1, max_cn0_val = 0;

	for(k = 0; k < gnssdata->noOfChannelsAv; k++)
	{
		if(gnssdata->OBS[k].prange_status == OBS_VALID && gnssdata->OBS[k].CN0_L1 > max_cn0_val)
		{
			max_cn0_ind = k;
			max_cn0_val = gnssdata->OBS[k].CN0_L1;
		}
	}
	return max_cn0_ind;
}

/* This function is specific to be used with Teseo II. It is known that the Teseo II clock bias includes
 * a scale factor which is function of the speed of light (multiplied by 0, 0.2 or 0.4). Moreover, the
 * evolution in time of the clock bias does not correspond to a linear behavior, but there also exists
 * a step in the clock drift, that shall be corrected prior to update the Kalman Filter, in order to
 * obtain an accurate propagation.
 * It has not been observed any deterministic pattern associated with this drift step, although the
 * values that this step can take do (in particular, three different values). Thus, the most suitable
 * step is calculated and corrected in the clock bias term, after propagating with drift.
 * The state vector shall be provided in ECEF coordinates.
 */
static double correct_drift_step(double *state_vec, obsdata_t *obs, char *stop_criterion)
{
const double max_step_discrepancy = 500;

	char use_isb_glo, use_isb_gal;
	double range, drift_step;
	double diff_0, diff_1, diff_2;

	float rsu[3];

	// Calculate distance between user and current satellite (PR estimated).
	range = GICSRxSatRange(state_vec, obs->sat.pos, rsu);

	use_isb_glo = (obs->sysFlag == SYS_GLO);
	use_isb_gal = (obs->sysFlag == SYS_GALOS || obs->sysFlag == SYS_GALPRS);

	// Obtain clock drift step value plus measurement noise.
	drift_step = obs->C1 - range - state_vec[6] - use_isb_glo*(FSTX_ISBIASGLO + state_vec[7]) - use_isb_gal*(FSTX_ISBIASGAL + state_vec[7+KAC]) + obs->sat.clk_bias;

	diff_0 = fabs(drift_step - drift_values[0]);
	diff_1 = fabs(drift_step - drift_values[1]);
	diff_2 = fabs(drift_step - drift_values[2]);

	// Check whether is the most suitable step value.
	if(diff_0 < diff_1 && diff_0 < diff_2)
	{
		drift_step = drift_values[0];
		if(diff_0 < max_step_discrepancy) { *stop_criterion = 1; }
	}
	else if(diff_1 < diff_0 && diff_1 < diff_2)
	{
		drift_step = drift_values[1];
		if(diff_1 < max_step_discrepancy) { *stop_criterion = 1; }
	}
	else
	{
		drift_step = drift_values[2];
		if(diff_2 < max_step_discrepancy) { *stop_criterion = 1; }
	}
	return drift_step;
}

/* Given a state vector of 9 components properly ordered, and which clock bias component has been linearly propagated,
 * updates the clock bias value to the current epoch, calculating the drift step, taken into account a pivot observation.
 */
void propagate_clock_bias(double *state_vec, gnssdata_t *gnssdata, lsq_state_t *lsq_state, int maxIter)
{
	char stop_criterion = 0, iterations = 0;
	int pivot = calculate_pivot_satellite(gnssdata);

	// Correct drift step with the help of the pivot satellite in the mask.
	while(pivot >= 0 && stop_criterion == 0 && iterations < maxIter)
	{
		state_vec[6] += correct_drift_step(state_vec, &gnssdata->OBS[pivot], &stop_criterion);
		iterations++;
	}

	if(iterations == maxIter && lsq_state->pos_valid == 1)
	{
		state_vec[6] = lsq_state->pos[3];
	}
}

#endif

/* This function check that the current observation is valid in terms of being locked, CN0 and elevation.
 * If the criteria are not met, the observation is set as invalid, with the corresponding information.
 */
char check_valid_obs(gicsrx_t *gicsrx, obsdata_t *obs)
{
	// Check whether the measurement is locked, tracked, and ephemeris are available.
	if(obs->lock == 0)
	{
		obs->prange_status  = REJ_EPHEM_NOT_LOCKED;
		obs->doppler_status = REJ_EPHEM_NOT_LOCKED;
		return 0;
	}

	// Check whether CN0 is above the configured mask value.
	if(obs->S1 < gicsrx->kconf.CN0_thr)
	{
		obs->prange_status  = REJ_LOW_CN0;
		obs->doppler_status = REJ_LOW_CN0;
		return 0;
	}

	// Check whether the elevation is above the configuration mask value.
	if(obs->elev < gicsrx->kconf.mask_angle)
	{
		obs->prange_status  = REJ_LOW_ELEV;
		obs->doppler_status = REJ_LOW_ELEV;
		return 0;
	}
	return 1;
}

/* This function calculates the pseudo-range innovation / residual for a specific observation, considering the state given by the sum
 * of the linear_state and delta_state vectors. In normal operation, the first vector shall provide a fix value of the state, and
 * the second the accumulated innovations throughout the update function in the corresponding Kalman Filter. The function can be
 * used both for the 18-states (Hybrid) and 9-states (GNSS-only) Kalman Filter, properly identified with num_states input parameter.
 * The only consideration required is that the satellite state (position, clock and velocity) has been previously computed.
 * Due to the possibility of working with synchronized IMU information, it may be possible that the solution time does not correspond
 * exactly with the GNSS time. This difference is obtained by the difference of the input variables gnsstow and kaltow, and is taken
 * into account in order to calculate the navigation solution consistently.
 */
void compute_prange_residual(const double *linear_state, const double *delta_state, const obsdata_t *obs, double gnsstow, double kaltow, float *residual, float *h_vec, int numStates)
{

	// Declaration of variables.
	int i;
	float rsu[3];

	double delta_t = gnsstow - kaltow;

	double userpos[3] = {(linear_state[0] + delta_state[0]) + delta_t*(linear_state[3] + delta_state[3]),
						 (linear_state[1] + delta_state[1]) + delta_t*(linear_state[4] + delta_state[4]),
						 (linear_state[2] + delta_state[2]) + delta_t*(linear_state[5] + delta_state[5])};

	// Initialize observation vector H.
	for (i = 0; i < numStates; i++) { h_vec[i] = 0.0; }

	// Calculate distance between user and current satellite (PR estimated).
	double range = GICSRxSatRange(userpos, obs->sat.pos, rsu);

	// Build observation vector.
	for (i = 0; i < 3; i++)
	{
		h_vec[i]   = -rsu[i];
		h_vec[i+3] = -rsu[i]*delta_t;
	}

	double rx_clock_bias = (linear_state[6] + delta_state[6]) + delta_t*(linear_state[6+NUM_SYS+NUM_FREQ-1] + delta_state[6+NUM_SYS+NUM_FREQ-1]);

	int sys = signal2system(obs->sigFlag);
	double isb = 0.0;
	if (sys > 0)
	{
		isb = FSTX_ISBIAS[sys] + linear_state[6+sys] + delta_state[6+sys];
	}

	// Clock, ISBs and Drift
	h_vec[6] = 1;
	if (sys > 0)
	{
		h_vec[6+sys] = 1;
	}
	h_vec[6+NUM_SYS+NUM_FREQ-1] = delta_t;

	// Inter-frequency bias
	int freq = signal2freq(obs->sigFlag);
	double ifreq_bias = 0.0;
	if (freq > 0)
	{
		h_vec[6+NUM_SYS+freq-1] = 1;
		ifreq_bias = linear_state[6+NUM_SYS+freq-1] + delta_state[6+NUM_SYS+freq-1];
	}

	if (TWO_RX == 1 && (sys == SYS_GLO || sys == SYS_BEI))
	{
		// Change reference clock
		h_vec[6]         = 0;
		h_vec[6+SYS_GLO] = 1;
		if (sys == SYS_GLO) { isb = 0.0;}

		// Recompute clock bias and isb with new reference clock
		rx_clock_bias = (linear_state[6+SYS_GLO] + delta_state[6+SYS_GLO]) + delta_t*(linear_state[6+NUM_SYS+NUM_FREQ-1+TWO_RX] + delta_state[6+NUM_SYS+NUM_FREQ-1+TWO_RX]);
	}

	*residual = obs->C1 - range - rx_clock_bias + obs->sat.clk_bias - isb - ifreq_bias;

} // END of function compute_prange_residual

/* This function calculates the doppler innovation / residual for a specific observation, considering the state given by the sum
 * of the linear_state and delta_state vectors. In normal operation, the first vector shall provide a fix value of the state, and
 * the second the accumulated innovations throughout the update function in the corresponding Kalman Filter. The function can be
 * used both for the 18-states (Hybrid) and 9-states (GNSS-only) Kalman Filter, properly identified with num_states input parameter.
 * The only consideration required is that the satellite state (position, clock and velocity) has been previously computed.
 */
void compute_doppler_residual(const double *linear_state, const double *delta_state, const obsdata_t *obs, float *residual, float *h_vec, int numStates)
{

	int i;
	float rsu[3], user_vel[3];

	double userpos[3] = {linear_state[0] + delta_state[0],
					 	 linear_state[1] + delta_state[1],
					 	 linear_state[2] + delta_state[2]};

	// Initialize observation vector H.
	for (i = 0; i < numStates; i++) { h_vec[i] = 0.0; }

	// Copy velocity vector.
	user_vel[0] = (linear_state[3] + delta_state[3]) - OMEGA_EARTH*userpos[1];
	user_vel[1] = (linear_state[4] + delta_state[4]) + OMEGA_EARTH*userpos[0];
	user_vel[2] = (linear_state[5] + delta_state[5]);

	// Calculate rsu vector (line of sight).
	GICSRxSatRange(userpos, obs->sat.pos, rsu);

	// Build observation vector.
	for (i = 0; i < 3; i++) { h_vec[i+3] = -rsu[i]; }

	// Calculate user and satellite Doppler speed.
	float user_doppler = -(rsu[0]*user_vel[0]     + rsu[1]*user_vel[1]     + rsu[2]*user_vel[2]);
	float sat_doppler  = -(rsu[0]*obs->sat.vel[0] + rsu[1]*obs->sat.vel[1] + rsu[2]*obs->sat.vel[2]);

	h_vec[6+NUM_SYS+NUM_FREQ-1] = 1;
	double clk_drift = linear_state[6+NUM_SYS+NUM_FREQ-1] + delta_state[6+NUM_SYS+NUM_FREQ-1];

	if (TWO_RX == 1)
	{
		int sys = signal2system(obs->sigFlag);
		if (sys == SYS_GLO || sys == SYS_BEI)
		{
			h_vec[6+NUM_SYS+NUM_FREQ-1]        = 0;
			h_vec[6+NUM_SYS+NUM_FREQ-1+TWO_RX] = 1;
			clk_drift = linear_state[6+NUM_SYS+NUM_FREQ-1+TWO_RX] + delta_state[6+NUM_SYS+NUM_FREQ-1+TWO_RX];
		}
	}

	*residual = obs->D1 + sat_doppler - user_doppler - clk_drift + obs->sat.clk_drift;

} // END of function compute_doppler_residual

