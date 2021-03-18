/*
 * GICSRxMechanization.h
 *
 *  Created on: 10/05/2013
 *      Author: cmvv
 */

#ifndef GICSRXMECHANIZATION_H_
#define GICSRXMECHANIZATION_H_

#include <stdio.h>

#include "dbtypedef.h"
#include "GICSRx_defines.h"

#if USE_IMU == 1

void read_sensor_file(double tow, sensor_t imu[MAX_IMU_SAMPLES], int *nsamples, kconf_t *kconf);

/* This function implements a Hybrid Kalman Filter and shall be invoked every epoch to perform properly the
 * propagation and update functions. In the Hybrid KF, the propagation functionality implements a mechanization
 * procedure, where data from inertial sensors (accelerometers and gyros) is incorporated into the navigation
 * system. The samples from IMU shall be read by an external function, since this function takes them as input.
 * Additional features are included, apart from navigation, as for instance consistency check of the filter,
 * distance computation or anti-spoofing algorithm (if activated).
 */
void kalman_hybrid(gicsrx_t *gicsrx, gnssdata_t *gnssdata, lsq_state_t *lsq_state, sensor_t imu[MAX_IMU_SAMPLES], int nsamples, nav_sol_t *nav_sol);

#endif /* GICSRXMECHANIZATION_H_ */

#endif
