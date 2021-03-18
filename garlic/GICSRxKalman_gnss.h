#include "dbtypedef.h"
#include "GICSRx_defines.h"

#ifndef GICSRXKALMAN2_H_
#define GICSRXKALMAN2_H_

/* This function implements a GNSS-only Kalman Filter, which shall be invoked to every epoch in order to perform properly
 * the propagation and update (if new measurements for the current time are provided) operations in the navigation algorithm.
 */
#if CALC_GNSS_ONLY == 1
void kalman_gnss_only(gicsrx_t *gicsrx, gnssdata_t *gnssdata, char result, lsq_state_t *lsq_state, nav_sol_t *nav_sol);
#endif

#endif /* GICSRXKALMAN2_H_ */
