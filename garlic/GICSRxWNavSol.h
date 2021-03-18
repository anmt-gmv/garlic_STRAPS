#ifndef GICSRXWNAVSOL_H
#define GICSRXWNAVSOL_H

#include "dbtypedef.h"
#include "GICSRx_defines.h"

/* This function calculates the GARLIC least-squares navigation solution and the IBPL integrity solution.
 * Two operations mode are available: calculating each epoch the satellite positions and velocities with the
 * ephemeris functions, which is time-consuming, or using ephemeris propagation, which provides a very accurate
 * approach for the satellite state vector but optimizing the CPU resources.
 */
#if USE_EPHEMERIS_PROP == 1
char GICSRxWNavSol(gnssdata_t *gnssdata, SP3STRUCT_D *pSp3, lsq_state_t *lsq_state, lsq_matrix_t *lsq_mat,
				   kconf_t *kconf, propephem_t *propephem, nav_sol_t *nav_sol, double init_pos[5], double init_vel[4]);
#else
char GICSRxWNavSol(gnssdata_t *gnssdata, SP3STRUCT_D *pSp3, lsq_state_t *lsq_state, lsq_matrix_t *lsq_mat,
				   kconf_t *kconf, nav_sol_t *nav_sol, double init_pos[5], double init_vel[4]);
#endif

#endif
