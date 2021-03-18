#ifndef GICSRXIBPL_H
#define GICSRXIBPL_H

#include "dbtypedef.h"
#include "GICSRx_defines.h"

void compute_dop(const double *position, int num_pars, const float *cov_mat, float *hdop, float *vdop);

#if CALCULATE_IBPL == 1

void initialize_kibpl_indicators(const kconf_t * conf, const lsq_state_t * lsq_state, double tow, kibpl_t * kibpl_state);

void compute_kibpl_indicators(
		double tow,
		unsigned int state_size,
		const kconf_t * kconf,
		float acc_indicator_factor,
		const double * state,
		const double * delta_state,
		const float  * kal_cov_mat,
		kibpl_variables_t * kibpl_state,
		gnssdata_t * gnssdata);

kibpl_code_t GICSRx_KF_IBPL(char flag, const kconf_t * conf, double tow, int kal_used_sats, kibpl_t * kibpl_state, nav_sol_t * nav_sol);
char GICSRx_LSQ_IBPL(lsq_state_t *lsq_state);

#endif

#endif
