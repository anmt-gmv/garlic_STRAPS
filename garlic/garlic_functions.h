#ifndef GARLIC_TASK_H
#define GARLIC_TASK_H

#include "dbtypedef.h"
#include "calendar.h"
#include "GICSRx_defines.h"
#include "Summary_Recorder.h"

#endif

/* This function is used for the initialization of the GARLIC navigation data.
 * Currently, the configuration parameters for the different algorithms is pre-set.
 */
#if USE_EPHEMERIS_PROP == 1
void initialize_garlic_data(gnssdata_t *gnssdata, gicsrx_t *gicsrxdata, propephem_t *propephem, lsq_state_t *lsq_state, lsq_matrix_t *lsq_mat);
#else
void initialize_garlic_data(gnssdata_t *gnssdata, gicsrx_t *gicsrxdata, lsq_state_t *lsq_state, lsq_matrix_t *lsq_mat, char **argv);
#endif

#if USE_IMU == 0
/* This function is used to calculate the traveled distance when the navigation algorithm
 * works in GNSS-only mode.
 */
void calculate_distance_metrics(gicsrx_t *gicsrx, nav_sol_t *nav_sol, double *distance);
#endif

/* This function is used to store the raw observation data, satellite positions and
 * least-squares navigation solution into the gicsrx_t data structure.
 */
void update_garlic_data(gnssdata_t *gnssdata, lsq_state_t *lsq_state, double geod_pos[3], kconf_t *conf);

/* This function is used to filter the list of measurements acordding to the configuration.
 */
void filter_obs_by_config(gnssdata_t *gnssdata, kconf_t *conf);

void combine_obs(gnssdata_t *gnssdata, kconf_t *conf);

/* This function is used to determine the system corresponding a to a given signal type
 */
int signal2system(sigFlag_t sigFlag);

/* This function is used to determine the frequency band corresponding a to a given signal type
 */
int signal2freq(sigFlag_t sigFlag);

/* This function is used to get internal index from a given prn number
 */
int sys2garlicIndex(int prnIndex, sysFlag_t sysFlag);

/* This function is used to get the prn number from the internal index
 */
int garlic2sysIndex(int prnIndex, sysFlag_t sysFlag);

/* This function is used to get the wavelength for a given signal type
 */
double get_wavelength(sigFlag_t sigFlag, int fbMode, int channelSlot);

char isPilotSignal(sigFlag_t sigFlag);
char isDataSignal(sigFlag_t sigFlag);
