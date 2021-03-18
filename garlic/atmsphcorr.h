#ifndef ATMSPHCORR_H_
#define ATMSPHCORR_H_

#include "dbtypedef.h"
#include "GICSRx_defines.h"

#if CALC_ATMSPHCORR == 1

void initialize_NeQuick_iono_model(char *mapping_path, ionoutc_gal_t *iono_gal);

char calculate_atmsphCorr(int month, double day, double tow, double geod_pos[3], ionoutc_t *iono_utc, obsdata_t *obsdata, ionoModel_t iono_model);
void update_NeQuick_table(gnssdata_t *gnssdata, ionoutc_gal_t *iono_utc_gal, double geod_pos[3], int month);
#endif

#endif /* ATMSPHCORR_H_ */
