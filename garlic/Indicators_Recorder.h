#ifndef INDICATORS_RECORDER_H_
#define INDICATORS_RECORDER_H_

#include "dbtypedef.h"

typedef struct {
	unsigned int num_code_obs;
	unsigned int num_doppler_obs;
	unsigned int num_outliers;
	double rms_code;
	double rms_doppler;
} Obs_Indicators;

typedef struct {
	Obs_Indicators obs_summary_gps;
	Obs_Indicators obs_summary_glonass;
	Obs_Indicators obs_summary_galileo;
	Obs_Indicators obs_summary_total;
} Indicators;

void reset_indicators(Indicators * indicators);
void add_code_observation(sysFlag_t type, float residual, Obs_Status prange_status, Indicators * indicators);
void add_doppler_observation(sysFlag_t type, float residual, Obs_Status doppler_status, Indicators * indicators);
void compute_indicators(Indicators * indicators);
void open_indicators_file(const char * filename, int index);
void write_indicators_epoch(double second, const Indicators * indicators, char status, int index);
void close_indicators_file(int index);


#endif /* INDICATORS_RECORDER_H_ */
