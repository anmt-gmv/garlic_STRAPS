#ifndef SUMMARY_RECORDER_H_
#define SUMMARY_RECORDER_H_

#include "dbtypedef.h"
#include "dbdefine.h"

typedef struct
{
	unsigned int rej_invalid_prn;
	unsigned int rej_sv_not_tracked;
	unsigned int rej_sv_not_locked;
	unsigned int rej_out_of_range;
	unsigned int rej_sv_state_not_available;
	unsigned int rej_ephem_not_locked;
	unsigned int rej_low_cn0;
	unsigned int rej_low_elevation;
	unsigned int rej_chisquare_test;

} Rejection_Summary;

typedef struct
{
	unsigned int num_obs;
	unsigned int num_rej_obs;
	double rms;
	double weighted_rms;
	Rejection_Summary rej_summary;

} Summary_Per_Observation_Type;

typedef struct
{
	Summary_Per_Observation_Type code_observations;
	Summary_Per_Observation_Type doppler_observations;

} Observation_Summary;

typedef struct
{
	Observation_Summary sat_summary[MAX_N_SAT_G];
	Observation_Summary global_summary;

} Algorithm_Summary;

void reset_algo_summary(Algorithm_Summary * summary);
void add_rejected_code_observation_to_summary(char SV_PRN, Obs_Status status, Algorithm_Summary * summary);
void add_rejected_doppler_observation_to_summary(char SV_PRN, Obs_Status status, Algorithm_Summary * summary);
void add_rejected_observation_to_summary(char SV_PRN, Obs_Status status, Algorithm_Summary * summary);
void add_code_observation_to_summary(char SV_PRN, float residual, float w_residual, Obs_Status status, Algorithm_Summary * summary);
void add_doppler_observation_to_summary(char SV_PRN, float residual, float w_residual, Obs_Status status, Algorithm_Summary * summary);
void compute_summary(Algorithm_Summary * summary);
void open_summary_file(const char * filename, int index);
void write_summary(double second, const Algorithm_Summary * summary, int index);
void close_summary_file(int index);


#endif /* SUMMARY_RECORDER_H_ */
