#include "Summary_Recorder.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
/*
static FILE * fid_[2] = {NULL};

static void add_rejected_observation(Obs_Status status, Summary_Per_Observation_Type * summary)
{
	if (status != OBS_VALID) {
		if (status == REJ_INVALID_PRN) {
			summary->rej_summary.rej_invalid_prn++;
		} else if (status == REJ_SV_NOT_TRACKED) {
			summary->rej_summary.rej_sv_not_tracked++;
		} else if (status == REJ_SV_NOT_LOCKED) {
			summary->rej_summary.rej_sv_not_locked++;
		} else if (status == REJ_OUT_OF_RANGE) {
			summary->rej_summary.rej_out_of_range++;
		} else if (status == REJ_SV_STATE_NOT_AVAILABLE) {
			summary->rej_summary.rej_sv_state_not_available++;
		} else if (status == REJ_EPHEM_NOT_LOCKED) {
			summary->rej_summary.rej_ephem_not_locked++;
		} else if (status == REJ_CHISQUARE_TEST) {
			summary->rej_summary.rej_chisquare_test++;
			summary->num_rej_obs++;
		} else if (status == REJ_LOW_CN0) {
			summary->rej_summary.rej_low_cn0++;
		} else if (status == REJ_LOW_ELEV) {
			summary->rej_summary.rej_low_elevation++;
		}
	}
}

static void add_observation(float residual, float w_residual, Obs_Status status, Summary_Per_Observation_Type * summary)
{
	if (status == OBS_VALID) {
		summary->num_obs++;
		summary->rms += residual * residual;
		summary->weighted_rms += w_residual * w_residual;
	} else {
		add_rejected_observation(status, summary);
	}
}

static void compute_summary_per_type(Summary_Per_Observation_Type * summary)
{
	if (summary->num_obs > 0) {
		summary->rms = sqrt(summary->rms / summary->num_obs);
		summary->weighted_rms = sqrt(summary->weighted_rms / summary->num_obs);
	}
}

static void reset_obs_summary_per_type(Summary_Per_Observation_Type * summary)
{
	summary->num_obs = 0;
	summary->num_rej_obs = 0;
	summary->rms = 0.0;
	summary->weighted_rms = 0.0;
	summary->rej_summary.rej_invalid_prn = 0;
	summary->rej_summary.rej_sv_not_tracked = 0;
	summary->rej_summary.rej_sv_not_locked = 0;
	summary->rej_summary.rej_out_of_range = 0;
	summary->rej_summary.rej_sv_state_not_available = 0;
	summary->rej_summary.rej_ephem_not_locked = 0;
	summary->rej_summary.rej_chisquare_test = 0;
	summary->rej_summary.rej_low_cn0 = 0;
	summary->rej_summary.rej_low_elevation = 0;
}

static void reset_obs_summary(Observation_Summary * summary)
{
	reset_obs_summary_per_type(&summary->code_observations);
	reset_obs_summary_per_type(&summary->doppler_observations);
}

void reset_algo_summary(Algorithm_Summary * summary)
{
	reset_obs_summary(&summary->global_summary);
	unsigned int i;
	for (i = 0; i < MAX_N_SAT_G; i++) {
		reset_obs_summary(&summary->sat_summary[i]);
	}
}

void add_rejected_code_observation_to_summary(char SV_PRN, Obs_Status status, Algorithm_Summary * summary)
{
	add_rejected_observation(status, &summary->global_summary.code_observations);
	if ((SV_PRN >= MINPRNGPS_G) && (SV_PRN <= MAXPRNGALPRS_G)) {
		add_rejected_observation(status, &summary->sat_summary[SV_PRN-1].code_observations);
	}
}

void add_rejected_doppler_observation_to_summary(char SV_PRN, Obs_Status status, Algorithm_Summary * summary)
{
	add_rejected_observation(status, &summary->global_summary.doppler_observations);
	if ((SV_PRN >= MINPRNGPS_G) && (SV_PRN <= MAXPRNGALPRS_G)) {
		add_rejected_observation(status, &summary->sat_summary[SV_PRN-1].doppler_observations);
	}
}

void add_rejected_observation_to_summary(char SV_PRN, Obs_Status status, Algorithm_Summary * summary)
{
	add_rejected_code_observation_to_summary(SV_PRN, status, summary);
	add_rejected_doppler_observation_to_summary(SV_PRN, status, summary);
}

void add_code_observation_to_summary(char SV_PRN, float residual, float w_residual, Obs_Status status, Algorithm_Summary * summary)
{
	if ((SV_PRN >= MINPRNGPS_G) && (SV_PRN <= MAXPRNGALPRS_G)) {
		add_observation(residual, w_residual, status, &summary->global_summary.code_observations);
		add_observation(residual, w_residual, status, &summary->sat_summary[SV_PRN-1].code_observations);
	}
}

void add_doppler_observation_to_summary(char SV_PRN, float residual, float w_residual, Obs_Status status, Algorithm_Summary * summary)
{
	if ((SV_PRN >= MINPRNGPS_G) && (SV_PRN <= MAXPRNGALPRS_G)) {
		add_observation(residual, w_residual, status, &summary->global_summary.doppler_observations);
		add_observation(residual, w_residual, status, &summary->sat_summary[SV_PRN-1].doppler_observations);
	}
}

void compute_summary(Algorithm_Summary * summary)
{
	compute_summary_per_type(&summary->global_summary.code_observations);
	compute_summary_per_type(&summary->global_summary.doppler_observations);
	unsigned int i;
	for (i = 0; i < MAX_N_SAT_G; i++) {
		compute_summary_per_type(&summary->sat_summary[i].code_observations);
		compute_summary_per_type(&summary->sat_summary[i].doppler_observations);
	}
}

void open_summary_file(const char * filename, int index)
{
	fid_[index] = fopen(filename, "wt");
}

void write_summary(double second, const Algorithm_Summary * summary, int index)
{
	if (fid_[index] != NULL) {
		fprintf(fid_[index], " Measurements statistics\n");
		fprintf(fid_[index], "########################\n\n");

		fprintf(fid_[index], " Sat Id   Num valid/rej code obs    Code RMS        Code WRMS       Num valid/rej doppler obs   Doppler RMS     Doppler WRMS \n");
		fprintf(fid_[index], " ------   ----------------------   ------------    ------------     -------------------------  --------------  --------------\n");

		int i;
		for (i = 0; i < MAX_N_SAT_G; i++) {
			if (	summary->sat_summary[i].code_observations.num_obs +
					summary->sat_summary[i].code_observations.num_rej_obs +
					summary->sat_summary[i].doppler_observations.num_obs +
					summary->sat_summary[i].doppler_observations.num_rej_obs > 0) {

				// Write Satellite ID
				int prn = i+1;
				if (prn < MINPRNGLO_G) {
					fprintf(fid_[index], "   G%02d   ", prn);
				} else if(prn < MAXPRNGALPRS_G) {
					fprintf(fid_[index], "   R%02d   ", prn-MINPRNGLO_G+1);
				} else {
					fprintf(fid_[index], "   E%02d   ", prn-MAXPRNGALPRS_G+1);
				}

				fprintf(fid_[index], "    %6u / %6u     %14.3f  %14.5f            %6u / %6u     %14.3f  %14.5f\n",
						summary->sat_summary[i].code_observations.num_obs,
						summary->sat_summary[i].code_observations.num_rej_obs,
						summary->sat_summary[i].code_observations.rms,
						summary->sat_summary[i].code_observations.weighted_rms,
						summary->sat_summary[i].doppler_observations.num_obs,
						summary->sat_summary[i].doppler_observations.num_rej_obs,
						summary->sat_summary[i].doppler_observations.rms,
						summary->sat_summary[i].doppler_observations.weighted_rms
				);
			}
		}

		fprintf(fid_[index], "\n ------   ----------------------   ------------    ------------     -------------------------  --------------  --------------\n");
		fprintf(fid_[index], " GLOBAL  ");
		fprintf(fid_[index], "    %6u / %6u     %14.3f  %14.5f            %6u / %6u     %14.3f  %14.5f\n",
				summary->global_summary.code_observations.num_obs,
				summary->global_summary.code_observations.num_rej_obs,
				summary->global_summary.code_observations.rms,
				summary->global_summary.code_observations.weighted_rms,
				summary->global_summary.doppler_observations.num_obs,
				summary->global_summary.doppler_observations.num_rej_obs,
				summary->global_summary.doppler_observations.rms,
				summary->global_summary.doppler_observations.weighted_rms
		);

		fprintf(fid_[index], "\n\n");
		fprintf(fid_[index], " Rejection summary               Code      Doppler\n");
		fprintf(fid_[index], " #################              ------     -------\n");
		fprintf(fid_[index], "\n");
		fprintf(fid_[index], "  Invalid PRN:                  %6u      %6u\n",
				summary->global_summary.code_observations.rej_summary.rej_invalid_prn,
				summary->global_summary.doppler_observations.rej_summary.rej_invalid_prn);
		fprintf(fid_[index], "  SV not tracked:               %6u      %6u\n",
				summary->global_summary.code_observations.rej_summary.rej_sv_not_tracked,
				summary->global_summary.doppler_observations.rej_summary.rej_sv_not_tracked);
		fprintf(fid_[index], "  SV not locked:                %6u      %6u\n",
				summary->global_summary.code_observations.rej_summary.rej_sv_not_locked,
				summary->global_summary.doppler_observations.rej_summary.rej_sv_not_locked);
		fprintf(fid_[index], "  Out of range:                 %6u      %6u\n",
				summary->global_summary.code_observations.rej_summary.rej_out_of_range,
				summary->global_summary.doppler_observations.rej_summary.rej_out_of_range);
		fprintf(fid_[index], "  SV state not available:       %6u      %6u\n",
				summary->global_summary.code_observations.rej_summary.rej_sv_state_not_available,
				summary->global_summary.doppler_observations.rej_summary.rej_sv_state_not_available);
		fprintf(fid_[index], "  Ephemeris not available       %6u      %6u\n",
				summary->global_summary.code_observations.rej_summary.rej_ephem_not_locked,
				summary->global_summary.doppler_observations.rej_summary.rej_ephem_not_locked);
		fprintf(fid_[index], "  Low SNR:                      %6u      %6u\n",
				summary->global_summary.code_observations.rej_summary.rej_low_cn0,
				summary->global_summary.doppler_observations.rej_summary.rej_low_cn0);
		fprintf(fid_[index], "  Low elevation:                %6u      %6u\n",
				summary->global_summary.code_observations.rej_summary.rej_low_elevation,
				summary->global_summary.doppler_observations.rej_summary.rej_low_elevation);
		fprintf(fid_[index], "  Chi-square test:              %6u      %6u\n",
				summary->global_summary.code_observations.rej_summary.rej_chisquare_test,
				summary->global_summary.doppler_observations.rej_summary.rej_chisquare_test);


	}
}

void close_summary_file(int index)
{
	if (fid_[index] != NULL) {
		fclose(fid_[index]);
	}
}
*/
