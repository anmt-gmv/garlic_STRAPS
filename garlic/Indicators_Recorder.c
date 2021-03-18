#include "Indicators_Recorder.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


static FILE * fid_[2] = {NULL};
static char write_header_[2] = {1};

static void write_header(int index)
{
	fprintf(fid_[index], "#                       *****************  TOTAL  ********************    **************   GPS   ************    *************  GLONASS  ***********\n");
	fprintf(fid_[index], "#    SECOND      Status    Ncode   RMScode   Ndopl   RMSdopl   Noutliers     Ncode   RMScode   Ndopl   RMSdopl      Ncode   RMScode   Ndopl   RMSdopl \n");
	fprintf(fid_[index], "# ------------   ------   -------  -------  -------  --------  ---------    -------  -------  -------  --------    -------  -------  -------  --------\n");
}

void reset_indicators(Indicators * indicators)
{
	indicators->obs_summary_total.num_code_obs = 0;
	indicators->obs_summary_total.num_doppler_obs = 0;
	indicators->obs_summary_total.num_outliers = 0;
	indicators->obs_summary_total.rms_code = 0.0;
	indicators->obs_summary_total.rms_doppler = 0.0;
	indicators->obs_summary_gps.num_code_obs = 0;
	indicators->obs_summary_gps.num_doppler_obs = 0;
	indicators->obs_summary_gps.num_outliers = 0;
	indicators->obs_summary_gps.rms_code = 0.0;
	indicators->obs_summary_gps.rms_doppler = 0.0;
	indicators->obs_summary_glonass.num_code_obs = 0;
	indicators->obs_summary_glonass.num_doppler_obs = 0;
	indicators->obs_summary_glonass.num_outliers = 0;
	indicators->obs_summary_glonass.rms_code = 0.0;
	indicators->obs_summary_glonass.rms_doppler = 0.0;
	indicators->obs_summary_galileo.num_code_obs    = 0;
	indicators->obs_summary_galileo.num_doppler_obs = 0;
	indicators->obs_summary_galileo.num_outliers    = 0;
	indicators->obs_summary_galileo.rms_code        = 0.0;
	indicators->obs_summary_galileo.rms_doppler     = 0.0;
}

void add_code_observation(sysFlag_t type, float residual, Obs_Status prange_status, Indicators * indicators)
{
	if (prange_status == OBS_VALID) {
		indicators->obs_summary_total.num_code_obs++;
		indicators->obs_summary_total.rms_code += residual * residual;
		switch (type)
		{
		case SYS_GPS:
			indicators->obs_summary_gps.num_code_obs++;
			indicators->obs_summary_gps.rms_code += residual * residual;
			break;
		case SYS_GLO:
			indicators->obs_summary_glonass.num_code_obs++;
			indicators->obs_summary_glonass.rms_code += residual * residual;
			break;
		case SYS_GAL:
			indicators->obs_summary_galileo.num_code_obs++;
			indicators->obs_summary_galileo.rms_code += residual * residual;
			break;
		default:
			// TODO include BEI
			break;
		}
	} else if (prange_status == REJ_CHISQUARE_TEST) {
		indicators->obs_summary_total.num_outliers++;
		switch (type)
		{
		case SYS_GPS:
			indicators->obs_summary_gps.num_outliers++; break;
		case SYS_GLO:
			indicators->obs_summary_glonass.num_outliers++; break;
		case SYS_GAL:
			indicators->obs_summary_glonass.num_outliers++; break;
		default:
			// TODO: Include BEI
			break;
		}
	}
}

void add_doppler_observation(sysFlag_t type, float residual, Obs_Status doppler_status, Indicators * indicators)
{
	if (doppler_status == OBS_VALID) {
		indicators->obs_summary_total.num_doppler_obs++;
		indicators->obs_summary_total.rms_doppler += residual * residual;
		switch (type)
		{
		case SYS_GPS:
			indicators->obs_summary_gps.num_doppler_obs++;
			indicators->obs_summary_gps.rms_doppler += residual * residual;
			break;
		case SYS_GLO:
			indicators->obs_summary_glonass.num_doppler_obs++;
			indicators->obs_summary_glonass.rms_doppler += residual * residual;
			break;
		case SYS_GAL:
			indicators->obs_summary_galileo.num_doppler_obs++;
			indicators->obs_summary_galileo.rms_doppler += residual * residual;
			break;
		default:
			// TODO: Include BEI
			break;
		}
	} else if (doppler_status == REJ_CHISQUARE_TEST) {
		indicators->obs_summary_total.num_outliers++;
		switch (type)
		{
		case  SYS_GPS:
			indicators->obs_summary_gps.num_outliers++; break;
		case SYS_GLO:
			indicators->obs_summary_glonass.num_outliers++; break;
		case SYS_GAL:
			indicators->obs_summary_galileo.num_outliers++; break;
		default:
			// TODO: Include BEI
			break;
		}
	}
}

void compute_indicators(Indicators * indicators)
{
	if (indicators->obs_summary_total.num_code_obs > 0) {
		indicators->obs_summary_total.rms_code = sqrt(indicators->obs_summary_total.rms_code / indicators->obs_summary_total.num_code_obs);
	}
	if (indicators->obs_summary_total.num_doppler_obs > 0) {
		indicators->obs_summary_total.rms_doppler = sqrt(indicators->obs_summary_total.rms_doppler / indicators->obs_summary_total.num_doppler_obs);
	}
	if (indicators->obs_summary_gps.num_code_obs > 0) {
		indicators->obs_summary_gps.rms_code = sqrt(indicators->obs_summary_gps.rms_code / indicators->obs_summary_gps.num_code_obs);
	}
	if (indicators->obs_summary_gps.num_doppler_obs > 0) {
		indicators->obs_summary_gps.rms_doppler = sqrt(indicators->obs_summary_gps.rms_doppler / indicators->obs_summary_gps.num_doppler_obs);
	}
	if (indicators->obs_summary_glonass.num_code_obs > 0) {
		indicators->obs_summary_glonass.rms_code = sqrt(indicators->obs_summary_glonass.rms_code / indicators->obs_summary_glonass.num_code_obs);
	}
	if (indicators->obs_summary_glonass.num_doppler_obs > 0) {
		indicators->obs_summary_glonass.rms_doppler = sqrt(indicators->obs_summary_glonass.rms_doppler / indicators->obs_summary_glonass.num_doppler_obs);
	}
	if (indicators->obs_summary_galileo.num_code_obs > 0) {
		indicators->obs_summary_galileo.rms_code = sqrt(indicators->obs_summary_galileo.rms_code / indicators->obs_summary_galileo.num_code_obs);
	}
	if (indicators->obs_summary_galileo.num_doppler_obs > 0) {
		indicators->obs_summary_galileo.rms_doppler = sqrt(indicators->obs_summary_galileo.rms_doppler / indicators->obs_summary_galileo.num_doppler_obs);
	}
}

void open_indicators_file(const char * filename, int index)
{
	fid_[index] = fopen(filename, "wt");
}

void write_indicators_epoch(double second, const Indicators * indicators, char status, int index)
{
	if (fid_[index] != NULL) {
		if (write_header_[index]) {
			write_header(index);
			write_header_[index] = 0;
		}

		fprintf(fid_[index], "  %12.4f      %1u     %7u  %7.3f  %7u  %7.3f    %6u     %7u  %7.3f  %7u  %7.3f   %7u  %7.3f  %7u  %7.3f     %7u  %7.3f  %7u  %7.3f\n",
				second, status,
				indicators->obs_summary_total.num_code_obs,
				indicators->obs_summary_total.rms_code,
				indicators->obs_summary_total.num_doppler_obs,
				indicators->obs_summary_total.rms_doppler,
				indicators->obs_summary_total.num_outliers,
				indicators->obs_summary_gps.num_code_obs,
				indicators->obs_summary_gps.rms_code,
				indicators->obs_summary_gps.num_doppler_obs,
				indicators->obs_summary_gps.rms_doppler,
				indicators->obs_summary_glonass.num_code_obs,
				indicators->obs_summary_glonass.rms_code,
				indicators->obs_summary_glonass.num_doppler_obs,
				indicators->obs_summary_glonass.rms_doppler,
				indicators->obs_summary_galileo.num_code_obs,
				indicators->obs_summary_galileo.rms_code,
				indicators->obs_summary_galileo.num_doppler_obs,
				indicators->obs_summary_galileo.rms_doppler);
	}
}

void close_indicators_file(int index)
{
	if (fid_[index] != NULL) {
		fclose(fid_[index]);
	}
}

