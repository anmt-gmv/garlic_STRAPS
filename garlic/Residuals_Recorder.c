#include "Residuals_Recorder.h"
#include "garlic_functions.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


static FILE * fid_[2] = {NULL};
static char write_header_[2] = {1};

static void write_header(int index)
{
	fprintf(fid_[index], "#    TOW      Sat  Type    Elev   CN0     Residual    Sigma    Status \n");
	fprintf(fid_[index], "# ----------- --- ------- ------ ------ ------------ -------- --------\n");
}

static void write_status(Obs_Status status, char * text)
{
	if (status == OBS_VALID) {
		sprintf(text, "OK              ");
	} else if (status == REJ_INVALID_PRN) {
		sprintf(text, "INVALID_PRN     ");
	} else if (status == REJ_SV_NOT_TRACKED) {
		sprintf(text, "SV_NOT_TRACKED  ");
	} else if (status == REJ_SV_NOT_LOCKED) {
		sprintf(text, "SV_NOT_LOCKED   ");
	} else if (status == REJ_OUT_OF_RANGE) {
		sprintf(text, "OUT_OF_RANGE    ");
	} else if (status == REJ_SV_STATE_NOT_AVAILABLE) {
		sprintf(text, "SV_STATE_NOT_AV ");
	} else if (status == REJ_EPHEM_NOT_LOCKED) {
		sprintf(text, "EPHEM_NOT_LOCKED");
	} else if (status == REJ_LOW_CN0) {
		sprintf(text, "LOW_CN0         ");
	} else if (status == REJ_LOW_ELEV) {
		sprintf(text, "LOW_ELEV        ");
	} else if (status == REJ_CHISQUARE_TEST) {
		sprintf(text, "CHISQUARE_TEST  ");
	} else if (status == REJ_FILTER_RESET) {
		sprintf(text, "FILTER_RESET    ");
	} else if (status == REJ_FILTER_INIT) {
		sprintf(text, "FILTER_INIT     ");
	} else {
		sprintf(text, "                ");
	}
}

void open_residuals_file(const char * filename, int index)
{
	fid_[index] = fopen(filename, "wt");
}

void write_residual(double tow, obsdata_t * obs, int index)
{
	if (fid_[index] != NULL) {
		if (write_header_[index]) {
			write_header(index);
			write_header_[index] = 0;
		}

		char sat_id[4];
		sigFlag_t sys = signal2system(obs->sigFlag);
		switch (sys)
		{
		case SYS_GPS:
			sprintf(sat_id, "G%02d", garlic2sysIndex(obs->PRN-1,sys)+1); break;
		case SYS_GLO:
			sprintf(sat_id, "R%02d", garlic2sysIndex(obs->PRN-1,sys)+1); break;
		case SYS_GAL:
			sprintf(sat_id, "E%02d", garlic2sysIndex(obs->PRN-1,sys)+1); break;
		case SYS_BEI:
			sprintf(sat_id, "C%02d", garlic2sysIndex(obs->PRN-1,sys)+1); break;
		default:
			sprintf(sat_id, "-%02d", garlic2sysIndex(obs->PRN-1,sys)+1); break;
		}

		float elevation = obs->elev;
		if (elevation < 0.0) {
			elevation = 0.0;
		}

		// CODE RESIDUAL
		char status[17];
		write_status(obs->prange_status, status);

		fprintf(fid_[index], "  %11.4f %3s %7s  %4.1f   %4.1f  %12.3f %8.4f  %16s\n",
				tow, sat_id, "CODE   ", obs->elev, obs->S1, obs->prange_residual, obs->prange_sigma, status);

		// DOPPLER RESIDUAL
		write_status(obs->doppler_status, status);

		fprintf(fid_[index], "  %11.4f %3s %7s  %4.1f   %4.1f  %12.3f %8.4f  %16s\n",
				tow, sat_id, "DOPPLER", obs->elev, obs->S1, obs->doppler_residual, obs->doppler_sigma, status);
	}
}

void close_residuals_file(int index)
{
	if (fid_[index] != NULL) {
		fclose(fid_[index]);
	}
}

