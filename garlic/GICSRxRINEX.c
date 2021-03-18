#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "GICSRxRINEX.h"

#include "calendar.h"
#include "atmsphcorr.h"

#include "garlic_functions.h"

#define MIN_EPHGPS_DISTANCE 3600
#define MIN_EPHGLO_DISTANCE 900
#define MIN_EPHGAL_DISTANCE 3600
#define MIN_EPHBEI_DISTANCE 3600

#define MAXLENGTHRV3 800

#define USE_GALPRS 0

#if USE_GALPRS == 0
#define C1X "C1B"
#define D1X "D1B"
#define L1X "L1B"
#define S1X "S1B"
#else
#define C1X "C1A"
#define D1X "D1A"
#define L1X "L1A"
#define S1X "S1A"
#endif

static rinex_header_t rinex_info;

static char      estimate_doppler = 0;
static obsdata_t last_obs[MAX_CHANNELS_G];

/* Copy the information gathered from the different RINEX headers (*.obs, *.nav) */
void get_rinex_info(rinex_header_t *info)
{
	memcpy(info, &rinex_info, sizeof(rinex_header_t));
}

/* Replace the 'D' exponential indicator by C standard 'E' to use sscanf */
static void replace_rinex_string(char *line)
{
	int k;
	for(k = 0; k < strlen(line); k++)
	{
		if(line[k] == 'D')
		{
			line[k] = 'E';
		}
	}
}

/* Calculate the time difference in seconds between two GPS times, expressed in (week,tow) */
static double calculate_time_difference(int WN_1, double ToW_1, int WN_2, double ToW_2)
{
	double gps_time_1 = WN_1*SECONDSONWEEK + ToW_1;
	double gps_time_2 = WN_2*SECONDSONWEEK + ToW_2;

	return (gps_time_1 - gps_time_2);
}

/* This function computes the leap second value for a given week number and TOW */
static int get_leap_second(long week, double time)
{
	long tmp_week;
	double tmp_time;
	int k, leap_second = 0;

	for(k = 0; k < MAXLEAPS; k++)
	{
		CalToGPS_G(LEAPS[k][0], LEAPS[k][1], LEAPS[k][2], LEAPS[k][3], LEAPS[k][4], LEAPS[k][5], &tmp_week, &tmp_time);
		if (calculate_time_difference(week, time, tmp_week, tmp_time) >= 0)
		{
			leap_second = LEAPS[k][6];
			break;
		}
	}
	return leap_second;
}

static int getConstellationFromRINEX(char sys_code)
{
	switch (sys_code)
	{
	case 'G': return SYS_GPS;
	case 'E': return SYS_GAL;
	case 'R': return SYS_GLO;
	case 'C': return SYS_BEI;
	default:  return -1;
	}
}

static int getServiceFromRINEX (char sys_code, char *obs_code)
{
	int sig = -1;

	switch (sys_code)
	{
	case 'G':
		if      (strncmp(obs_code, "1C", 2) == 0) { sig = SIG_GPSL1;  }
		else if (strncmp(obs_code, "5I", 2) == 0) { sig = SIG_GPSL5I; }
		else if (strncmp(obs_code, "5Q", 2) == 0) { sig = SIG_GPSL5Q; }
		break;

	case 'E':
		if      (strncmp(obs_code, "1A", 2) == 0) { sig = SIG_GALE1AD; }
		else if (strncmp(obs_code, "1T", 2) == 0) { sig = SIG_GALE1AP; }
		else if (strncmp(obs_code, "1B", 2) == 0) { sig = SIG_GALE1B;  }
		else if (strncmp(obs_code, "1C", 2) == 0) { sig = SIG_GALE1C;  }
		else if (strncmp(obs_code, "1X", 2) == 0) { sig = SIG_GALE1B;  }
		else if (strncmp(obs_code, "5I", 2) == 0) { sig = SIG_GALE5aI; }
		else if (strncmp(obs_code, "5Q", 2) == 0) { sig = SIG_GALE5aQ; }
		else if (strncmp(obs_code, "7I", 2) == 0) { sig = SIG_GALE5bI; }
		else if (strncmp(obs_code, "7Q", 2) == 0) { sig = SIG_GALE5bQ; }
		break;

	case 'R':
		if      (strncmp(obs_code, "1C", 2) == 0) { sig = SIG_GLOL1; }
//		if      (strncmp(obs_code, "1P", 2) == 0) { sig = SIG_GLOL1; }
		// To support RINEX v2
		else if (obs_code[0] == '1') { sig = SIG_GLOL1; }
		break;

	case 'C':
		if      (strncmp(obs_code, "2I", 2) == 0) { sig = SIG_BEIB1; }
		if      (strncmp(obs_code, "1I", 2) == 0) { sig = SIG_BEIB1; }
		break;
	}
	return sig;
}

/* Update the observation data structure with the field read from RINEX file */
static void process_obs_field(char *field, char obs_type, obsdata_t *obsdata)
{
	// These SV types is not being accounted in the current version.
	if(obsdata->sigFlag < 0) { return; }

	switch (obs_type)
	{
	case 'C': // Pseudorange measurement
		if(obsdata->lock == 1 && (obsdata->C1 = atof(field)) != 0)
		{
			obsdata->prange_status = OBS_VALID;
		}
		break;

	case 'L': // Phase measurement
		if(obsdata->lock == 1 && (obsdata->L1 = atof(field)) != 0)
		{
			obsdata->cphase_status = OBS_VALID;

			// TODO: think in a better method to estimate doppler from the carrier phase measurements.
			if(estimate_doppler == 1)
			{
				int j;
				for(j = 0; j < MAX_CHANNELS_G; j++)
				{
					if(last_obs[j].PRN == obsdata->PRN && last_obs[j].cphase_status == OBS_VALID && last_obs[j].sigFlag == obsdata->sigFlag)
					{
						obsdata->D1  = (obsdata->L1 - last_obs[j].L1) * obsdata->lambda;

						if(last_obs[j].doppler_status == OBS_VALID && fabs(obsdata->D1 - last_obs[j].D1) < 0.5)
						{
							obsdata->D1 = 2*obsdata->D1 - last_obs[j].D1;
						}
						break;
					}
				}
			}
		}
		break;

	case 'D': // Doppler measurement
		if(obsdata->lock == 1 && (obsdata->D1 = -atof(field)) != 0)
		{
			obsdata->D1 *= obsdata->lambda;
			obsdata->doppler_status = OBS_VALID;
		}
		break;

	case 'S':
		if(obsdata->lock == 1)
		{
			obsdata->S1 = atof(field);
		}
	}
}

static void add_new_observation(gnssdata_t *gnssdata, sigFlag_t sig, int sat_no)
{

	int ind = gnssdata->noOfChannelsAv;
	int sys = signal2system(sig);
	int prnIndex = sys2garlicIndex(sat_no - 1, signal2system(sig));
	char  vflg = 0;

	switch (sys)
	{
	case SYS_GPS: vflg = gnssdata->EPH[prnIndex].GPS.vflg; break;
	case SYS_GAL: vflg = gnssdata->EPH[prnIndex].GAL.vflg; break;
	case SYS_GLO: vflg = gnssdata->EPH[prnIndex].GLO.vflg; break;
	case SYS_BEI: vflg = gnssdata->EPH[prnIndex].BEI.vflg; break;
	}

	gnssdata->OBS[ind].channel   = 0;
	gnssdata->OBS[ind].lock      = (vflg || gnssdata->SP3_ON);
	gnssdata->OBS[ind].sigFlag   = sig;
	gnssdata->OBS[ind].iono_free = 0;

	gnssdata->OBS[ind].S1 = 40;

	gnssdata->OBS[ind].prange_status  = REJ_SV_NOT_LOCKED;
	gnssdata->OBS[ind].cphase_status  = REJ_SV_NOT_LOCKED;
	gnssdata->OBS[ind].doppler_status = REJ_SV_NOT_LOCKED;

	gnssdata->OBS[ind].PRN  = prnIndex + 1;
	gnssdata->SAT_ID[ind]   = prnIndex + 1;

	if (sys == SYS_GLO)
	{
		gnssdata->OBS[ind].channel = (gnssdata->EPH[prnIndex].GLO.freqslot);
		gnssdata->OBS[ind].lambda = get_wavelength(sig, 1, gnssdata->OBS[ind].channel+7);
	} else {
		gnssdata->OBS[ind].lambda = get_wavelength(sig, 1, 0);
	}

	gnssdata->noOfChannelsAv++;
}

void read_obs_header_RINEX(FILE *fp_obs, gnssdata_t *gnssdata, lsq_state_t *lsq_state)
{
	char current_line[100];
	int updated_ls = 0;

	rinex_info.epochs_interval         = 1.0;

	while(1)
	{
		fgets(current_line, 100, fp_obs);

		if(strncmp(&current_line[60], "RINEX VERSION / TYPE", 20) == 0)
		{
			if(current_line[5] == '2')
			{
				rinex_info.obs_version = RINEX_v2;
			}
			else if(current_line[5] == '3')
			{
				rinex_info.obs_version = RINEX_v3;
			}
		}

		if(strncmp(&current_line[60], "APPROX POSITION XYZ", 19) == 0)
		{
			sscanf(current_line, "%lf %lf %lf %*s", &rinex_info.XYZ_antenna[0],
					&rinex_info.XYZ_antenna[1],  &rinex_info.XYZ_antenna[2]);
		}

		if(strncmp(&current_line[60], "ANTENNA: DELTA H/E/N", 20) == 0)
		{
			sscanf(current_line, "%lf %lf %lf %*s",  &rinex_info.delta_antenna[0],
					&rinex_info.delta_antenna[1], &rinex_info.delta_antenna[2]);
		}

		if(strncmp(&current_line[60], "INTERVAL", 8) == 0)
		{
			sscanf(current_line, "%lf %*s", &rinex_info.epochs_interval);
		}

		// Observation type in RINEX v2
		if(strncmp(&current_line[60], "# / TYPES OF OBSERV", 19) == 0)
		{

			sscanf(current_line, "%d %c%c %c%c %c%c %c%c %c%c %c%c %c%c %c%c %c%c",
					&rinex_info.number_of_observations[0],
					&rinex_info.observation_list[0][0][0], &rinex_info.observation_list[0][0][1],
					&rinex_info.observation_list[0][1][0], &rinex_info.observation_list[0][1][1],
					&rinex_info.observation_list[0][2][0], &rinex_info.observation_list[0][2][1],
					&rinex_info.observation_list[0][3][0], &rinex_info.observation_list[0][3][1],
					&rinex_info.observation_list[0][4][0], &rinex_info.observation_list[0][4][1],
					&rinex_info.observation_list[0][5][0], &rinex_info.observation_list[0][5][1],
					&rinex_info.observation_list[0][6][0], &rinex_info.observation_list[0][6][1],
					&rinex_info.observation_list[0][7][0], &rinex_info.observation_list[0][7][1],
					&rinex_info.observation_list[0][8][0], &rinex_info.observation_list[0][8][1]);

			if(rinex_info.number_of_observations[0] > 9)
			{
				fgets(current_line, 100, fp_obs);

				sscanf(current_line, "          %c%c %c%c %c%c %c%c %c%c %c%c %c%c %c%c %c%c",
						&rinex_info.observation_list[0][ 9][0], &rinex_info.observation_list[0][ 9][1],
						&rinex_info.observation_list[0][10][0], &rinex_info.observation_list[0][10][1],
						&rinex_info.observation_list[0][11][0], &rinex_info.observation_list[0][11][1],
						&rinex_info.observation_list[0][12][0], &rinex_info.observation_list[0][12][1],
						&rinex_info.observation_list[0][13][0], &rinex_info.observation_list[0][13][1],
						&rinex_info.observation_list[0][14][0], &rinex_info.observation_list[0][14][1],
						&rinex_info.observation_list[0][15][0], &rinex_info.observation_list[0][15][1],
						&rinex_info.observation_list[0][16][0], &rinex_info.observation_list[0][16][1],
						&rinex_info.observation_list[0][17][0], &rinex_info.observation_list[0][17][1]);
			}

			if(rinex_info.number_of_observations[0] > 18)
			{
				fgets(current_line, 100, fp_obs);

				sscanf(current_line, "          %c%c %c%c %c%c %c%c %c%c %c%c %c%c %c%c",
						&rinex_info.observation_list[0][18][0], &rinex_info.observation_list[0][18][1],
						&rinex_info.observation_list[0][19][0], &rinex_info.observation_list[0][19][1],
						&rinex_info.observation_list[0][20][0], &rinex_info.observation_list[0][20][1],
						&rinex_info.observation_list[0][21][0], &rinex_info.observation_list[0][21][1],
						&rinex_info.observation_list[0][22][0], &rinex_info.observation_list[0][22][1],
						&rinex_info.observation_list[0][23][0], &rinex_info.observation_list[0][23][1],
						&rinex_info.observation_list[0][24][0], &rinex_info.observation_list[0][24][1],
						&rinex_info.observation_list[0][25][0], &rinex_info.observation_list[0][25][1]);
			}

			int k;
			for(k = 1; k < NUM_SYS; k++)
			{
				rinex_info.number_of_observations[k] = rinex_info.number_of_observations[0];
				memcpy(rinex_info.observation_list[k], rinex_info.observation_list[0], sizeof(rinex_info.observation_list[0]));
			}

			if(rinex_info.obs_version == RINEX_v3)
			{
				printf("Warning: header label \"# / TYPES OF OBSERV\" not supported in RINEX v3. Trying to set RINEX v2...\n");
				rinex_info.obs_version = RINEX_v2;
			}

			char doppler_found = 0, phase_found = 0;
			for(k = 0; k < rinex_info.number_of_observations[0]; k++)
			{
				phase_found   |= (rinex_info.observation_list[0][k][0] == 'L' && rinex_info.observation_list[0][k][1] == '1');
				doppler_found |= (rinex_info.observation_list[0][k][0] == 'D' && rinex_info.observation_list[0][k][1] == '1');
			}
			estimate_doppler = (phase_found && !doppler_found);
		}

		// Observation type in RINEX v3
		if(strncmp(&current_line[60], "SYS / # / OBS TYPES", 19) == 0)
		{
			char sys_code;
			int total_num_obs = 0, obs = 0, idx = 7;
			int sys;

			sscanf(current_line,"%c %d",&sys_code, &total_num_obs);

			sys = getConstellationFromRINEX(sys_code);

			if(sys < 0) { continue; }

			rinex_info.number_of_observations[sys] = total_num_obs;

			for (obs = 0; obs < total_num_obs; obs++)
			{
				if (obs != 0 && obs % 13 == 0) {
					fgets(current_line, 100, fp_obs);
					idx = 7;
				}
				sscanf(current_line+idx,"%3s", rinex_info.observation_list[sys][obs]);
				idx += 4;
			}
		}

		if(strncmp(&current_line[60], "LEAP SECONDS", 12) == 0)
		{
			sscanf(current_line, "%d %*s", &rinex_info.leap_seconds);
			gnssdata->iono_utc.gps.dtls  = rinex_info.leap_seconds;
			gnssdata->iono_utc.gal.dtlsf = rinex_info.leap_seconds;
			updated_ls = 1;
		}

		if(strncmp(&current_line[60], "TIME OF FIRST OBS", 17) == 0)
		{
			int ye, mo, da, ho, mi;
			double sec;
			sscanf(current_line, "%d %d %d %d %d %lf %*s", &ye, &mo, &da, &ho, &mi, &sec);

			CalToGPS_G(ye, mo, da, ho, mi, sec, &rinex_info.first_week, &rinex_info.first_tow);
		}

		if(strncmp(&current_line[60], "END OF HEADER", 13) == 0)
		{
			break;
		}
	}

	gnssdata->week = rinex_info.first_week;
	gnssdata->tow  = rinex_info.first_tow - rinex_info.epochs_interval;

	if (!updated_ls && gnssdata->week > 0) {
		int lp = get_leap_second(gnssdata->week,gnssdata->tow);
		gnssdata->iono_utc.gps.dtls =  lp;
		gnssdata->iono_utc.gal.dtlsf = lp;
	}


	memset(lsq_state->pos, 0, sizeof(lsq_state->pos));
	lsq_state->pos[0] = rinex_info.XYZ_antenna[0];
	lsq_state->pos[1] = rinex_info.XYZ_antenna[1];
	lsq_state->pos[2] = rinex_info.XYZ_antenna[2];

	memset(lsq_state->vel, 0, sizeof(lsq_state->vel));
}

static void read_nav_header_RINEX_v2(FILE *fp_nav, gnssdata_t *gnssdata)
{
	char current_line[100];

	char iono_model_ready = -1;

	while(1)
	{
		fgets(current_line, 100, fp_nav);

#if CALC_ATMSPHCORR == 1
		if(strncmp(&current_line[60], "ION ALPHA", 8) == 0)
		{
			double ion_alpha[4];
			replace_rinex_string(current_line);
			sscanf(current_line," %le %le %le %le %*s",	&ion_alpha[0], &ion_alpha[1], &ion_alpha[2], &ion_alpha[3]);

			gnssdata->iono_utc.gps.alpha[0] = ion_alpha[0];
			gnssdata->iono_utc.gps.alpha[1] = ion_alpha[1];
			gnssdata->iono_utc.gps.alpha[2] = ion_alpha[2];
			gnssdata->iono_utc.gps.alpha[3] = ion_alpha[3];

			iono_model_ready++;
		}
		if(strncmp(&current_line[60], "ION BETA", 7) == 0)
		{
			double ion_beta[4];
			replace_rinex_string(current_line);
			sscanf(current_line," %le %le %le %le %*s",	&ion_beta[0], &ion_beta[1], &ion_beta[2], &ion_beta[3]);

			gnssdata->iono_utc.gps.beta[0] = ion_beta[0];
			gnssdata->iono_utc.gps.beta[1] = ion_beta[1];
			gnssdata->iono_utc.gps.beta[2] = ion_beta[2];
			gnssdata->iono_utc.gps.beta[3] = ion_beta[3];

			iono_model_ready++;
		}
		if(strncmp(&current_line[0], "GAL", 3) == 0)
		{
			double iono_gal[4];
			replace_rinex_string(current_line);
			sscanf(current_line,"%*s %le %le %le %le %*s", &iono_gal[0], &iono_gal[1], &iono_gal[2], &iono_gal[3]);

			gnssdata->iono_utc.gal.ai[0]  = iono_gal[0];
			gnssdata->iono_utc.gal.ai[1]  = iono_gal[1];
			gnssdata->iono_utc.gal.ai[2]  = iono_gal[2];

			gnssdata->iono_utc.gal.vflg[0] = 1;
			gnssdata->iono_utc.gal.vflg[1] = 1;
		}
#endif
		if(strncmp(&current_line[60], "DELTA-UTC: A0,A1,T,W", 20) == 0)
		{
			double A0, A1;
			int dtot;

			replace_rinex_string(current_line);
			sscanf(current_line, "%le%le %d %d %*s", &A0, &A1, &dtot, &gnssdata->iono_utc.gps.WNt);

			gnssdata->iono_utc.gps.A0  = A0;
			gnssdata->iono_utc.gps.A1  = A1;
			gnssdata->iono_utc.gps.tot = dtot;

			gnssdata->iono_utc.gps.GPSweek = gnssdata->iono_utc.gps.WNt;
		}
		if(strncmp(&current_line[60], "LEAP SECONDS", 12) == 0)
		{
			if(gnssdata->iono_utc.gps.dtls == 0)
			{
				sscanf(current_line, "%d %*s", &rinex_info.leap_seconds);
				gnssdata->iono_utc.gps.dtls  = rinex_info.leap_seconds;
				gnssdata->iono_utc.gal.dtlsf = rinex_info.leap_seconds;
			}
			else
			{
				rinex_info.leap_seconds = gnssdata->iono_utc.gps.dtls;
			}
		}
		if(strncmp(&current_line[60], "CORR TO SYSTEM TIME", 19) == 0)
		{
			replace_rinex_string(current_line);

			double tau_c;
			sscanf(current_line, "%*d %*d %*d %lf", &tau_c);

			int k;
			for(k = 0; k < MAXNUMGLO_G; k++)
			{
				gnssdata->EPH[sys2garlicIndex(k,SYS_GLO)].GLO.tau_c = tau_c;
			}
		}
		if(strncmp(&current_line[60], "PGM / RUN BY / DATE", 19) == 0)
		{
			// Ignore line.
		}
		if(strncmp(&current_line[60], "COMMENT", 7) == 0)
		{
			// Ignore line.
		}
		if(strncmp(&current_line[60], "END OF HEADER", 13) == 0)
		{
			break;
		}
	}
	gnssdata->iono_utc.gps.vflg |= (iono_model_ready >= 1);
}

static void read_nav_header_RINEX_v3(FILE *fp_nav, gnssdata_t *gnssdata)
{
	char current_line[100];

	char iono_model_ready = -1;

	while(1)
	{
		fgets(current_line, 100, fp_nav);

#if CALC_ATMSPHCORR == 1
		if(strncmp(&current_line[60], "IONOSPHERIC CORR", 16) == 0)
		{
			if(strncmp(&current_line[0], "GPSA", 4) == 0)
			{
				double ion_alpha[4];
				replace_rinex_string(current_line);
				sscanf(current_line,"%*s %le %le %le %le %*s", &ion_alpha[0], &ion_alpha[1], &ion_alpha[2], &ion_alpha[3]);

				gnssdata->iono_utc.gps.alpha[0] = ion_alpha[0];
				gnssdata->iono_utc.gps.alpha[1] = ion_alpha[1];
				gnssdata->iono_utc.gps.alpha[2] = ion_alpha[2];
				gnssdata->iono_utc.gps.alpha[3] = ion_alpha[3];

				iono_model_ready++;
			}
			if(strncmp(&current_line[0], "GPSB", 4) == 0)
			{
				double ion_beta[4];
				replace_rinex_string(current_line);
				sscanf(current_line,"%*s %le %le %le %le %*s", &ion_beta[0], &ion_beta[1], &ion_beta[2], &ion_beta[3]);

				gnssdata->iono_utc.gps.beta[0] = ion_beta[0];
				gnssdata->iono_utc.gps.beta[1] = ion_beta[1];
				gnssdata->iono_utc.gps.beta[2] = ion_beta[2];
				gnssdata->iono_utc.gps.beta[3] = ion_beta[3];

				iono_model_ready++;
			}
			if(strncmp(&current_line[0], "GAL", 3) == 0)
			{
				double iono_gal[4];
				replace_rinex_string(current_line);
				sscanf(current_line,"%*s %le %le %le %le %*s", &iono_gal[0], &iono_gal[1], &iono_gal[2], &iono_gal[3]);

				gnssdata->iono_utc.gal.ai[0] = iono_gal[0];
				gnssdata->iono_utc.gal.ai[1] = iono_gal[1];
				gnssdata->iono_utc.gal.ai[2] = iono_gal[2];

				gnssdata->iono_utc.gal.vflg[0] = 1;
				gnssdata->iono_utc.gal.vflg[1] = 1;
			}
		}
#endif

		if(strncmp(&current_line[60], "LEAP SECONDS", 12) == 0)
		{
			if(gnssdata->iono_utc.gps.dtls == 0)
			{
				sscanf(current_line, "%d %*s", &rinex_info.leap_seconds);
				gnssdata->iono_utc.gps.dtls  = rinex_info.leap_seconds;
				gnssdata->iono_utc.gal.dtlsf = rinex_info.leap_seconds;
			}
		}
		if((strncmp(&current_line[60], "CORR TO SYSTEM TIME", 19) == 0) ||
		   (strncmp(&current_line[60],    "TIME SYSTEM CORR", 16) == 0))
		{
			double A0, A1;
			int dtot;

			replace_rinex_string(current_line);

			if(strncmp(&current_line[0], "GPUT", 4) == 0)
			{
				sscanf(current_line, "%le%le %d %d %*s", &A0, &A1, &dtot, &gnssdata->iono_utc.gps.WNt);

				gnssdata->iono_utc.gps.A0  = A0;
				gnssdata->iono_utc.gps.A1  = A1;
				gnssdata->iono_utc.gps.tot = dtot;

				gnssdata->iono_utc.gps.GPSweek = gnssdata->iono_utc.gps.WNt;
			}
			if(strncmp(&current_line[0], "GAUT", 4) == 0)
			{
				sscanf(current_line, "GAUT %le%le %d %d %*s", &A0, &A1, &dtot, &gnssdata->iono_utc.gal.wnt);

				gnssdata->iono_utc.gal.A0  = A0;
				gnssdata->iono_utc.gal.A1  = A1;
				gnssdata->iono_utc.gal.tot = dtot;
			}
			if(strncmp(&current_line[0], "GPGA", 4) == 0)
			{
				sscanf(current_line, "GPGA %le%le %d %d %*s", &A0, &A1, &dtot, &gnssdata->iono_utc.gal.wng);

				gnssdata->iono_utc.gal.A0g = A0;
				gnssdata->iono_utc.gal.A1g = A1;
				gnssdata->iono_utc.gal.tog = dtot;
			}
		}
		if(strncmp(&current_line[60], "PGM / RUN BY / DATE", 19) == 0)
		{
			// Ignore line.
		}
		if(strncmp(&current_line[60], "COMMENT", 7) == 0)
		{
			// Ignore line.
		}
		if(strncmp(&current_line[60], "END OF HEADER", 13) == 0)
		{
			break;
		}
	}
	gnssdata->iono_utc.gps.vflg |= (iono_model_ready >= 1);

#if CALC_ATMSPHCORR == 1
	if(gnssdata->iono_utc.gal.vflg)
	{
		initialize_NeQuick_iono_model("nQfiles\\", &gnssdata->iono_utc.gal);
	}
#endif
}

char read_nav_header_RINEX(FILE *fp_nav, gnssdata_t *gnssdata)
{
	char constellation_type = -1;

	char current_line[100];

	fgets(current_line, 100, fp_nav);

	if(strncmp(&current_line[60], "RINEX VERSION / TYPE", 20) == 0)
	{
		if(current_line[5] == '2')
		{
			if(strncmp(&current_line[20], "NAVIGATION DATA", 15) == 0 ||
			   strncmp(&current_line[20], "N: GPS NAV DATA", 15) == 0)
			{
				rinex_info.navgps_version = RINEX_v2;
				constellation_type = 0;
			}
			if(strncmp(&current_line[20], "GLONASS NAV DATA", 16) == 0)
			{
				rinex_info.navglo_version = RINEX_v2;
				constellation_type = 1;
			}
			if(strncmp(&current_line[20], "E: GALILEO NAV DATA", 19) == 0)
			{
				rinex_info.navgal_version = RINEX_v2;
				constellation_type = 2;
			}
			read_nav_header_RINEX_v2(fp_nav, gnssdata);
		}

		if(current_line[5] == '3')
		{
			if(current_line[40] == 'G')
			{
				rinex_info.navgps_version   = RINEX_v3;
				rinex_info.mixed_navigation = 0;
				constellation_type = 0;
			}
			if(current_line[40] == 'R')
			{
				rinex_info.navglo_version   = RINEX_v3;
				rinex_info.mixed_navigation = 0;
				constellation_type = 1;
			}
			if(current_line[40] == 'E')
			{
				rinex_info.navgal_version   = RINEX_v3;
				rinex_info.mixed_navigation = 0;
				constellation_type = 2;
			}
			if(current_line[40] == 'C')
			{
				rinex_info.navbei_version   = RINEX_v3;
				rinex_info.mixed_navigation = 0;
				constellation_type = 3;
			}
			if(current_line[40] == 'M')
			{
				rinex_info.navgps_version   = RINEX_v3;
				rinex_info.navglo_version   = RINEX_v3;
				rinex_info.navgal_version   = RINEX_v3;
				rinex_info.navbei_version   = RINEX_v3;
				rinex_info.mixed_navigation = 1;
				constellation_type = 4;
			}
			read_nav_header_RINEX_v3(fp_nav, gnssdata);
		}
	}
	else
	{
		printf("Wrong format in input RINEX navigation files\n");
	}
	return constellation_type;
}

static void parse_nav_gps_RINEX(FILE *fp_nav, long *last_valid_cursor, int week, double gps_time, int prn, ephgps_t *ephgps)
{
	int ye, mo, da, ho, mi;
	double sec;

	double field_1, field_2, field_3, field_4;

	char current_line[100];

	ephgps_t tmpEph;

	double min_tow_distance = calculate_time_difference(week, gps_time, ephgps->nWeekNumber, ephgps->dtoe);
	while(!feof(fp_nav))
	{
		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);
		if(rinex_info.navgps_version == RINEX_v3)
		{
			char sys_type;
			sscanf(current_line, "%c%d %d %d %d %d %d %lf%le%le%le", &sys_type, &tmpEph.prn, &ye, &mo, &da, &ho, &mi, &sec, &field_1, &field_2, &field_3);
			if(sys_type != 'G')
			{
				int k;
				for(k = 0; k < ((sys_type == 'R' || sys_type == 'S') ? 3 : 7); k++)
				{
					fgets(current_line, 100, fp_nav);
				}
				continue;
			}
		}
		else
		{
			sscanf(current_line, "%d %d %d %d %d %d %lf%le%le%le", &tmpEph.prn, &ye, &mo, &da, &ho, &mi, &sec, &field_1, &field_2, &field_3);
			ye += 2000;
		}

		tmpEph.dAfo = field_1;
		tmpEph.dAf1 = field_2;
		tmpEph.dAf2 = field_3;

		CalToGPS_G(ye, mo, da, ho, mi, sec, (long int *)&tmpEph.nWeekCLK, &tmpEph.dTOC);

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.nIODE = field_1;
		tmpEph.dCrs  = field_2;
		tmpEph.dDn   = field_3;
		tmpEph.dMo   = field_4;

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.dCuc   = field_1;
		tmpEph.dEcc   = field_2;
		tmpEph.dCus   = field_3;
		tmpEph.dAroot = field_4;

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.dtoe    = field_1;
		tmpEph.dCic    = field_2;
		tmpEph.dOMEGAo = field_3;
		tmpEph.dCis    = field_4;

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.d_io      = field_1;
		tmpEph.dCrc      = field_2;
		tmpEph.d_w       = field_3;
		tmpEph.dOMEGADOT = field_4;

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.dIDOT       = field_1;
		tmpEph.CodeOnL2    = field_2;
		tmpEph.nWeekNumber = field_3;
		tmpEph.L2PFlag     = field_4;

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.accuracy = field_1;
		tmpEph.nHealth  = field_2;
		tmpEph.dTGD     = field_3;
		tmpEph.nIODC    = field_4;

		fgets(current_line, 100, fp_nav);

		tmpEph.vflg = (tmpEph.nHealth == 0);

		tmpEph.updateTime = 0;

		double time_difference = calculate_time_difference(tmpEph.nWeekNumber, tmpEph.dtoe, week, gps_time);
		if(!feof(fp_nav) && tmpEph.prn == prn)
		{
			if(fabs(time_difference) <= min_tow_distance)
			{
				min_tow_distance = fabs(time_difference);
				memcpy(ephgps, &tmpEph, sizeof(ephgps_t));

				if(rinex_info.mixed_navigation == 1)
				{
					last_valid_cursor[sys2garlicIndex(prn-1, SYS_GPS)] = ftell(fp_nav);
				}
				else
				{
					last_valid_cursor[prn-1] = ftell(fp_nav);
				}

				if(min_tow_distance < MIN_EPHGPS_DISTANCE)
				{
					break;
				}
			}
			else if(tmpEph.dtoe > gps_time)
			{
				break;
			}
		}
	}
}

static void parse_nav_glo_RINEX(FILE *fp_nav, long *last_valid_cursor, int week, double gps_time, int prn, ephglo_t *ephglo, ionoutc_t *iono_utc)
{
	int   ye, mo, da, ho, mi;
	double sec;

	double field_1, field_2, field_3, field_4;

	char current_line[100];

	ephglo_t tmpEph;

	double min_tow_distance = calculate_time_difference(week, gps_time, ephglo->nWeekNumber, ephglo->dtoe);
	while(!feof(fp_nav))
	{
		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		if(rinex_info.navglo_version == RINEX_v3)
		{
			char sys_type;
			sscanf(current_line, "%c%d %d %d %d %d %d %lf%le%le%le", &sys_type, &tmpEph.prn, &ye, &mo, &da, &ho, &mi, &sec, &field_1, &field_2, &field_3);
			if(sys_type != 'R')
			{
				int k;
				for(k = 0; k < (sys_type == 'S' ? 3 : 7); k++)
				{
					fgets(current_line, 100, fp_nav);
				}
				continue;
			}
		}
		else
		{
			sscanf(current_line, "%d %d %d %d %d %d %lf%le%le%le", &tmpEph.prn, &ye, &mo, &da, &ho, &mi, &sec, &field_1, &field_2, &field_3);
			ye += 2000;
		}

		tmpEph.SV_timedev = -field_1;
		tmpEph.SV_freqdev = field_2;
		tmpEph.tk         = field_3;

		CalToGPS_G(ye, mo, da, ho, mi, sec, (long int *)&tmpEph.nWeekNumber, &tmpEph.dtoe);
		tmpEph.dtoe = tmpEph.dtoe + iono_utc->gps.dtls;

		tmpEph.tb = fmod(tmpEph.dtoe, DAYSECS);

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.X_tb     = field_1;
		tmpEph.Xdot_tb  = field_2;
		tmpEph.Xdot2_tb = field_3;
		tmpEph.nHealth  = field_4;

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.Y_tb     = field_1;
		tmpEph.Ydot_tb  = field_2;
		tmpEph.Ydot2_tb = field_3;
		tmpEph.freqslot = field_4;

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.Z_tb     = field_1;
		tmpEph.Zdot_tb  = field_2;
		tmpEph.Zdot2_tb = field_3;
		tmpEph.ageDays  = field_4;

		tmpEph.timeMark  = 0;
		tmpEph.dTGD      = 0;
		tmpEph.vflg = (tmpEph.nHealth == 0);

		// Unused parameters.
		tmpEph.P  = 0;
		tmpEph.P1 = 0;
		tmpEph.P2 = 0;
		tmpEph.P3 = 0;
		tmpEph.P4 = 0;

		double time_difference = calculate_time_difference(tmpEph.nWeekNumber, tmpEph.dtoe, week, gps_time);
		if(!feof(fp_nav) && tmpEph.prn == prn)
		{
			if(fabs(time_difference) <= min_tow_distance)
			{
				double tmp_tau_c = ephglo->tau_c;

				min_tow_distance = fabs(time_difference);
				memcpy(ephglo, &tmpEph, sizeof(ephglo_t));
				ephglo->tau_c = tmp_tau_c;

				if(rinex_info.mixed_navigation == 1)
				{
					last_valid_cursor[sys2garlicIndex(prn-1, SYS_GLO)] = ftell(fp_nav);
				}
				else
				{
					last_valid_cursor[prn-1] = ftell(fp_nav);
				}

				if(min_tow_distance < MIN_EPHGLO_DISTANCE)
				{
					break;
				}
			}
			else if(tmpEph.dtoe > gps_time && tmpEph.dtoe < 604800)
			{
				break;
			}
		}
	}
}

static void parse_nav_gal_RINEX(FILE *fp_nav, long *last_valid_cursor, int week, double gps_time, int prn, ephgal_t *ephgal)
{
	int ye, mo, da, ho, mi;
	double sec;

	double field_1, field_2, field_3, field_4;

	char current_line[100];

	ephgal_t tmpEph;

	double min_tow_distance = calculate_time_difference(week, gps_time, ephgal->nWeekNumber, ephgal->dtoe);
	while(!feof(fp_nav))
	{
		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);
		if(rinex_info.navgal_version == RINEX_v3)
		{
			char sys_type;
			sscanf(current_line, "%c%d %d %d %d %d %d %lf%le%le%le", &sys_type, &tmpEph.prn, &ye, &mo, &da, &ho, &mi, &sec, &field_1, &field_2, &field_3);
			if(sys_type != 'E')
			{
				int k;
				for(k = 0; k < ((sys_type == 'R' || sys_type == 'S') ? 3 : 7); k++)
				{
					fgets(current_line, 100, fp_nav);
				}
				continue;
			}
		}
		else
		{
			sscanf(current_line, "%d %d %d %d %d %d %lf%le%le%le", &tmpEph.prn, &ye, &mo, &da, &ho, &mi, &sec, &field_1, &field_2, &field_3);
			ye += 2000;
		}

		tmpEph.dAfo = field_1;
		tmpEph.dAf1 = field_2;
		tmpEph.dAf2 = field_3;

		CalToGPS_G(ye, mo, da, ho, mi, sec, (long int *)&tmpEph.nWeekCLK, &tmpEph.dTOC);

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.IODnav = field_1;
		tmpEph.dCrs   = field_2;
		tmpEph.dDn    = field_3;
		tmpEph.dMo    = field_4;

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.dCuc   = field_1;
		tmpEph.dEcc   = field_2;
		tmpEph.dCus   = field_3;
		tmpEph.dAroot = field_4;

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.dtoe    = field_1;
		tmpEph.dCic    = field_2;
		tmpEph.dOMEGAo = field_3;
		tmpEph.dCis    = field_4;

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.d_io      = field_1;
		tmpEph.dCrc      = field_2;
		tmpEph.d_w       = field_3;
		tmpEph.dOMEGADOT = field_4;

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.dIDOT       = field_1;
		int sources        = field_2;
		tmpEph.nWeekNumber = field_3;

		if(tmpEph.nWeekNumber < 1024)
		{
			tmpEph.nWeekNumber += 1024;
		}

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.sSISA = field_1;
		int health   = field_2;
		double bgd1  = field_3;
		double bgd2  = field_4;

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		field_1 = 0; field_2 = 0; field_3 = 0; field_4 = 0;
		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);

		if(field_2 == 0)
		{
			if (sources == 0 || (sources & 0x0001) || (sources & 0x0004)){
				tmpEph.srv = GAL_INAV;
			}
			else if (sources & 0x0002) {
				tmpEph.srv = GAL_FNAV;
			}

			tmpEph.E1B_dvs = ((health & 0x0001) >> 0);
			tmpEph.E1B_hst = ((health & 0x0006) >> 1);
			tmpEph.E5a_dvs = ((health & 0x0008) >> 3);
			tmpEph.E5a_hst = ((health & 0x0030) >> 4);
			tmpEph.E5b_dvs = ((health & 0x0040) >> 6);
			tmpEph.E5b_hst = ((health & 0x0180) >> 7);

			tmpEph.bgd_E1E5a  = bgd1;
			tmpEph.bgd_E1E5b  = bgd2;
		}
		else
		{
			tmpEph.srv = GAL_PRS;
			tmpEph.E1A_hst = health;
			tmpEph.bgd_E1E6A  = bgd1;
		}

		tmpEph.toeWeek = tmpEph.nWeekNumber;
		tmpEph.ToW     = tmpEph.dTOC;

		// FIXME: Define when ephemeris are considered valid based on HS and DVS flags
		tmpEph.vflg = 1;

		double time_difference = calculate_time_difference(tmpEph.nWeekNumber, tmpEph.dtoe, week, gps_time);
		if(!feof(fp_nav) && tmpEph.prn == prn && ((USE_GALPRS && tmpEph.srv == GAL_PRS) || (!USE_GALPRS && tmpEph.srv == GAL_INAV)))
		{
			
			if(fabs(time_difference) <= min_tow_distance)
			{
				min_tow_distance = fabs(time_difference);
				memcpy(ephgal, &tmpEph, sizeof(ephgal_t));

				if(rinex_info.mixed_navigation == 1)
				{
					last_valid_cursor[sys2garlicIndex(prn-1, SYS_GAL)] = ftell(fp_nav);
				}
				else
				{
					last_valid_cursor[prn-1] = ftell(fp_nav);
				}

				if(min_tow_distance < MIN_EPHGAL_DISTANCE)
				{
					break;
				}
			}
		}
	}

}

static void parse_nav_bei_RINEX(FILE *fp_nav, long *last_valid_cursor, int week, double gps_time, int prn, ephbei_t *ephbei)
{
	int ye, mo, da, ho, mi;
	double sec;

	double field_1, field_2, field_3, field_4;

	char current_line[100];

	ephbei_t tmpEph;

	double min_tow_distance = calculate_time_difference(week, gps_time, ephbei->nWeekNumber, ephbei->dtoe);
	while(!feof(fp_nav))
	{
		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);
		if(rinex_info.navbei_version == RINEX_v3)
		{
			char sys_type;
			sscanf(current_line, "%c%d %d %d %d %d %d %lf%le%le%le", &sys_type, &tmpEph.prn, &ye, &mo, &da, &ho, &mi, &sec, &field_1, &field_2, &field_3);
			if(sys_type != 'C')
			{
				int k;
				for(k = 0; k < ((sys_type == 'R' || sys_type == 'S') ? 3 : 7); k++)
				{
					fgets(current_line, 100, fp_nav);
				}
				continue;
			}
		}
		else
		{
			sscanf(current_line, "%d %d %d %d %d %d %lf%le%le%le", &tmpEph.prn, &ye, &mo, &da, &ho, &mi, &sec, &field_1, &field_2, &field_3);
			ye += 2000;
		}

		tmpEph.dAfo = field_1;
		tmpEph.dAf1 = field_2;
		tmpEph.dAf2 = field_3;

		CalToGPS_G(ye, mo, da, ho, mi, sec, (long int *)&tmpEph.nWeekCLK, &tmpEph.dTOC);

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.nIODE = field_1;
		tmpEph.dCrs  = field_2;
		tmpEph.dDn   = field_3;
		tmpEph.dMo   = field_4;

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.dCuc   = field_1;
		tmpEph.dEcc   = field_2;
		tmpEph.dCus   = field_3;
		tmpEph.dAroot = field_4;

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.dtoe_bei= field_1;
		tmpEph.dCic    = field_2;
		tmpEph.dOMEGAo = field_3;
		tmpEph.dCis    = field_4;

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.d_io      = field_1;
		tmpEph.dCrc      = field_2;
		tmpEph.d_w       = field_3;
		tmpEph.dOMEGADOT = field_4;

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.dIDOT       = field_1;
		tmpEph.nWeekNumber = field_3;
		if(tmpEph.nWeekNumber < 1356)
		{
			tmpEph.nWeekNumber += 1356;
		}

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.accuracy  = field_1;
		tmpEph.nHealth   = field_2;
		tmpEph.dTGD_B1B3 = field_3;
		tmpEph.dTGD_B2B3 = field_4;

		fgets(current_line, 100, fp_nav);
		replace_rinex_string(current_line);

		sscanf(current_line, "%le%le%le%le", &field_1, &field_2, &field_3, &field_4);
		tmpEph.ToW    = field_1;
		tmpEph.nIODC  = field_2;

		// Translate time variables from UTC (BDT) to GPS time
		long   bei0_week;
		double bei0_time;
		CalToGPS_G(BDT0[0], BDT0[1], BDT0[2], BDT0[3], BDT0[4], BDT0[5], &bei0_week, &bei0_time);
		int bei_lp  = get_leap_second(bei0_week, bei0_time);
		tmpEph.ToW  += bei_lp;
		tmpEph.dtoe  = tmpEph.dtoe_bei + bei_lp;
		tmpEph.dTOC += bei_lp;

		tmpEph.vflg = (tmpEph.nHealth == 0);

		double time_difference = calculate_time_difference(tmpEph.nWeekNumber, tmpEph.dtoe, week, gps_time);
		if(!feof(fp_nav) && tmpEph.prn == prn)
		{
			if(fabs(time_difference) <= min_tow_distance)
			{
				min_tow_distance = fabs(time_difference);
				memcpy(ephbei, &tmpEph, sizeof(ephbei_t));

				if(rinex_info.mixed_navigation == 1)
				{
					last_valid_cursor[sys2garlicIndex(prn-1, SYS_BEI)] = ftell(fp_nav);
				}
				else
				{
					last_valid_cursor[prn-1] = ftell(fp_nav);
				}

				if(min_tow_distance < MIN_EPHBEI_DISTANCE)
				{
					break;
				}
			}
			else if(tmpEph.dtoe > gps_time)
			{
				break;
			}
		}
	}
}

static void read_next_nav_gps_RINEX(FILE *fp_nav, int week, double gps_time, int prn, ephgps_t *ephgps)
{
	static long int last_valid_cursor[MAXNUMGPS_G] = {0};
	static char first_call = 1;

	if(first_call)
	{
		long int current_cursor = ftell(fp_nav);

		int k;
		for(k = 0; k < MAXNUMGPS_G; k++)
		{
			last_valid_cursor[k] = current_cursor;
		}
		first_call = 0;
	}
	else
	{
		fseek(fp_nav, last_valid_cursor[prn-1], SEEK_SET);
	}
	parse_nav_gps_RINEX(fp_nav, last_valid_cursor, week, gps_time, prn, ephgps);
}

static void read_next_nav_glo_RINEX(FILE *fp_nav, int week, double gps_time, int prn, ephglo_t *ephglo, ionoutc_t *iono_utc)
{
	static long int last_valid_cursor[MAXNUMGLO_G] = {0};
	static char first_call = 1;

	if(first_call)
	{
		long int current_cursor = ftell(fp_nav);

		int k;
		for(k = 0; k < MAXNUMGLO_G; k++)
		{
			last_valid_cursor[k] = current_cursor;
		}
		first_call = 0;
	}
	else
	{
		fseek(fp_nav, last_valid_cursor[prn-1], SEEK_SET);
	}
	parse_nav_glo_RINEX(fp_nav, last_valid_cursor, week, gps_time, prn, ephglo, iono_utc);
}

static void read_next_nav_gal_RINEX(FILE *fp_nav, int week, double gps_time, int prn, ephgal_t *ephgal)
{
	static long int last_valid_cursor[MAXNUMGAL_G] = {0};
	static char first_call = 1;

	if(first_call)
	{
		long int current_cursor = ftell(fp_nav);

		int k;
		for(k = 0; k < MAXNUMGAL_G; k++)
		{
			last_valid_cursor[k] = current_cursor;
		}
		first_call = 0;
	}
	else
	{
		fseek(fp_nav, last_valid_cursor[prn-1], SEEK_SET);
	}
	parse_nav_gal_RINEX(fp_nav, last_valid_cursor, week, gps_time, prn, ephgal);
}

static void read_next_nav_bei_RINEX(FILE *fp_nav, int week, double gps_time, int prn, ephbei_t *ephbei)
{
	static long int last_valid_cursor[MAXNUMBEI_G] = {0};
	static char first_call = 1;

	if(first_call)
	{
		long int current_cursor = ftell(fp_nav);

		int k;
		for(k = 0; k < MAXNUMBEI_G; k++)
		{
			last_valid_cursor[k] = current_cursor;
		}
		first_call = 0;
	}
	else
	{
		fseek(fp_nav, last_valid_cursor[prn-1], SEEK_SET);
	}
	parse_nav_bei_RINEX(fp_nav, last_valid_cursor, week, gps_time, prn, ephbei);
}

static void read_next_nav_multi_RINEX(FILE *fp_nav, int week, double gps_time, int prn, sysFlag_t sys, ephdata_t *ephdata, ionoutc_t *iono_utc)
{
	static long int last_valid_cursor[MAX_N_SAT_G] = {0};

	static char first_call = 1;

	if(first_call)
	{
		long int current_cursor = ftell(fp_nav);

		int k;
		for(k = 0; k < MAX_N_SAT_G; k++)
		{
			last_valid_cursor[k] = current_cursor;
		}
		first_call = 0;
	}
	else
	{
		fseek(fp_nav, last_valid_cursor[sys2garlicIndex(prn-1, sys)], SEEK_SET);
	}

	switch (sys)
	{
	case SYS_GPS:
		parse_nav_gps_RINEX(fp_nav, last_valid_cursor, week, gps_time, prn, &ephdata->GPS);
		break;
	case SYS_GAL:
		parse_nav_gal_RINEX(fp_nav, last_valid_cursor, week, gps_time, prn, &ephdata->GAL);
		break;
	case SYS_GLO:
		parse_nav_glo_RINEX(fp_nav, last_valid_cursor, week, gps_time, prn, &ephdata->GLO, iono_utc);
		break;
	case SYS_BEI:
		parse_nav_bei_RINEX(fp_nav, last_valid_cursor, week, gps_time, prn, &ephdata->BEI);
		break;
	default:
		break;
	}
}

void initialize_ephemeris(FILE *fp_ngps, FILE *fp_nglo, FILE *fp_ngal, FILE *fp_nbei, gnssdata_t *gnssdata)
{
	int prnIndex, ephIndex;
	for(prnIndex = 0; fp_ngps != NULL && prnIndex < MAXNUMGPS_G; prnIndex++)
	{
		ephIndex = sys2garlicIndex(prnIndex, SYS_GPS);
		if(rinex_info.mixed_navigation == 0)
		{
			read_next_nav_gps_RINEX(fp_ngps, rinex_info.first_week, rinex_info.first_tow, prnIndex+1, &gnssdata->EPH[ephIndex].GPS);
		}
		else
		{
			read_next_nav_multi_RINEX(fp_ngps, rinex_info.first_week, rinex_info.first_tow, prnIndex+1, SYS_GPS, &gnssdata->EPH[ephIndex], &gnssdata->iono_utc);
		}
	}
	for(prnIndex = 0; fp_ngal != NULL && prnIndex < MAXNUMGAL_G; prnIndex++)
	{
		ephIndex = sys2garlicIndex(prnIndex, SYS_GAL);
		if(rinex_info.mixed_navigation == 0)
		{
			read_next_nav_gal_RINEX(fp_ngal, rinex_info.first_week, rinex_info.first_tow, prnIndex+1, &gnssdata->EPH[ephIndex].GAL);
		}
		else
		{
			read_next_nav_multi_RINEX(fp_ngal, rinex_info.first_week, rinex_info.first_tow, prnIndex+1, SYS_GAL, &gnssdata->EPH[ephIndex], &gnssdata->iono_utc);
		}
	}
	for(prnIndex = 0; fp_nglo != NULL && prnIndex < MAXNUMGLO_G; prnIndex++)
	{
		ephIndex = sys2garlicIndex(prnIndex, SYS_GLO);
		if(rinex_info.mixed_navigation == 0)
		{
			read_next_nav_glo_RINEX(fp_nglo, rinex_info.first_week, rinex_info.first_tow, prnIndex+1, &gnssdata->EPH[ephIndex].GLO, &gnssdata->iono_utc);
		}
		else
		{
			read_next_nav_multi_RINEX(fp_nglo, rinex_info.first_week, rinex_info.first_tow, prnIndex+1, SYS_GLO, &gnssdata->EPH[ephIndex], &gnssdata->iono_utc);
		}
	}
	for(prnIndex = 0; fp_nbei != NULL && prnIndex < MAXNUMBEI_G; prnIndex++)
	{
		ephIndex = sys2garlicIndex(prnIndex, SYS_BEI);
		if(rinex_info.mixed_navigation == 0)
		{
			read_next_nav_bei_RINEX(fp_nbei, rinex_info.first_week, rinex_info.first_tow, prnIndex+1, &gnssdata->EPH[ephIndex].BEI);
		}
		else
		{
			read_next_nav_multi_RINEX(fp_nbei, rinex_info.first_week, rinex_info.first_tow, prnIndex+1, SYS_BEI, &gnssdata->EPH[ephIndex], &gnssdata->iono_utc);
		}
	}
}

void read_next_obs_RINEX_v2(FILE *fp_obs, gnssdata_t *gnssdata)
{
	int ye, mo, da, ho, mi;
	double sec;

	int  number_of_measurements;
	int  epoch_flag;

	char sat_id[MAX_CHANNELS_G] = {0};
	int  sat_no[MAX_CHANNELS_G] = {0};

	char current_line[100], current_field[16];

	memset(current_line, 0, 100);
	fgets(current_line, 100, fp_obs);

	sscanf(current_line, "%d %d %d %d %d %lf %d %d %c%d%c%d%c%d%c%d%c%d%c%d%c%d%c%d%c%d%c%d%c%d%c%d\n",
			&ye, &mo, &da, &ho, &mi, &sec, &epoch_flag, &number_of_measurements,
			&sat_id[0], &sat_no[0], &sat_id [1], &sat_no [1], &sat_id [2], &sat_no [2],
			&sat_id[3], &sat_no[3], &sat_id [4], &sat_no [4], &sat_id [5], &sat_no [5],
			&sat_id[6], &sat_no[6], &sat_id [7], &sat_no [7], &sat_id [8], &sat_no [8],
			&sat_id[9], &sat_no[9], &sat_id[10], &sat_no[10], &sat_id[11], &sat_no[11]);

	if(number_of_measurements > 12)
	{
		memset(current_line, 0, 100);
		fgets(current_line, 100, fp_obs);
		sscanf(current_line, " %c%d%c%d%c%d%c%d%c%d%c%d%c%d%c%d%c%d%c%d%c%d%c%d\n",
				&sat_id[12], &sat_no[12], &sat_id[13], &sat_no[13], &sat_id[14], &sat_no[14],
				&sat_id[15], &sat_no[15], &sat_id[16], &sat_no[16], &sat_id[17], &sat_no[17],
				&sat_id[18], &sat_no[18], &sat_id[19], &sat_no[19], &sat_id[20], &sat_no[20],
				&sat_id[21], &sat_no[21], &sat_id[22], &sat_no[22], &sat_id[23], &sat_no[23]);
	}
#if MAX_CHANNELS_G == 36
	if(number_of_measurements > 24)
	{
		memset(current_line, 0, 100);
		fgets(current_line, 100, fp_obs);
		sscanf(current_line, " %c%d%c%d%c%d%c%d%c%d%c%d%c%d%c%d%c%d%c%d%c%d%c%d\n",
				&sat_id[24], &sat_no[24], &sat_id[25], &sat_no[25], &sat_id[26], &sat_no[26],
				&sat_id[27], &sat_no[27], &sat_id[28], &sat_no[28], &sat_id[29], &sat_no[29],
				&sat_id[30], &sat_no[30], &sat_id[31], &sat_no[31], &sat_id[32], &sat_no[32],
				&sat_id[33], &sat_no[33], &sat_id[34], &sat_no[34], &sat_id[35], &sat_no[35]);
	}
#endif

	CalToGPS_G(2000+ye, mo, da, ho, mi, sec, (long int *)&gnssdata->week, &gnssdata->tow);

	if(!feof(fp_obs))
	{
		int j, k;
		int sys;
		memset(gnssdata->OBS, 0, sizeof(gnssdata->OBS));

		for(j = 0; j < number_of_measurements; j++)
		{
			sys = getConstellationFromRINEX(sat_id[j]);

			int sig_last = -1, sig;
			for(k = 0; k < rinex_info.number_of_observations[sys]; k++)
			{
				char obs_code[3];
				obs_code[0] = rinex_info.observation_list[sys][k][1];
				if (obs_code[0] == '1') { obs_code[1] = 'C'; }
				else if (obs_code[0] == '5' || obs_code[0] == '7') { obs_code[1] = 'Q'; }
				obs_code[2] = '\0';
				sig = getServiceFromRINEX(sat_id[j], obs_code);
				if (sig < 0) { continue; }

				// To support multi-frequency add new observation in the list if the observation type has changed
				if (sig != sig_last)
				{
					add_new_observation(gnssdata, sig, sat_no[j]);
					sig_last = sig;
				}

				if(k % 5 == 0)
				{
					memset(current_line, 0, 100);
					fgets(current_line, 100, fp_obs);
				}

				// The total size of the field is 16, although 2 chars are dedicated to monitoring info.
				strncpy(current_field, &current_line[16*(k % 5)], 16);

				process_obs_field(current_field, rinex_info.observation_list[sys][k][0], &gnssdata->OBS[gnssdata->noOfChannelsAv-1]);
			}
		}
		memcpy(&last_obs, &gnssdata->OBS, sizeof(last_obs));
	}
}

void read_next_obs_RINEX_v3(FILE *fp_obs, gnssdata_t *gnssdata)
{
	int ye, mo, da, ho, mi;
	int  number_of_satellites, epoch_flag, sat_no;
	double sec, clock_offset;
	char sys_code,current_line[MAXLENGTHRV3], current_field[16];

	// Read first line of the epoch including the time and number of meas
	memset(gnssdata->OBS, 0, sizeof(gnssdata->OBS));
	fgets(current_line, MAXLENGTHRV3, fp_obs);
	sscanf(current_line, "> %d %d %d %d %d %lf %d %d %lf\n", &ye, &mo, &da, &ho, &mi, &sec, &epoch_flag, &number_of_satellites, &clock_offset);
	CalToGPS_G(ye, mo, da, ho, mi, sec, (long int *)&gnssdata->week, &gnssdata->tow);

	// Read observation lines
	if(!feof(fp_obs))
	{
		int j, k;
		int sys;

		for(j = 0; j < number_of_satellites; j++)
		{
			memset(current_line, 0, MAXLENGTHRV3);
			fgets (current_line, MAXLENGTHRV3, fp_obs);

			sscanf(current_line,"%c%d", &sys_code, &sat_no);

			sys = getConstellationFromRINEX(sys_code);

			// FIXME: skipping GEO BeiDou Satellites
			if (sys == SYS_BEI && ( (sat_no >= 1 && sat_no <= 5) || sat_no == 13 || sat_no > MAXNUMBEI_G)) { continue; }

			if (sys >= 0)
			{
				int sig_last = -1, sig;
				for(k = 0; k < rinex_info.number_of_observations[sys]; k++)
				{
					char obs_code[3];
					obs_code[0] = rinex_info.observation_list[sys][k][1];
					obs_code[1] = rinex_info.observation_list[sys][k][2];
					obs_code[2] = '\0';
					sig = getServiceFromRINEX(sys_code, obs_code);
					if (sig < 0) { continue; }

					// To support multi-frequency add new observation in the list if the observation type has changed
					if (sig != sig_last)
					{
						add_new_observation(gnssdata, sig, sat_no);
						sig_last = sig;
					}

					// The total size of the field is 16, although 2 chars are dedicated to monitoring info.
					strncpy(current_field, &current_line[4 + 16*k], 16);
					process_obs_field(current_field, rinex_info.observation_list[sys][k][0], &gnssdata->OBS[gnssdata->noOfChannelsAv-1]);
				}
			}
		}
		memcpy(&last_obs, &gnssdata->OBS, sizeof(last_obs));
	}
}

void check_ephemeris_status(FILE *fp_ngps, FILE *fp_nglo, FILE *fp_ngal, FILE *fp_nbei, gnssdata_t *gnssdata,
		int sat_new_in, int epoch_update_eph, int new_sat[MAX_CHANNELS_G])
{
	static char unvalid_ephemeris[MAX_N_SAT_G] = {0};

	int k, prn, i, sat_found;
	for(k = 0; k < gnssdata->noOfChannelsAv; k++)
	{

		/* if this function is executed as a new sat came in (and not just because it's the routinely per-minute ephemeris update),
		update onyl ephemeris of new satellites */
		if( sat_new_in > 0 && !epoch_update_eph ){
			sat_found = 0;
			for( i = 0; i < sat_new_in; i++ ){
				if( gnssdata->OBS[k].PRN == new_sat[i] ){
					sat_found = 1;
				}
			}
			if( !sat_found ){
				continue;
			}

		}

		prn = gnssdata->OBS[k].PRN - 1;
		if(signal2system(gnssdata->OBS[k].sigFlag) == SYS_GPS && fp_ngps != NULL && unvalid_ephemeris[prn] == 0)
		{
			double diff_tow = calculate_time_difference(gnssdata->week, gnssdata->tow, gnssdata->EPH[prn].GPS.nWeekNumber, gnssdata->EPH[prn].GPS.dtoe);

			if(diff_tow > MIN_EPHGPS_DISTANCE)
			{
				if(rinex_info.mixed_navigation == 0)
				{
					read_next_nav_gps_RINEX(fp_ngps, gnssdata->week, gnssdata->tow, garlic2sysIndex(prn,SYS_GPS)+1, &gnssdata->EPH[prn].GPS);
				}
				else
				{
					read_next_nav_multi_RINEX(fp_ngps, gnssdata->week, gnssdata->tow, garlic2sysIndex(prn,SYS_GPS)+1, SYS_GPS, &gnssdata->EPH[prn], &gnssdata->iono_utc);
				}
				diff_tow = calculate_time_difference(gnssdata->week, gnssdata->tow, gnssdata->EPH[prn].GPS.nWeekNumber, gnssdata->EPH[prn].GPS.dtoe);
				if(diff_tow > 4*MIN_EPHGPS_DISTANCE)
				{
					unvalid_ephemeris[prn] = 1;

					gnssdata->EPH[prn].GPS.vflg = 0;
					gnssdata->OBS[k].lock = 0;

					printf("WARNING: Use of PRN G%02d has been disabled due to lapsed ephemeris at TOW=%d.\n", garlic2sysIndex(prn,SYS_GPS)+1, (int)gnssdata->tow);
				}
			}
		}

		else if(signal2system(gnssdata->OBS[k].sigFlag) == SYS_GAL && fp_ngal != NULL && unvalid_ephemeris[prn] == 0)
		{
			double diff_tow = calculate_time_difference(gnssdata->week, gnssdata->tow, gnssdata->EPH[prn].GAL.nWeekNumber, gnssdata->EPH[prn].GAL.dtoe);

			if(diff_tow > MIN_EPHGAL_DISTANCE)
			{

				if(rinex_info.mixed_navigation == 0)
				{
					read_next_nav_gal_RINEX(fp_ngal, gnssdata->week, gnssdata->tow, garlic2sysIndex(prn,SYS_GAL)+1, &gnssdata->EPH[prn].GAL);
				}
				else
				{
					read_next_nav_multi_RINEX(fp_ngal, gnssdata->week, gnssdata->tow, garlic2sysIndex(prn,SYS_GAL)+1, SYS_GAL, &gnssdata->EPH[prn], &gnssdata->iono_utc);
				}
				diff_tow = calculate_time_difference(gnssdata->week, gnssdata->tow, gnssdata->EPH[prn].GAL.nWeekNumber, gnssdata->EPH[prn].GAL.dtoe);
				if(diff_tow > 4*MIN_EPHGAL_DISTANCE)
				{
					unvalid_ephemeris[prn] = 1;

					gnssdata->EPH[prn].GAL.vflg = 0;
					gnssdata->OBS[k].lock = 0;

					printf("WARNING: Use of PRN E%02d has been disabled due to lapsed ephemeris at TOW=%d.\n", garlic2sysIndex(prn,SYS_GAL)+1, (int)gnssdata->tow);
				}
			}
		}
		else if(signal2system(gnssdata->OBS[k].sigFlag) ==  SYS_GLO && fp_nglo != NULL && unvalid_ephemeris[prn] == 0)
		{
			double diff_tow = calculate_time_difference(gnssdata->week, gnssdata->tow, gnssdata->EPH[prn].GLO.nWeekNumber, gnssdata->EPH[prn].GLO.dtoe);

			if(diff_tow > MIN_EPHGLO_DISTANCE)
			{
				if(rinex_info.mixed_navigation == 0)
				{
					read_next_nav_glo_RINEX(fp_nglo, gnssdata->week, gnssdata->tow, garlic2sysIndex(prn,SYS_GLO)+1, &gnssdata->EPH[prn].GLO, &gnssdata->iono_utc);
				}
				else
				{
					read_next_nav_multi_RINEX(fp_nglo, gnssdata->week, gnssdata->tow, garlic2sysIndex(prn,SYS_GLO)+1, SYS_GLO, &gnssdata->EPH[prn], &gnssdata->iono_utc);
				}
				diff_tow = calculate_time_difference(gnssdata->week, gnssdata->tow, gnssdata->EPH[prn].GLO.nWeekNumber, gnssdata->EPH[prn].GLO.dtoe);
				if(diff_tow > 4*MIN_EPHGLO_DISTANCE)
				{
					unvalid_ephemeris[prn] = 1;

					gnssdata->EPH[prn].GLO.vflg = 0;
					gnssdata->OBS[k].lock = 0;

					printf("WARNING: Use of PRN R%02d has been disabled due to lapsed ephemeris at TOW=%d.\n", garlic2sysIndex(prn,SYS_GLO)+1, (int)gnssdata->tow);
				}
			}
		}

		else if(signal2system(gnssdata->OBS[k].sigFlag) ==  SYS_BEI && fp_nbei != NULL && unvalid_ephemeris[prn] == 0)
		{
			double diff_tow = calculate_time_difference(gnssdata->week, gnssdata->tow, gnssdata->EPH[prn].BEI.nWeekNumber, gnssdata->EPH[prn].BEI.dtoe);

			if(diff_tow > MIN_EPHBEI_DISTANCE)
			{
				if(rinex_info.mixed_navigation == 0)
				{
					read_next_nav_bei_RINEX(fp_ngps, gnssdata->week, gnssdata->tow, garlic2sysIndex(prn,SYS_BEI)+1, &gnssdata->EPH[prn].BEI);
				}
				else
				{
					read_next_nav_multi_RINEX(fp_ngps, gnssdata->week, gnssdata->tow, garlic2sysIndex(prn,SYS_BEI)+1, SYS_BEI, &gnssdata->EPH[prn], &gnssdata->iono_utc);
				}
				diff_tow = calculate_time_difference(gnssdata->week, gnssdata->tow, gnssdata->EPH[prn].BEI.nWeekNumber, gnssdata->EPH[prn].BEI.dtoe);
				if(diff_tow > 4*MIN_EPHBEI_DISTANCE)
				{
					unvalid_ephemeris[prn] = 1;

					gnssdata->EPH[prn].BEI.vflg = 0;
					gnssdata->OBS[k].lock = 0;

					printf("WARNING: Use of PRN C%02d has been disabled due to lapsed ephemeris at TOW=%d.\n", garlic2sysIndex(prn,SYS_BEI)+1, (int)gnssdata->tow);
				}
			}
		}
	}
}

//void check_ephemeris_status(FILE *fp_ngps, FILE *fp_nglo, FILE *fp_ngal, FILE *fp_nbei, gnssdata_t *gnssdata)
//{
//	static char unvalid_ephemeris[MAX_N_SAT_G] = {0};
//
//	int k, prn;
//	for(k = 0; k < gnssdata->noOfChannelsAv; k++)
//	{
//		prn = gnssdata->OBS[k].PRN - 1;
//		if(signal2system(gnssdata->OBS[k].sigFlag) == SYS_GPS && fp_ngps != NULL && unvalid_ephemeris[prn] == 0)
//		{
//			double diff_tow = calculate_time_difference(gnssdata->week, gnssdata->tow, gnssdata->EPH[prn].GPS.nWeekNumber, gnssdata->EPH[prn].GPS.dtoe);
//
//			if(diff_tow > MIN_EPHGPS_DISTANCE)
//			{
//				if(rinex_info.mixed_navigation == 0)
//				{
//					read_next_nav_gps_RINEX(fp_ngps, gnssdata->week, gnssdata->tow, garlic2sysIndex(prn,SYS_GPS)+1, &gnssdata->EPH[prn].GPS);
//				}
//				else
//				{
//					read_next_nav_multi_RINEX(fp_ngps, gnssdata->week, gnssdata->tow, garlic2sysIndex(prn,SYS_GPS)+1, SYS_GPS, &gnssdata->EPH[prn], &gnssdata->iono_utc);
//				}
//				diff_tow = calculate_time_difference(gnssdata->week, gnssdata->tow, gnssdata->EPH[prn].GPS.nWeekNumber, gnssdata->EPH[prn].GPS.dtoe);
//				if(diff_tow > 4*MIN_EPHGPS_DISTANCE)
//				{
//					unvalid_ephemeris[prn] = 1;
//
//					gnssdata->EPH[prn].GPS.vflg = 0;
//					gnssdata->OBS[k].lock = 0;
//
//					printf("WARNING: Use of PRN G%02d has been disabled due to lapsed ephemeris at TOW=%d.\n", garlic2sysIndex(prn,SYS_GPS)+1, (int)gnssdata->tow);
//				}
//			}
//		}
//
//		else if(signal2system(gnssdata->OBS[k].sigFlag) == SYS_GAL && fp_ngal != NULL && unvalid_ephemeris[prn] == 0)
//		{
//			double diff_tow = calculate_time_difference(gnssdata->week, gnssdata->tow, gnssdata->EPH[prn].GAL.nWeekNumber, gnssdata->EPH[prn].GAL.dtoe);
//
//			if(diff_tow > MIN_EPHGAL_DISTANCE)
//			{
//
//				if(rinex_info.mixed_navigation == 0)
//				{
//					read_next_nav_gal_RINEX(fp_ngal, gnssdata->week, gnssdata->tow, garlic2sysIndex(prn,SYS_GAL)+1, &gnssdata->EPH[prn].GAL);
//				}
//				else
//				{
//					read_next_nav_multi_RINEX(fp_ngal, gnssdata->week, gnssdata->tow, garlic2sysIndex(prn,SYS_GAL)+1, SYS_GAL, &gnssdata->EPH[prn], &gnssdata->iono_utc);
//				}
//				diff_tow = calculate_time_difference(gnssdata->week, gnssdata->tow, gnssdata->EPH[prn].GAL.nWeekNumber, gnssdata->EPH[prn].GAL.dtoe);
//				if(diff_tow > 4*MIN_EPHGAL_DISTANCE)
//				{
//					unvalid_ephemeris[prn] = 1;
//
//					gnssdata->EPH[prn].GAL.vflg = 0;
//					gnssdata->OBS[k].lock = 0;
//
//					printf("WARNING: Use of PRN E%02d has been disabled due to lapsed ephemeris at TOW=%d.\n", garlic2sysIndex(prn,SYS_GAL)+1, (int)gnssdata->tow);
//				}
//			}
//		}
//		else if(signal2system(gnssdata->OBS[k].sigFlag) ==  SYS_GLO && fp_nglo != NULL && unvalid_ephemeris[prn] == 0)
//		{
//			double diff_tow = calculate_time_difference(gnssdata->week, gnssdata->tow, gnssdata->EPH[prn].GLO.nWeekNumber, gnssdata->EPH[prn].GLO.dtoe);
//
//			if(diff_tow > MIN_EPHGLO_DISTANCE)
//			{
//				if(rinex_info.mixed_navigation == 0)
//				{
//					read_next_nav_glo_RINEX(fp_nglo, gnssdata->week, gnssdata->tow, garlic2sysIndex(prn,SYS_GLO)+1, &gnssdata->EPH[prn].GLO, &gnssdata->iono_utc);
//				}
//				else
//				{
//					read_next_nav_multi_RINEX(fp_nglo, gnssdata->week, gnssdata->tow, garlic2sysIndex(prn,SYS_GLO)+1, SYS_GLO, &gnssdata->EPH[prn], &gnssdata->iono_utc);
//				}
//				diff_tow = calculate_time_difference(gnssdata->week, gnssdata->tow, gnssdata->EPH[prn].GLO.nWeekNumber, gnssdata->EPH[prn].GLO.dtoe);
//				if(diff_tow > 4*MIN_EPHGLO_DISTANCE)
//				{
//					unvalid_ephemeris[prn] = 1;
//
//					gnssdata->EPH[prn].GLO.vflg = 0;
//					gnssdata->OBS[k].lock = 0;
//
//					printf("WARNING: Use of PRN R%02d has been disabled due to lapsed ephemeris at TOW=%d.\n", garlic2sysIndex(prn,SYS_GLO)+1, (int)gnssdata->tow);
//				}
//			}
//		}
//
//		else if(signal2system(gnssdata->OBS[k].sigFlag) ==  SYS_BEI && fp_nbei != NULL && unvalid_ephemeris[prn] == 0)
//		{
//			double diff_tow = calculate_time_difference(gnssdata->week, gnssdata->tow, gnssdata->EPH[prn].BEI.nWeekNumber, gnssdata->EPH[prn].BEI.dtoe);
//
//			if(diff_tow > MIN_EPHBEI_DISTANCE)
//			{
//				if(rinex_info.mixed_navigation == 0)
//				{
//					read_next_nav_bei_RINEX(fp_ngps, gnssdata->week, gnssdata->tow, garlic2sysIndex(prn,SYS_BEI)+1, &gnssdata->EPH[prn].BEI);
//				}
//				else
//				{
//					read_next_nav_multi_RINEX(fp_ngps, gnssdata->week, gnssdata->tow, garlic2sysIndex(prn,SYS_BEI)+1, SYS_BEI, &gnssdata->EPH[prn], &gnssdata->iono_utc);
//				}
//				diff_tow = calculate_time_difference(gnssdata->week, gnssdata->tow, gnssdata->EPH[prn].BEI.nWeekNumber, gnssdata->EPH[prn].BEI.dtoe);
//				if(diff_tow > 4*MIN_EPHBEI_DISTANCE)
//				{
//					unvalid_ephemeris[prn] = 1;
//
//					gnssdata->EPH[prn].BEI.vflg = 0;
//					gnssdata->OBS[k].lock = 0;
//
//					printf("WARNING: Use of PRN C%02d has been disabled due to lapsed ephemeris at TOW=%d.\n", garlic2sysIndex(prn,SYS_BEI)+1, (int)gnssdata->tow);
//				}
//			}
//		}
//	}
//}
