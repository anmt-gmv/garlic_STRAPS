///////////////////////////////////////////////////////////////////////////////
//
// File     : garlic_functions.h.c
// Purpose  : The file includes several functions to define an independent
//            task in Teseo II. The task is able to handle with GARLIC project
//			  navigation structures, and calculate a proprietary navigation
//            solution.
//
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE                |
//          |         |          |      |                                     |
//          |   1.0   | 12/11/29 | cmvv |            First version            |
//
///////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "dbtypedef.h"
#include "atmsphcorr.h"
#include "garlic_functions.h"
#include "GICSRxAux.h"
#include "GICSRxObs.h"

/* This function is used for the initialization of the GARLIC navigation data.
 * Currently, the configuration parameters for the different algorithms is pre-set.
 */
#if USE_EPHEMERIS_PROP == 1
void initialize_garlic_data(gnssdata_t *gnssdata, gicsrx_t *gicsrxdata, propephem_t *propephem, lsq_state_t *lsq_state, lsq_matrix_t *lsq_mat)
#else
void initialize_garlic_data(gnssdata_t *gnssdata, gicsrx_t *gicsrxdata, lsq_state_t *lsq_state, lsq_matrix_t *lsq_mat, char **argv)
#endif
{
	// Loop variable.
	int k;

	memset(gicsrxdata, 0, sizeof(gicsrx_t));
	memset(gnssdata,   0, sizeof(gnssdata_t));

	// Set time of last reset to an invalid tow.
	gicsrxdata->exchange.timeOfLastReset = -1;

	memset(lsq_state, 0, sizeof(lsq_state_t));
	memset(lsq_mat,   0, sizeof(lsq_matrix_t));

#if USE_EPHEMERIS_PROP == 1
	propephem->currChannel = 0;
#endif
	for(k = 0; k < MAX_N_SAT_G; k++)
	{
		if(k >= 0 && k < MAXNUMGPS_G)
		{
			gnssdata->EPH[k].GPS.nIODE   = -1;
			gnssdata->EPH[k].GPS.nIODC   = -1;
			gnssdata->EPH[k].GPS.nHealth = -1;
		}
		else if(k >= MAXNUMGPS_G && k < MAXNUMGPS_G + MAXNUMGAL_G)
		{
			gnssdata->EPH[k].GAL.IODnav  = -1;
			gnssdata->EPH[k].GAL.E1B_hst = 1;
			gnssdata->EPH[k].GAL.E1A_hst = 1;
			gnssdata->EPH[k].GAL.E5a_hst = 1;
			gnssdata->EPH[k].GAL.E5b_hst = 1;
			gnssdata->EPH[k].GAL.E6A_hst = 1;
		}
		else if(k >= MAXNUMGPS_G + MAXNUMGAL_G && k < MAXNUMGPS_G + MAXNUMGAL_G + MAXNUMGLO_G)
		{
			gnssdata->EPH[k].GLO.tk = 86400;
			gnssdata->EPH[k].GLO.nHealth = -1;
		}
		else if(k >= MAXNUMGPS_G + MAXNUMGAL_G + MAXNUMGLO_G && k < MAXNUMGPS_G + MAXNUMGAL_G + MAXNUMGLO_G + MAXNUMBEI_G)
		{
			gnssdata->EPH[k].BEI.nIODE   = -1;
			gnssdata->EPH[k].BEI.nIODC   = -1;
			gnssdata->EPH[k].BEI.nHealth = -1;
		}

#if USE_EPHEMERIS_PROP == 1
		// Mark propagated satellite ephemeris as unavailable.
		propephem->tow_node[k] = -DAYSECS;
#endif
	}

	// Set default leap seconds.
	gnssdata->iono_utc.gps.dtls = LEAPSECONDS;
	gnssdata->iono_utc.gal.dtls = LEAPSECONDS;

	// In this version, the configuration parameters of the algorithms are pre-defined.
	//
	// Mask angle.
	gicsrxdata->kconf.mask_angle = 5.0;
	//
	// Minimum threshold for carrier to noise ratio in measurements.
	gicsrxdata->kconf.CN0_thr = 25.0;  // this was 25 in princple
//	printf("ATTENTION: CN0 threshold changed!!\n");

	// Config the type of measurement to be used, Pilot/Data
	gicsrxdata->kconf.use_pilot            = 1;
	gicsrxdata->kconf.use_data             = 1;

	// Config constellations to be used in the PVT
	gicsrxdata->kconf.use_sys[SYS_GPS]     = atoi(argv[3]);
	gicsrxdata->kconf.use_sys[SYS_GAL]     = atoi(argv[4]);
	gicsrxdata->kconf.use_sys[SYS_GLO]     = atoi(argv[5]);
	gicsrxdata->kconf.use_sys[SYS_BEI]     = atoi(argv[6]);

	// Config the frequency bands to be used in the PVT
	gicsrxdata->kconf.use_freq[FREQ_L1E1]  = 1;
	gicsrxdata->kconf.use_freq[FREQ_L5]    = 0;
	gicsrxdata->kconf.use_freq[FREQ_E5a]   = 0;
	gicsrxdata->kconf.use_freq[FREQ_E5b]   = 0;

	// Config the type of ionospheric model to be applied
	gicsrxdata->kconf.iono_model = ION_NONE;
//	gicsrxdata->kconf.iono_model = ION_FREE;

	// Size of observation window.
	gicsrxdata->kconf.win_PR = 10;
	//
	// Tolerance in failed number of epochs.
	gicsrxdata->kconf.tol_PR = 2;
	//
	// Minimum percentage of PR measurements that must be accepted to consider valid epoch.
	gicsrxdata->kconf.per_PR = 30.0;
	//
	// Maximum number of epochs during which propagation can be performed. Reset if higher.
	gicsrxdata->kconf.maxKalx = 30;

#if (USE_IMU == 0 || CALC_GNSS_ONLY == 1)
	// Process noise variance for speed.
	gicsrxdata->kconf.sigma2_vel = 1e-2;
#endif
#define CLOCK_SIGMA 100
	// Process noise variance for drift in clock.
	gicsrxdata->kconf.sigma2_clk = CLOCK_SIGMA;

	int i;
	// Process noise variance for inter-system biases.
	for(i=0;i<NUM_SYS-1;i++) {
		gicsrxdata->kconf.sigma2_isb[i] = CLOCK_SIGMA;
	}

	// Process noise variance for inter-frequency biases.
	for(i=0;i<NUM_FREQ-1;i++) {
		gicsrxdata->kconf.sigma2_ifreqb[i] = 100;
	}

	// Value for the threshold in chi-square test.
	gicsrxdata->kconf.chi_test1 = 100;

#if CALC_GNSS_ONLY == 1
	gicsrxdata->kconf.chi_test2 = 25;
#endif

	// Minimum required velocity to initialize GNSS+IMU KF.
	gicsrxdata->kconf.hkvel0 = 0.00;

	// Minimum number of satellites required to initialize GNSS-only KF.
	// The GNSS+IMU KF is initialized with 2 additional satellites in view.
	gicsrxdata->kconf.noSats_min = 5;

	// Number of epochs during which speed smoothing is not applied for better filter convergence.
	gicsrxdata->kconf.ktimeNoSmooth = 60;

	// Threshold velocity to consider the vehicle moving: used for weighting ranges.
	gicsrxdata->kconf.kvelStopThr = 1;

	// Factor used to smooth the effect of the pseudo-range measurements while stopped.
	gicsrxdata->kconf.kvelStopFactor = 1;

#if USE_WEIGHTING_MODEL == 0
	// PR measurement noise variance in Kalman filter.
	gicsrxdata->kconf.sigma_prange = 5.00;
	//
	// Doppler measurement noise variance in Kalman filter.
	gicsrxdata->kconf.sigma_doppler = 0.50;
#endif

	// Indicate the frequency that GNSS data is dumped (in seconds).
	gicsrxdata->kconf.GNSS_rate = 1;

#if USE_IMU == 1

	gicsrxdata->kconf.IMU_rate     = 1;
	gicsrxdata->kconf.sigma2_bacc  = 1e-2;
	gicsrxdata->kconf.sigma2_bgyr  = 1e-4;
	gicsrxdata->kconf.sigma2_wnacc = 1e-2;
	gicsrxdata->kconf.sigma2_wngyr = 1e-4;
	gicsrxdata->kconf.tau_acc      = 1e4;
	gicsrxdata->kconf.tau_gyr      = 1e4;

	double drot_x, drot_y, drot_z;

#if CALC_INSTALL_MAT == 1
	drot_x = drot_y = drot_z = 0;
#else
	char set_config = 'B';

	switch(set_config)
	{
	case 'A':
		// Matriz de instalacion para placa post-proceso.
		drot_x =  +2.312*DEGTORAD;
		drot_y =  -9.263*DEGTORAD;
		drot_z = +11.891*DEGTORAD;
		break;

	case 'C':
		// Matriz de instalacion para placa real-time.
		drot_x = +0.634*DEGTORAD;
		drot_y = +8.523*DEGTORAD;
		drot_z = +9.924*DEGTORAD;
		break;

	case 'B':
	default:
		drot_x = drot_y = drot_z = 0;
	}
#endif

	// Configure installation matrix.
	gicsrxdata->kconf.Imat[0][0] = +cos(drot_y)*cos(drot_z);
	gicsrxdata->kconf.Imat[0][1] = -cos(drot_x)*sin(drot_z) + sin(drot_x)*sin(drot_y)*cos(drot_z);
	gicsrxdata->kconf.Imat[0][2] = +sin(drot_x)*sin(drot_z) + cos(drot_x)*sin(drot_y)*cos(drot_z);
	gicsrxdata->kconf.Imat[1][0] = +cos(drot_y)*sin(drot_z);
	gicsrxdata->kconf.Imat[1][1] = +cos(drot_x)*cos(drot_z) + sin(drot_x)*sin(drot_y)*sin(drot_z);
	gicsrxdata->kconf.Imat[1][2] = -sin(drot_x)*cos(drot_z) + cos(drot_x)*sin(drot_y)*sin(drot_z);
	gicsrxdata->kconf.Imat[2][0] = -sin(drot_y);
	gicsrxdata->kconf.Imat[2][1] = +sin(drot_x)*cos(drot_y);
	gicsrxdata->kconf.Imat[2][2] = +cos(drot_x)*cos(drot_y);
#endif

#if CALCULATE_IBPL == 1
	// IBPL
	gicsrxdata->kconf.errors_absolute_sum       = 1;
	gicsrxdata->kconf.gamma_corr_prange         = 0.9;
	gicsrxdata->kconf.gamma_corr_doppler        = 0.5;
	gicsrxdata->kconf.dof_tail_depth            = 3;
	gicsrxdata->kconf.max_dof                   = 16;
	gicsrxdata->kconf.prange_dop_factor         = 1;
	gicsrxdata->kconf.doppler_dop_factor        = 1;
	gicsrxdata->kconf.gamma_acc_indicator       = 0.8;
	gicsrxdata->kconf.acc_indicator_factor_gnss = 0.150;
	gicsrxdata->kconf.acc_indicator_factor_imu  = 0.050;
#endif

#if USE_IMU == 1
	gicsrxdata->exchange.quaternion[3] = 1;
#endif
}

void filter_obs_by_config(gnssdata_t *gnssdata, kconf_t *conf)
{
	int i, ind = 0;
	for(i = 0; i < gnssdata->noOfChannelsAv; i++)
	{
//		if( gnssdata->OBS[i].PRN == 77 ){
//			continue;
//		}

		int  sys = signal2system(gnssdata->OBS[i].sigFlag);
		int freq = signal2freq(gnssdata->OBS[i].sigFlag);

		if ( sys  >= 0 && conf->use_sys[sys] &&									// Constellation configured
			 freq >= 0 && conf->use_freq[freq] && 								// Frequency configured
			 ((conf->use_pilot && isPilotSignal(gnssdata->OBS[i].sigFlag)) ||	// Pilot signals configured
			  (conf->use_data  && isDataSignal (gnssdata->OBS[i].sigFlag)) ))	// Data  signals configured
		{
	    	gnssdata->OBS[ind] = gnssdata->OBS[i];
	    	gnssdata->SAT_ID[ind++] = gnssdata->SAT_ID[i];
		}
		if (ind-1 != i)
		{
			memset(&gnssdata->OBS[i], 0, sizeof(obsdata_t));
			gnssdata->SAT_ID[i] = 0;
		}
	}
	gnssdata->noOfChannelsAv = ind;
}

void combine_obs(gnssdata_t *gnssdata, kconf_t *conf)
{
	int i, index_map[NUM_FREQ][MAX_N_SAT_G];


	if (conf->iono_model != ION_FREE) { return; }

	// Get the list of available frequencies per satellite storing the index to observation list.
	memset(index_map, -1, sizeof(index_map));

	for(i = 0; i < gnssdata->noOfChannelsAv; i++)
	{
		if (gnssdata->OBS[i].lock == 1 && gnssdata->OBS[i].prange_status == OBS_VALID)
		{
			int prn, freq;
			// if pilot signals usage is activated, it is preferred over the data
			// components to perform iono-free combination
			if (conf->use_pilot && isPilotSignal(gnssdata->OBS[i].sigFlag))
			{
				prn  = gnssdata->OBS[i].PRN-1;
				freq = signal2freq(gnssdata->OBS[i].sigFlag);
				index_map[freq][prn] = i;
			}
			else if (!conf->use_pilot && conf->use_data && isDataSignal(gnssdata->OBS[i].sigFlag))
			{
				prn  = gnssdata->OBS[i].PRN-1;
				freq = signal2freq(gnssdata->OBS[i].sigFlag);
				index_map[freq][prn] = i;
			}
		}
	}

	// it is computed iono-free combination taking a frequency as reference for every combination
	// By default, GPS L1 / GAL E1 is considered the reference frequency
	int prn, freq, ref_freq = FREQ_L1E1;

	for(prn = 0; prn < MAX_N_SAT_G; prn++)
	{
		for(freq = 0; freq < NUM_FREQ; freq++)
		{
			if (freq == ref_freq) { continue; }

			if (index_map[ref_freq][prn] != -1 && index_map[freq][prn] != -1)
			{
				obsdata_t *obsf1 = &gnssdata->OBS[index_map[ref_freq][prn]];
				obsdata_t *obsf2 = &gnssdata->OBS[index_map[freq][prn]];

				// The combined measurement is stored in the second frequency to allow the possible
				// combination of the reference frequency with other ones
				double gamma = (obsf2->lambda * obsf2->lambda) / (obsf1->lambda*obsf1->lambda);
				obsf2->C1 = (obsf2->C1 - gamma * obsf1->C1) / (1 - gamma);
				obsf2->L1 = (obsf2->L1 - gamma * obsf1->L1) / (1 - gamma);
				obsf2->sigFlag = obsf1->sigFlag;
				obsf2->iono_free = 1;

				// Remove the reference observation since already available iono-free
				obsf1->lock = 0;
			}
		}
	}
}

#if USE_IMU == 0
/* This function is used to calculate the traveled distance when the navigation algorithm
 * works in GNSS-only mode.
 */
void calculate_distance_metrics(gicsrx_t *gicsrx, nav_sol_t *nav_sol, double *distance)
{
	static char init_distance = 1;
	if(init_distance == 1 && nav_sol->quality >= KALMAN_ONLY_1)
	{
		gicsrx->exchange.initial_tow      = gicsrx->exchange.tow;
		gicsrx->exchange.initial_distance = 0;
		init_distance = 0;
	}

	// Used for distance computation.
	char prev_state = (nav_sol->prev_used_sats > 0 && nav_sol->kalman_valid == 1),
		 curr_state = (nav_sol->curr_used_sats > 0 && nav_sol->kalman_valid == 1);

	// Update the traveled distance using the Kalman Filter velocity.
	update_distance(nav_sol, nav_sol->prop_time, prev_state, curr_state, distance);

	nav_sol->t_dist_0   = gicsrx->exchange.initial_tow;
	nav_sol->distance_0 = gicsrx->exchange.initial_distance;
	nav_sol->distance_t = *distance;
}
#endif

/* This function is used to store the raw observation data, satellite positions and
 * least-squares navigation solution into the gicsrx_t data structure.
 */
void update_garlic_data(gnssdata_t *gnssdata, lsq_state_t *lsq_state, double geod_pos[3], kconf_t *conf)
{
#if CALC_ATMSPHCORR == 1
	char update_eph = 0;
	int    iyear, imonth, iday, ihour, iminute;
	double day, second;

	GPSToCal_G(gnssdata->week, gnssdata->tow, &iyear, &imonth, &iday, &ihour, &iminute, &second);
	day = DateToDay(iyear, imonth, iday);
#endif

	// Correct clock bias in the pseudo-range measurements to linearize them.
	int i;
	for(i = 0; i < gnssdata->noOfChannelsAv; i++)
	{
		if(gnssdata->OBS[i].prange_status == OBS_VALID)
		{
#if CALC_ATMSPHCORR == 1
			update_eph |= calculate_atmsphCorr(imonth, day, gnssdata->tow, geod_pos, &gnssdata->iono_utc, &gnssdata->OBS[i], conf->iono_model);
			gnssdata->OBS[i].C1 -= (lsq_state->scaled_bias[i] + gnssdata->OBS[i].atmsphCorr);
#endif
		}
	}
#if CALC_ATMSPHCORR == 1
	if(update_eph == 0)
	{
		update_NeQuick_table(gnssdata, &gnssdata->iono_utc.gal, geod_pos, imonth);
	}
#endif
}

int signal2system(sigFlag_t sigFlag)
{
	int sys = -1;
	switch(sigFlag)
	{
	case SIG_GPSL1:
	case SIG_GPSL5I:
	case SIG_GPSL5Q:
		sys = SYS_GPS;
		break;

	case SIG_GALE1B:
	case SIG_GALE1C:
	case SIG_GALE1AD:
	case SIG_GALE1AP:
	case SIG_GALE5aI:
	case SIG_GALE5aQ:
	case SIG_GALE5bI:
	case SIG_GALE5bQ:
	case SIG_GALE6AD:
	case SIG_GALE6AP:
		sys = SYS_GAL;
		break;

	case SIG_GLOL1:
		sys = SYS_GLO;
		break;

	case SIG_BEIB1:
		sys = SYS_BEI;
		break;
	}
	return sys;
}

int signal2freq(sigFlag_t sigFlag)
{
	int freq = -1;
	switch(sigFlag)
	{
	case SIG_GPSL1:
	case SIG_GALE1B:
	case SIG_GALE1C:
	case SIG_GALE1AD:
	case SIG_GALE1AP:
	case SIG_GLOL1:
	case SIG_BEIB1:
		freq = FREQ_L1E1;
		break;

	case SIG_GPSL5I:
	case SIG_GPSL5Q:
		freq = FREQ_L5;
		break;

	case SIG_GALE5aI:
	case SIG_GALE5aQ:
		freq = FREQ_E5a;
		break;

	case SIG_GALE5bI:
	case SIG_GALE5bQ:
		freq = FREQ_E5b;
		break;

	case SIG_GALE6AD:
	case SIG_GALE6AP:
		freq = FREQ_E6;
		break;
	}
	return freq;
}


int sys2garlicIndex(int prnIndex, sysFlag_t sysFlag)
{
	switch(sysFlag)
	{
	case SYS_BEI:
		prnIndex += MAXNUMGLO_G;
	case SYS_GLO:
		prnIndex += MAXNUMGAL_G;
	case SYS_GAL:
		prnIndex += MAXNUMGPS_G;
	case SYS_GPS:
		break;
	default:
		prnIndex = -1;
	}
	return prnIndex;
}

int garlic2sysIndex(int prnIndex, sysFlag_t sysFlag)
{
	switch(sysFlag)
	{
	case SYS_BEI:
		prnIndex -= MAXNUMGLO_G;
	case SYS_GLO:
		prnIndex -= MAXNUMGAL_G;
	case SYS_GAL:
		prnIndex -= MAXNUMGPS_G;
	case SYS_GPS:
		break;
	default:
		prnIndex = -1;
	}
	return prnIndex;
}

char isPilotSignal(sigFlag_t sigFlag)
{
	switch(sigFlag)
	{
	case SIG_GPSL1:
	case SIG_GPSL5Q:
	case SIG_GALE1C:
	case SIG_GALE5aQ:
	case SIG_GALE5bQ:
	case SIG_GALE1AP:
	case SIG_GALE6AP:
	case SIG_GLOL1:
	case SIG_BEIB1:
		return 1;
	default:
		return 0;
	}
}

char isDataSignal(sigFlag_t sigFlag)
{
	switch(sigFlag)
	{
	case SIG_GPSL1:
	case SIG_GPSL5I:
	case SIG_GALE1B:
	case SIG_GALE5aI:
	case SIG_GALE5bI:
	case SIG_GALE1AD:
	case SIG_GALE6AD:
	case SIG_GLOL1:
	case SIG_BEIB1:
		return 1;
	default:
		return 0;
	}
}

double get_wavelength(sigFlag_t sigFlag, int fbMode, int channelSlot)
{
	double lambda = LAMBDA_GPS_L1;

	switch(sigFlag)
	{
	case SIG_GPSL1:
	default:
		break;

	case SIG_GALE1B:
	case SIG_GALE1C:
		lambda /= ((fbMode == 0) ? SBFACTOR_E1BC : 1.);
		break;

	case SIG_GALE1AD:
	case SIG_GALE1AP:
		lambda /= ((fbMode == 0) ? SBFACTOR_E1A : 1.);
		break;

	case SIG_GPSL5I:
	case SIG_GPSL5Q:
	case SIG_GALE5aI:
	case SIG_GALE5aQ:
		lambda = LAMBDA_GPS_L5;
		break;

	case SIG_GALE5bI:
	case SIG_GALE5bQ:
		lambda = LAMBDA_GAL_E5b;
		break;

	case SIG_GALE6AD:
	case SIG_GALE6AP:
		lambda = LAMBDA_GAL_E6;
		break;

	case SIG_GLOL1:
		lambda = LAMBDA_GLO_L1[channelSlot];
		break;

	case SIG_BEIB1:
		lambda = LAMBDA_BEI_B1;
		break;
	}

	return lambda;
}

