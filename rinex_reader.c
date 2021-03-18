/* ATTENTION!!! This version has been used to output the following data:
 * - carrier-phase measurements: file "carrier_phase.txt"
 * - observation matrixes: "obs_matrix.txt"
 * - clock biases: "clk.txt"
 * - satellites positions. "sat_pos.txt"
 * ...
 * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "GICSRxRINEX.h"

#include "matrix.h"
#include "algebra.h"
#include "GICSRxObs.h"
#include "GICSRxPosition.h"
#include "garlic_functions.h"

#include "GICSRxWNavSol.h"
#include "GICSRxKalman.h"

#include "ARAIM_RINEX.h"

#if USE_IMU == 1
#include "GICSRxMechanization.h"
#endif

// Declaration of file pointers
static FILE *fp_obs  = NULL;
static FILE *fp_ngps = NULL;
static FILE *fp_nglo = NULL;
static FILE *fp_ngal = NULL;
static FILE *fp_nbei = NULL;
static FILE *fp_lsq  = NULL;
//static FILE *fp_kal  = NULL;
//static FILE *fp_res  = NULL;
//static FILE *fp_cp = NULL;
//static FILE *fp_obs_matrix = NULL;
//static FILE *fp_obs_matrix_ECEF = NULL;
//static FILE *fp_clk = NULL;
//static FILE *fp_sat_pos = NULL;
//static FILE *fp_true_pos = NULL;
//static FILE *fp_residuals  = NULL;
//static FILE *fp_atm  = NULL;
//static FILE *fp_weight  = NULL;
static FILE *fp_output_all  = NULL;
static FILE *fp_satview = NULL;

FILE *fp_ref = NULL;

// Declaration of internal data structures
//
#if USE_EPHEMERIS_PROP == 1
// Propagated ephemeris support. Stores the node satellite positions
static propephem_t propephem;
#endif

// LSQ state information
static lsq_state_t lsq_state;

// LSQ geometry and residuals information
static lsq_matrix_t lsq_mat;

// Compiled information from the different RINEX sources (*.obs, *.nav)
static rinex_header_t rinex_info;

// Data read from SP3 input file
static SP3STRUCT_D pSp3;

// GNSS observations information
static gnssdata_t gnssdata;

// GNL/garlic algorithms information
static gicsrx_t gicsrx;

#if USE_IMU == 1
// Data block used to contain the IMU samples within a second
int nsamples = 0;
sensor_t sensor_data[MAX_IMU_SAMPLES];
#endif

/* Open output files. If Hybrid Kalman Filter is compiled, the IMU data samples is also open */
static void open_io_files(int argc, char **argv)
{
#if USE_IMU == 1
	// Open IMU file if Hybrid Kalman Filter is compiled
	gicsrx.kconf.fp_imu = fopen(argv[3], "rt");
#endif

	// create a new folder with the name contained in variable argv[4]
	char copyname[300];

	sprintf(copyname, "Output\\%s", argv[2] );
	mkdir(copyname);

	// Open file descriptors for output files (DEBUG). Different file names could
	// be obtained from *argv[] instead of default names lsq_state and kal_state

	sprintf(copyname, "Output\\%s\\lsq_state.txt", argv[2]);
	fp_lsq = fopen(copyname,"wt");

	sprintf(copyname, "Output\\%s\\sat_view.txt", argv[2]);
	fp_satview = fopen(copyname,"wt");

//	sprintf(copyname, "Output\\%s\\kal_state.txt", argv[2]);
//	fp_kal = fopen(copyname,"wt");
//
//	sprintf(copyname, "Output\\%s\\residuals.txt", argv[2]);
//	fp_res = fopen(copyname,"wt");
//
//	// save carrier-phase measurements to file
//	sprintf(copyname, "Output\\%s\\carrier_phase.txt", argv[2]);
//	fp_cp = fopen(copyname,"wt");
//
	// save observation matrix to file
//	sprintf(copyname, "Output\\%s\\obs_matrix.txt", argv[2]);
//	fp_obs_matrix = fopen(copyname,"wt");
//
//	// save observation matrix to file
//	sprintf(copyname, "Output\\%s\\obs_matrix_ECEF.txt", argv[2]);
//	fp_obs_matrix_ECEF = fopen(copyname,"wt");
//
//	// save clock biases
//	sprintf(copyname, "Output\\%s\\clk.txt", argv[2]);
//	fp_clk = fopen(copyname,"wt");
//
//	// save satellites' position
//	sprintf(copyname, "Output\\%s\\sat_pos.txt", argv[2]);
//	fp_sat_pos = fopen(copyname,"wt");
//
//	// save true position
//	sprintf(copyname, "Output\\%s\\true_pos.txt", argv[2]);
//	fp_true_pos = fopen(copyname,"wt");
//
	// save residuals
//	sprintf(copyname, "Output\\%s\\residuals_ibpl.txt", argv[2]);
//	fp_residuals = fopen(copyname,"wt");
//
	// save output_all
	sprintf(copyname, "Output\\%s\\output_all.txt", argv[2]);
	fp_output_all = fopen(copyname,"wt");

//	// save atmospheric correction
//	sprintf(copyname, "Output\\%s\\atmosph_corr.txt", argv[2]);
//	fp_atm = fopen(copyname,"wt");
//
//	// save weigthing matrix
//	sprintf(copyname, "Output\\%s\\weighting.txt", argv[2]);
//	fp_weight = fopen(copyname,"wt");

}

static int getObsIndex(obsdata_t *obslist, int prn)
{
	int k, index;
	for(k = 0, index = -1; index < 0 && k < MAX_CHANNELS_G; k++)
	{
		if(obslist[k].PRN == prn && obslist[k].lock)
		{
			index = k;
		}
	}
	return index;
}

static gnssdata_t gnssdata_prev[3] = {{0}};

/* Close all files. This function is not exactly the opposite from open_io_files, since it also tries
 * to close the RINEX files used during the execution, while the first function only opens the output
 * files and, in case Hybrid Kalman Filter is compiled, the IMU data samples file
 */
static void close_io_files()
{
	if(fp_lsq  != NULL) { fclose(fp_lsq); }
//	if(fp_kal  != NULL) { fclose(fp_kal); }
//	if(fp_res  != NULL) { fclose(fp_res); }

	if(fp_obs  != NULL) { fclose( fp_obs); }
	if(fp_ngps != NULL) { fclose(fp_ngps); }
	if(fp_nglo != NULL) { fclose(fp_nglo); }
	if(fp_ngal != NULL) { fclose(fp_ngal); }

	if(fp_satview  != NULL) { fclose( fp_satview); }

//	if(fp_cp != NULL) { fclose(fp_cp); }
//	if(fp_obs_matrix != NULL) { fclose(fp_obs_matrix); }
//	if(fp_obs_matrix_ECEF != NULL) { fclose(fp_obs_matrix_ECEF); }
//	if(fp_clk != NULL) { fclose(fp_clk); }
//	if(fp_sat_pos != NULL) { fclose(fp_sat_pos); }
//	if(fp_true_pos != NULL) { fclose(fp_true_pos); }
//	if(fp_residuals != NULL) { fclose(fp_residuals); }
//	if(fp_atm != NULL) { fclose(fp_atm); }
//	if(fp_weight != NULL) { fclose(fp_weight); }
//	if(fp_output_all != NULL) { fclose(fp_output_all); }

#if USE_IMU == 1
	if(gicsrx.kconf.fp_imu != NULL) { fclose(gicsrx.kconf.fp_imu); }
#endif
}

/* Write two different output files containing the Least Squares and Kalman Filter solutions */
static void write_output_files(lsq_state_t *lsq_state, gnssdata_t *gnssdata, nav_sol_t *nav_solution)
{
	// Least Squares file
#if CALCULATE_IBPL == 1
//	fprintf(fp_lsq, "%.3f, %d, %d, %d, %d, %d, %d, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n",
//						nav_solution->tow, lsq_state->ibpl_valid, lsq_state->num_iter, lsq_state->num_LOS[SYS_GPS], lsq_state->num_LOS[SYS_GLO], lsq_state->num_LOS[SYS_GAL], lsq_state->num_LOS[SYS_BEI],
//						lsq_state->pos[0], lsq_state->pos[1], lsq_state->pos[2], lsq_state->pos[3], lsq_state->pos[4], lsq_state->pos[5], lsq_state->pos[6], lsq_state->pos[7], lsq_state->pos[8],
//						lsq_state->vel[0], lsq_state->vel[1], lsq_state->vel[2], lsq_state->vel[3],
//						lsq_state->hibpl[0], lsq_state->hibpl[1], lsq_state->hibpl[2], lsq_state->hibpl[3], lsq_state->hibpl[4], lsq_state->hibpl[5], lsq_state->hibpl[6]);

	fprintf(fp_lsq, "%.3f, %d, %d, %d, %d, %d, %d, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n",
						nav_solution->tow, nav_solution->week, lsq_state->ibpl_valid,
						lsq_state->num_LOS[SYS_GPS], lsq_state->num_LOS[SYS_GLO], lsq_state->num_LOS[SYS_GAL], lsq_state->num_LOS[SYS_BEI],
						lsq_state->pos[0], lsq_state->pos[1], lsq_state->pos[2],
						lsq_state->hibpl[0], lsq_state->hibpl[1], lsq_state->hibpl[2], lsq_state->hibpl[3], lsq_state->hibpl[4], lsq_state->hibpl[5], lsq_state->hibpl[6]);
//	// Kalman Filter file
//	fprintf(fp_kal, "%.3f, %d, %d, %d, %d, %d, %d, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n",
//						nav_solution->tow, nav_solution->kalman_valid, nav_solution->kalman_ibpl_valid, gnssdata->noOfChannelsAv,
//						nav_solution->position[0], nav_solution->position[1], nav_solution->position[2], nav_solution->clock_bias, nav_solution->isb_g1,
//						nav_solution->velocity[0], nav_solution->velocity[1], nav_solution->velocity[2], nav_solution->clock_drift,
//						nav_solution->hibpl[0], nav_solution->hibpl[1], nav_solution->hibpl[2], nav_solution->hibpl[3], nav_solution->hibpl[4], nav_solution->hibpl[5], nav_solution->hibpl[6]);
#else
	fprintf(fp_lsq, "%.3f, %d, %d, %d, %d, %d, %d, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n",
						nav_solution->tow, lsq_state->pos_valid, lsq_state->num_iter, lsq_state->num_LOS[SYS_GPS], lsq_state->num_LOS[SYS_GLO], lsq_state->num_LOS[SYS_GAL], lsq_state->num_LOS[SYS_BEI],
						lsq_state->pos[0], lsq_state->pos[1], lsq_state->pos[2], lsq_state->pos[3], lsq_state->pos[4], lsq_state->pos[5], lsq_state->pos[6], lsq_state->pos[7], lsq_state->pos[8], lsq_state->pos[9],
						lsq_state->vel[0], lsq_state->vel[1], lsq_state->vel[2], lsq_state->vel[3], (TWO_RX == 1 ? lsq_state->vel[4] : 0.0));

	// Kalman Filter file
	fprintf(fp_kal, "%.3f, %d, %d, %d, %d, %d, %d, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n",
						nav_solution->tow, nav_solution->kalman_valid, gnssdata->noOfChannelsAv, nav_solution->num_LOS[SYS_GPS], nav_solution->num_LOS[SYS_GLO], nav_solution->num_LOS[SYS_GAL], nav_solution->num_LOS[SYS_BEI],
						nav_solution->position[0], nav_solution->position[1], nav_solution->position[2], nav_solution->clock_bias, nav_solution->isb[0], nav_solution->isb[1], nav_solution->isb[2], nav_solution->ifreqBias[0], nav_solution->ifreqBias[1], nav_solution->ifreqBias[2],
						nav_solution->velocity[0], nav_solution->velocity[1], nav_solution->velocity[2], nav_solution->clock_drift[0], (TWO_RX == 1 ? nav_solution->clock_drift[1] : 0.0));
#endif
	fflush(fp_lsq);
//	fflush(fp_kal);


//	float rsu[3];
//
//	int k, ind;
//	for(k = 0; k < gnssdata->noOfChannelsAv; k++)
//	{
//		double rx_clock_bias;
//
//		// Calculate distance between user and current satellite (PR estimated).
//		double range = GICSRxSatRange(rinex_info.XYZ_antenna, gnssdata->OBS[k].sat.pos, rsu);
//
//		int sys = signal2system(gnssdata->OBS[k].sigFlag);
//		if (sys >= 0)
//		{
//			rx_clock_bias = FSTX_ISBIAS[sys] + lsq_state->pos[3+sys];
//		}
//
//		gnssdata->OBS[k].prange_residual = gnssdata->OBS[k].C1 - range - rx_clock_bias + gnssdata->OBS[k].sat.clk_bias;
//	}
//
//	float median_list_pr[15][MAX_CHANNELS_G];
//	int num[15] = {0};
//	for(k = 0, ind = 0; k < gnssdata->noOfChannelsAv; k++)
//	{
//		if(gnssdata->OBS[k].prange_status == OBS_VALID)
//		{
//			median_list_pr[gnssdata->OBS[k].sigFlag][ind] = gnssdata->OBS[k].prange_residual;
//			ind++;
//			num[gnssdata->OBS[k].sigFlag]++;
//		}
//	}
//
//	float median_val_pr[15];
//	int i;
//	for (i=0; i<15; i++)
//	{
//		compute_median(median_list_pr[i], num[i], &median_val_pr[i]);
//	}
//
//	for(k = 0; k < gnssdata->noOfChannelsAv; k++)
//	{
//		gnssdata->OBS[k].prange_residual  -= median_val_pr[gnssdata->OBS[k].sigFlag];
//	}
//
//	fprintf(fp_res, "%.3f, %.3f, %.3f, ", nav_solution->tow, lsq_state->hdop, lsq_state->vdop);
//	for(k = 0, ind = 0; k < gnssdata->noOfChannelsAv; k++)
//	{
//		if(gnssdata->OBS[k].prange_status == OBS_VALID)
//		{
//			int ind2 = getObsIndex(gnssdata_prev[2].OBS, gnssdata->OBS[k].PRN);
//			int ind1 = getObsIndex(gnssdata_prev[1].OBS, gnssdata->OBS[k].PRN);
//			int ind0 = getObsIndex(gnssdata_prev[0].OBS, gnssdata->OBS[k].PRN);
//
//			int valid_3deriv = ind0 >= 0 && ind1 >= 0 && ind2 >= 0 && (gnssdata_prev[0].OBS[ind0].prange_status == OBS_VALID && gnssdata_prev[1].OBS[ind1].prange_status == OBS_VALID && gnssdata_prev[2].OBS[ind2].prange_status == OBS_VALID);
//			double prange_deriv = (gnssdata->OBS[k].C1 - 3*gnssdata_prev[0].OBS[ind0].C1 + 3*gnssdata_prev[1].OBS[ind1].C1 - gnssdata_prev[2].OBS[ind2].C1)/sqrt(20);
//
//			int prn;
//			// 1:32->GPS_L1 33:68->GAL_E1 69:100->GPS_L5 101:136->GAL_E5a 137:172->GAL_E5b
//			if (signal2system(gnssdata->OBS[k].sigFlag) == SYS_GPS)
//			{
//				prn = signal2freq(gnssdata->OBS[k].sigFlag) == FREQ_L1E1 ? gnssdata->OBS[k].PRN : gnssdata->OBS[k].PRN +32+36;
//			}
//			else if (signal2system(gnssdata->OBS[k].sigFlag) == SYS_GAL)
//			{
//				switch (signal2freq(gnssdata->OBS[k].sigFlag))
//				{
//					case FREQ_L1E1:
//						prn = gnssdata->OBS[k].PRN; break;
//					case FREQ_E5a:
//						prn = gnssdata->OBS[k].PRN +32+36; break;
//					case FREQ_E5b:
//						prn = gnssdata->OBS[k].PRN +32+36+36; break;
//				}
//			}
//
//			fprintf(fp_res, "%2d, %.4f, %.4f, %.4f, %.4f, ", prn , gnssdata->OBS[k].prange_residual,
//					gnssdata->OBS[k].D1,
//					prange_deriv, (double)valid_3deriv);
//			ind++;
//		}
//	}
//	for(k = ind; ind < MAX_CHANNELS_G; ind++)
//	{
//		fprintf(fp_res, "00, 0.0000, 0.0000, 0.0000, 0.0000, ");
//	}
//	fprintf(fp_res, "\n");

	gnssdata_prev[2] = gnssdata_prev[1];
	gnssdata_prev[1] = gnssdata_prev[0];
	gnssdata_prev[0] = *gnssdata;
}

/* This function sets to zero the number of channels available, such that the measurements of the following epoch are properly stored */
static void prepare_next_iteration(gnssdata_t *gnssdata)
{
	// Set to zero the number of channels available
	gnssdata->noOfChannelsAv = 0;

	// Clear satellite mask
	memset(gnssdata->SAT_ID, 0, sizeof(gnssdata->SAT_ID));

	// Reset all the channels
	int k;
	for(k = 0; k < MAX_CHANNELS_G; k++)
	{
		gnssdata->OBS[k].lock      = 0;
		gnssdata->OBS[k].atmsph_av = 0;
	}
}

#if RECEIVER_MODEL == 1
/* During the processing of the NVS data recorded by IFSTTAR, some instabilities in the clock bias behavior has
 * been detected. In particular, the NVS receiver fails when selecting the proper millisecond and the bias drift
 * this amount. This step shall be corrected for an optimal clock bias propagation in the Kalman Filter */
static void handle_nvs_clk_steps(gicsrx_t *gicsrx, gnssdata_t *gnssdata, lsq_state_t *lsq_state)
{
	static double previous_clock_bias = 0.0;

	double bias_correction  = (previous_clock_bias - lsq_state->pos[3]);
	if(previous_clock_bias == 0 || fabs(bias_correction) < (SPEED_OF_LIGHT/1000)/2 ||
	   gicsrx->exchange.tow == gicsrx->exchange.timeOfLastReset)
	{
		bias_correction     = 0;
		previous_clock_bias = lsq_state->pos[3];
	}

	int i;
	for(i = 0; i < gnssdata->noOfChannelsAv; i++)
	{
		gnssdata->OBS[i].C1 += bias_correction;
	}
}
#endif

/* This function is used to read all the navigation files from the input, without knowing the a-priori order of them.
 * Once the files are opened and the header read, we can assign to which flow (GPS, GLONASS, Galileo, mixed) the current
 * file belongs */
static void assign_navfile_pointers(char argc, char **argv)
{
	FILE *fp1 = NULL, *fp2 = NULL, *fp3 = NULL;
	char sat_typ1 = 0, sat_typ2 = 0, sat_typ3 = 0;

	// Open first input navigation file.
	fp1 = fopen(argv[8+USE_IMU], "rt");    // *CHANGE* 3->8)

	// Read the header of the file and get the type of constellation it points to.
	// 0: GPS, 1: GLONASS, 2: Galileo, 3: Mixed
	sat_typ1 = read_nav_header_RINEX(fp1, &gnssdata);

	// Check if the navigation file contains ephemeris information of all constellations (mixed) or not
	if(sat_typ1 < 4)
	{
		// RINEX files including only single constellations. If possible, read more navigation files (default: GLONASS)
		if(argc > (9 + USE_IMU))     // *CHANGE* 4->9
		{
			// Open input navigation file
			fp2 = fopen(argv[9+USE_IMU], "rt");  // *CHANGE* 4->9

			// Read header and get type of constellation. Mixed data should not be got at this point.
			sat_typ2 = read_nav_header_RINEX(fp2, &gnssdata);

			// If possible, read more navigation files (default: Galileo)
			if(argc > (USE_IMU + 10))   // *CHANGE* 5->10
			{
				// Open input navigation file
				fp3 = fopen(argv[10+USE_IMU], "rt");  // *CHANGE* 5->10

				// Read header and get type of constellation. Mixed data should not be got at this point.
				sat_typ3 = read_nav_header_RINEX(fp3, &gnssdata);
			}
		}

		// Assign flows taking into account the type of constellation for each file
		fp_ngps = (sat_typ1 == 0 ? fp1 : (sat_typ2 == 0 ? fp2 : (sat_typ3 == 0 ? fp3 : NULL)));
		fp_nglo = (sat_typ1 == 1 ? fp1 : (sat_typ2 == 1 ? fp2 : (sat_typ3 == 1 ? fp3 : NULL)));
		fp_ngal = (sat_typ1 == 2 ? fp1 : (sat_typ2 == 2 ? fp2 : (sat_typ3 == 2 ? fp3 : NULL)));
		// TODO: Add support for BeiDou only navigation rinex
		fp_nbei = NULL;
	}
	else
	{
		// In mixed mode, all descriptors shall point the same file
		fp_ngps = fp_nglo = fp_ngal = fp_nbei = fp1;
	}
	// Get compiled information of all the RINEX headers read at this point
	get_rinex_info(&rinex_info);
}

/* Assign initial position and velocity to be used in the LSQ algorithm */
static void set_initial_lsq_position(double *init_pos, double *init_vel)
{
	int k;
	// State vector for position and velocity.
	// Origin from KF state vector only if initialized
	if (gicsrx.exchange.kalman_on)
	{
		// Initial position
		memcpy(init_pos, gicsrx.exchange.ins_state_ecef, 3*sizeof(double));
		// Initial clocks
		init_pos[3] = gicsrx.exchange.ins_state_ecef[6];
		for (k = 0; k < NUM_SYS-1; k++) { init_pos[4+k] = gicsrx.exchange.ins_state_ecef[6] + gicsrx.exchange.ins_state_ecef[7+k]; }
		// Initial inter-frequency bias
		memcpy(init_pos+3+NUM_SYS, gicsrx.exchange.ins_state_ecef+6+NUM_SYS, (NUM_FREQ-1)*sizeof(double));

		// Initial velocity
		memcpy(init_vel, gicsrx.exchange.ins_state_ecef+3, 3*sizeof(double));
		// Initial drift
		init_vel[3] =  gicsrx.exchange.ins_state_ecef[6+NUM_SYS+NUM_FREQ-1];
		if (TWO_RX == 1) {init_vel[4] =  gicsrx.exchange.ins_state_ecef[6+NUM_SYS+NUM_FREQ+TWO_RX-1];}
	}
	else
	{
		memcpy(init_pos, lsq_state.pos, sizeof(lsq_state.pos));
		memcpy(init_vel, lsq_state.vel, sizeof(lsq_state.vel));
	}
}

int main(int argc, char **argv)
{

	// keep track of the duration of the simulation
	double test_duration = 0; clock_t begin = clock();

	// Check number of input parameters
	if(argc < (9 + USE_IMU))  // *CHANGE* 4->9
	{
		printf("Invalid number of input parameters.\n");
		return 0;
	}

	// create output folder
	int check_output_folder;
	char* folder_name = "Output";

	check_output_folder = mkdir(folder_name, 0777);

//	int whole_loop_cnt = 0;
//	while( whole_loop_cnt < 2 ){

		// Initialization of all the data structures in garlic: GNSS and algorithm information.
		// This is performed prior to process the navigation and observation files.
#if USE_EPHEMERIS_PROP == 1
		initialize_garlic_data(&gnssdata, &gicsrx, &propephem, &lsq_state, &lsq_mat);
#else
		initialize_garlic_data(&gnssdata, &gicsrx, &lsq_state, &lsq_mat, argv);
#endif

		// Open I/O files
		open_io_files(argc, argv);

		// Get navigation mode: RINEX ephemeris or SP3
		gnssdata.SP3_ON = (argv[1][0] != '0');

		// Open observation RINEX file
		fp_obs = fopen(argv[7], "rt");               // *CHANGE* 2->7!!!!!!!!!!!

		// Read the file header. RINEX v2 and v3 are supported
		read_obs_header_RINEX(fp_obs, &gnssdata, &lsq_state);

		// Check if SP3 mode is activated
		if(gnssdata.SP3_ON == 1)
		{
			// Read the corresponding SP3 orbit file
			GICSRxReadSP3File(&pSp3, argv[3+USE_IMU]);
			// Get compiled information of all the RINEX headers read at this point
			get_rinex_info(&rinex_info);
		}
		else
		{
			// Navigation with RINEX ephemeris
			//
			// Initialize file descriptors for the navigation data
			assign_navfile_pointers(argc, argv);
			// Update GNL/garlic data structures to first available ephemeris
			initialize_ephemeris(fp_ngps, fp_nglo, fp_ngal, fp_nbei, &gnssdata);
		}

		// Ancillary data structure to store the last valid GNSS information read from
		// observation RINEX file. This is useful when epochs in RINEX are missing and
		// state and covariance propagation in the navigation algorithms is required
		gnssdata_t last_gnss = {0};

		// Initial position and velocity state vectors to be used in LSQ
		double init_pos[NUM_PARS], init_vel[4+TWO_RX], geod_pos[3];

		fp_ref = fopen("ref_file.txt", "rt");

		// count the number of epochs
		int n_epochs = 0, max_epochs = 1e9;

		// get number of sat to see if new sats came in in the last epoch
		int sat_iw_previous[MAX_CHANNELS_G] = {0}, sat_iw[MAX_CHANNELS_G] = {0}, new_sat[MAX_CHANNELS_G] = {0};
		int nsat_previous = 0, nsat = 0;
		int sat_new_in, epoch_update_eph;

		// Iterate until EOF of observation file
		while( !feof(fp_obs) && n_epochs < max_epochs )
		{

			// number of epochs
			n_epochs++;

			// Initialization of the output data structure for the current epoch
			nav_sol_t nav_solution = {0};

			// Store last GNSS available data structure and estimate the next TOW
			last_gnss      = gnssdata;
			last_gnss.tow += rinex_info.epochs_interval;

			// Read current observation set from RINEX file. Different functions are used
			// depending on the version of the input file (v2 and v3 are supported)
			if(rinex_info.obs_version == RINEX_v2)
			{
				read_next_obs_RINEX_v2(fp_obs, &gnssdata);
			}
			else
			{
				read_next_obs_RINEX_v3(fp_obs, &gnssdata);
			}

			// Remove observation discarded by configuration choice
			filter_obs_by_config(&gnssdata, &gicsrx.kconf);

			// Combine observations
			combine_obs(&gnssdata, &gicsrx.kconf);

			// check if new sat came in
			memset( new_sat , -1 , MAX_CHANNELS_G*4 );
			sat_new_in = check_new_sat_in( &gnssdata, nsat, &nsat_previous, sat_iw, sat_iw_previous, new_sat );

			// Check if stored ephemeris are valid. If they are out-of-date, try to refresh
			epoch_update_eph = fmod(round(gnssdata.tow), 60) == 0;
			if(gnssdata.SP3_ON == 0 && ( epoch_update_eph || sat_new_in > 0) )
//			if(gnssdata.SP3_ON == 0)
			{
				check_ephemeris_status(fp_ngps, fp_nglo, fp_ngal, fp_nbei, &gnssdata, sat_new_in, epoch_update_eph, new_sat);
			}

			while(!feof(fp_obs))
			{
				// Check whether the time of the week value has been updated to a value close to the
				// expected one (obtained by adding to the previous epoch the INTERVAL field of the
				// observation RINEX header)
				char tow_updated = (last_gnss.tow > gnssdata.tow  + rinex_info.epochs_interval/2 ||
						fabs(last_gnss.tow - gnssdata.tow) < rinex_info.epochs_interval/2);

				// Get pointer to GNSS data structure: if tow is updated as expected, we select the nominal structure.
				// In case tow is not updated, it is possible that epochs are missing in observation RINEX file. Select
				// in this case ancillary structure to properly propagate the navigation solutions
				gnssdata_t *gnss_source = (tow_updated ? &gnssdata : &last_gnss);

				// Set initial position for LSQ (in steady state, the last position computed by Kalman Filter)
				set_initial_lsq_position(init_pos, init_vel);

				double previous_clock = lsq_state.pos[3];


				// Call Least Squares algorithm
#if USE_EPHEMERIS_PROP == 1
				GICSRxWNavSol(gnss_source, gnssdata.SP3_ON ? &pSp3 : NULL, &lsq_state, &lsq_mat, &gicsrx.kconf, &propephem, &nav_solution, init_pos, init_vel);
#else
				GICSRxWNavSol(gnss_source, gnssdata.SP3_ON ? &pSp3 : NULL, &lsq_state, &lsq_mat, &gicsrx.kconf, &nav_solution, init_pos, init_vel);
#endif


				if(lsq_state.num_LOS[0] != 0){
					gnss_source->tow -= lsq_state.pos[3]/SPEED_OF_LIGHT;
				}
				else if(lsq_state.num_LOS[1] != 0){
					gnss_source->tow -= lsq_state.pos[4]/SPEED_OF_LIGHT;
				}
				else if(lsq_state.num_LOS[2] != 0){
					gnss_source->tow -= lsq_state.pos[5]/SPEED_OF_LIGHT;
				}
				else if(lsq_state.num_LOS[3] != 0){
					gnss_source->tow -= lsq_state.pos[6]/SPEED_OF_LIGHT;
				}


				/*-------------------------------------------------------------------------------------------------------------------------------------------------
	RUN ARAIM USER ALGORITHM
-------------------------------------------------------------------------------------------------------------------------------------------------*/
#if ARAIM_ALG > 0

				// 3) Import data from RINEX structures
				configuration ARAIM; load_ARAIM_parameters( &ARAIM );
				int Nsat = lsq_state.total_obs, i, Nconst = 0; for( i = 0; i < NUM_SYS; i++ ) Nconst += ( lsq_state.num_LOS[i] != 0 );
				int Nx = 3 + Nconst, lsq_num_LOS[NUM_SYS], PRN_used[Nsat]; double G[Nsat][Nx];
				double w_ura[Nsat][Nsat], w_ure[Nsat][Nsat], c_ure[Nsat][Nsat], pr_res[Nsat];
				from_rinex_structures( &lsq_state, &lsq_mat, &gnssdata, &ARAIM, Nsat, Nx, lsq_num_LOS, PRN_used, &G[0][0], &w_ura[0][0], &w_ure[0][0], &c_ure[0][0], pr_res);


				// --------------------------------- print information Miguel ---------------------------------------
				fprintf( fp_output_all, "$POS, %.5f, %d, %d, %.3f, %.3f, %.3f\n",
						nav_solution.tow, nav_solution.week, lsq_state.ibpl_valid, lsq_state.pos[0], lsq_state.pos[1], lsq_state.pos[2]);

				int isat, iobs;
				for( isat = 0; isat < Nsat; isat++ ){
					fprintf( fp_output_all, "$SAT, %.5f, %d, %d, %.5f, %.5f, %.5f, ", nav_solution.tow, nav_solution.week, PRN_used[isat], G[isat][0], G[isat][1], G[isat][2]);

					for( iobs = 0; iobs < gnss_source->noOfChannelsAv; iobs++ ){
						if( gnss_source->OBS[iobs].PRN == PRN_used[isat] ){
							fprintf( fp_output_all, "%.5f, %.5f\n", gnss_source->OBS[iobs].prange_residual, gnss_source->OBS[iobs].S1);
						}
					}

				}

				// --------------------------------- sav sat in view ---------------------------------------
				fprintf( fp_satview, "%.5f ", nav_solution.tow);
				for( isat = 0; isat < Nsat; isat++ ){
					fprintf( fp_satview, "%d ", PRN_used[isat]);
				}
				for( isat = Nsat; isat < 40; isat++ ){
					fprintf( fp_satview, "%d ", 0);
				}
				fprintf( fp_satview, "\n");


////				// --------------------------------- save to file observation matrix (ENU) ---------------------------------------
//				int i_obs, j_obs;
//				printf( "%d, %.5f, %d, %d, %d, %d \n" , Nsat, nav_solution.tow, 0, 0, 0, 0 );
//				for( i_obs = 0; i_obs < Nsat; i_obs++ ){
//					if( PRN_used[i_obs] <= MAXNUMGPS_G ){
//						printf( "%d, ", PRN_used[i_obs] );
//						for( j_obs = 0; j_obs < Nx; j_obs++ ){
//							printf( "%.5f ," , G[i_obs][j_obs] );
//						}
//						printf( "\n" );
//					}
//				}
//				for( i_obs = 0; i_obs < Nsat; i_obs++ ){
//					if( PRN_used[i_obs] > MAXNUMGPS_G && PRN_used[i_obs] <= MAXNUMGPS_G + MAXNUMGAL_G ){
//						printf( "%d, ", PRN_used[i_obs] );
//						for( j_obs = 0; j_obs < Nx; j_obs++ ){
//							printf( "%.5f ," , G[i_obs][j_obs] );
//						}
//						printf( "\n" );
//					}
//				}
//				for( i_obs = 0; i_obs < Nsat; i_obs++ ){
//					if( PRN_used[i_obs] > MAXNUMGPS_G + MAXNUMGAL_G && PRN_used[i_obs] <= MAXNUMGPS_G + MAXNUMGAL_G + MAXNUMGLO_G ){
//						printf( "%d, ", PRN_used[i_obs] );
//						for( j_obs = 0; j_obs < Nx; j_obs++ ){
//							printf( "%.5f ," , G[i_obs][j_obs] );
//						}
//						printf( "\n" );
//					}
//				}
//				for( i_obs = 0; i_obs < Nsat; i_obs++ ){
//					if( PRN_used[i_obs] > MAXNUMGPS_G + MAXNUMGAL_G + MAXNUMGLO_G && PRN_used[i_obs] ){
//						printf( "%d, ", PRN_used[i_obs] );
//						for( j_obs = 0; j_obs < Nx; j_obs++ ){
//							printf( "%.5f ," , G[i_obs][j_obs] );
//						}
//						printf( "\n" );
//					}
//				}
//				printf( "\n\n");
//
//				for( i_obs = 0; i_obs < Nsat; i_obs++ ){
//					printf( "%d, ", PRN_used[i_obs] );
//					for( j_obs = 0; j_obs < Nx; j_obs++ ){
//						printf( "%.5f ," , G[i_obs][j_obs] );
//					}
//					printf( "\n" );
//				}
//				printf( "\n\n");
//				//
//				exit(0);


//				// --------------------------------- save to file observation matrix (ENU) ---------------------------------------
//				int i_obs, j_obs;
//				fprintf( fp_obs_matrix , "%d, %.5f, %d, %d, %d, %d \n" , Nsat, nav_solution.tow, 0, 0, 0, 0 );
//				for( i_obs = 0; i_obs < Nsat; i_obs++ ){
//					fprintf( fp_obs_matrix , "%d, ", PRN_used[i_obs] );
//					for( j_obs = 0; j_obs < Nx; j_obs++ ){
//						fprintf( fp_obs_matrix, "%.5f ," , G[i_obs][j_obs] );
//					}
//					fprintf( fp_obs_matrix , "\n" );
//				}
//				fprintf( fp_obs_matrix, "\n\n");

				// --------------------------------- save to file carrier-phase measurements ---------------------------------------
//				if( gnss_source->tow >= 286000 && gnss_source->tow <= 286260 ){
//					printf( "%f, %d\n" , gnss_source->tow, gnss_source->noOfChannelsAv);
//					for( isat = 0; isat < gnss_source->noOfChannelsAv ; isat ++ ){
//						printf( "%10.5f, %d\n", gnss_source->OBS[isat].C1 , gnss_source->OBS[isat].PRN);
//					}
//					printf( "\n\n" );
//				}

				// --------------------------------- save to file weighting matrix ---------------------------------------

				//			fprintf( fp_weight , "%d, %d\n" , Nsat , 0 );
				//			for( isat = 0; isat < Nsat ; isat ++ ){
				//				fprintf( fp_weight , "%d, %10.5f\n", PRN_used[isat] , lsq_mat.Wpos[isat] );
				//			}
				//			fprintf( fp_weight , "\n\n" );

				//// --------------------------------- save to file LS residuals ---------------------------------------
//				fprintf( fp_residuals , "%d, %.5f\n" , gnss_source->noOfChannelsAv , nav_solution.tow );
//				for( isat = 0; isat < gnss_source->noOfChannelsAv ; isat ++ ){
//					fprintf( fp_residuals , "%d, %10.5f\n", gnss_source->OBS[isat].PRN, gnss_source->OBS[isat].prange_residual );
//				}
//				fprintf( fp_residuals , "\n\n" );

				// --------------------------------- save to file clock biases ---------------------------------------

				//			fprintf( fp_clk , "%d, %d, %d, %d\n" , gnss_source->noOfChannelsAv , 0 , 0 , 0 );
				//			for( isat = 0; isat < gnss_source->noOfChannelsAv ; isat ++ ){
				//				fprintf( fp_clk , "%d, %10.5f, %10.5f, %10.5f\n", gnss_source->OBS[isat].PRN,
				//						gnss_source->OBS[isat].sat.clk_bias , gnss_source->OBS[isat].sat.clk_drift, gnss_source->OBS[isat].sat.clk_driftrate );
				//			}
				//			fprintf( fp_clk , "\n\n" );

				// --------------------------------- save to file atmospheric corrections ---------------------------------------

				//			fprintf( fp_atm , "%d, %d\n" , gnss_source->noOfChannelsAv , 0 );
				//			for( isat = 0; isat < gnss_source->noOfChannelsAv ; isat ++ ){
				//				fprintf( fp_atm , "%d, %10.5f\n", gnss_source->OBS[isat].PRN, gnss_source->OBS[isat].atmsphCorr );
				//			}
				//			fprintf( fp_atm , "\n\n" );

				// --------------------------------- save observation matrix ECEF ---------------------------------------


				//			fprintf( fp_obs_matrix_ECEF , "%d, %d, %d, %d, %d, %d \n" , Nsat, 0, 0, 0, 0, 0 );
				//			for( i_obs = 0; i_obs < Nsat; i_obs++ ){
				//				fprintf( fp_obs_matrix_ECEF , "%d, ", PRN_used[i_obs] );
				//				for( j_obs = 0; j_obs < Nx; j_obs++ ){
				//					fprintf( fp_obs_matrix_ECEF, "%.5f ," , lsq_mat.G[i_obs][j_obs] );
				//				}
				//				fprintf( fp_obs_matrix_ECEF , "\n" );
				//			}
				//			fprintf( fp_obs_matrix_ECEF, "\n\n");

				// --------------------------------- save to file satellites position ---------------------------------------

				//			fprintf( fp_sat_pos , "%d, %d, %d, %d\n" , gnss_source->noOfChannelsAv , 0 , 0 , 0 );
				//			for( isat = 0; isat < gnss_source->noOfChannelsAv ; isat ++ ){
				//				fprintf( fp_sat_pos , "%d, %10.5f, %10.5f, %10.5f\n", gnss_source->OBS[isat].PRN,
				//						gnss_source->OBS[isat].sat.pos[0] , gnss_source->OBS[isat].sat.pos[1], gnss_source->OBS[isat].sat.pos[2] );
				//			}
				//			fprintf( fp_sat_pos , "\n\n" );

				// --------------------------------- save to file true position ------------------------------------------------------------------------------

				//			fprintf( fp_true_pos , "%10.5f, %10.5f, %10.5f\n", lsq_state.pos[0] , lsq_state.pos[1], lsq_state.pos[2] );

				// -------------------------------------------------------------------------------------------------------------------------------------


#endif
				/*-------------------------------------------------------------------------------------------------------------------------------------------------
	END OF ARAIM USER ALGORITHM
-------------------------------------------------------------------------------------------------------------------------------------------------*/

				// If IMU is being used (post-processing mode), read the set of samples for the current epoch from the
				// selected input file.
#if USE_IMU == 1
				read_sensor_file(gnss_source->tow, sensor_data, &nsamples, &gicsrx.kconf);
#endif
				// Calculate atmospheric corrections. For the tropospheric term, the Neill model is proposed.
				// For the ionospheric term, NeQuick is used if available by default. On the contrary, Klobuchar model is set
				ECEFtoNAV_pos(init_pos, geod_pos);
				update_garlic_data(gnss_source, &lsq_state, geod_pos, &gicsrx.kconf);

#if RECEIVER_MODEL == 1
				// Check clock bias steps if NVS receiver is selected
				handle_nvs_clk_steps(&gicsrx, gnss_source, &lsq_state);
#endif

				if(fabs(lsq_state.pos[3] - previous_clock) > SPEED_OF_LIGHT/2000.0)
				{
					gicsrx.exchange.ins_state_ecef[6] += (lsq_state.pos[3] > previous_clock ? +1 : -1)*SPEED_OF_LIGHT/1000.0;
				}

				// Call Kalman Filter algorithm
#if USE_IMU == 1
				kalman_hybrid(&gicsrx, gnss_source, &lsq_state, sensor_data, nsamples, &nav_solution);
#else
				kalman_gnss(&gicsrx, gnss_source, lsq_state.pos_valid, &lsq_state, &nav_solution);
#endif
				// Dump the two available navigation solutions
				write_output_files(&lsq_state, gnss_source, &nav_solution);


				// If current TOW is equal to the last available epoch in RINEX observation file, go to the
				// next one
				if(tow_updated == 1)
				{
					// Reset channels
					prepare_next_iteration(gnss_source);
					break;
				}
				else
				{
					// If there exist a gap between epochs in RINEX file, propagate time in the ancillary data structure
					gnss_source->tow += rinex_info.epochs_interval;
				}
			}

		}

		// Close all open files
		close_io_files();
		fclose(fp_ref);

		/* keep track of the duration of the simulation */
		clock_t end = clock(); test_duration = (double)( end - begin ) / CLOCKS_PER_SEC;

		printf("RINEX file correctly read in %.2f seconds for %d epochs!", test_duration, n_epochs );
//		whole_loop_cnt++;
//
//	}

	return 0;
}
