//////////////////////////////////////////////////////////////////////////////
// Project:                     GICSRX
// Purpose:                     Standard functions for GICSRx algorithms
// File:                        GICSRx_defines.h
// Version:                     1.0
// Author:                      GMV, S.A.
// Language:                    c
// Compiler:                    ANSI
// Creation date:               -
// Last change date:            -
// Changes:
// +--------------------------------------------------------------------------+
// | Version |   Date   | Author       | Change reason                        |
// |---------+----------+--------------+--------------------------------------|
// |   1.00  | 13/12/12 | cmvv         | Initial Release                      |
// +--------------------------------------------------------------------------+
//
// DESCRIPTION:
// Definition of basic GICSRx data structure and general settings for GICSRx
// algorithm.
//
///////////////////////////////////////////////////////////////////////////////
#ifndef GICSRx_definesH
#define GICSRx_definesH

#include "dbdefine.h"

#include <stdio.h>

#define RECEIVER_MODEL 		   3 	// 0: Teseo-II, 1: NVS, 2: Other
#define ARAIM_ALG			   1
#define USE_IMU                0	// Use inertial navigation or GNSS-only.
#define USE_EPHEMERIS_PROP     0	// Use ephemeris propagation.
#define USE_WEIGHTING_MODEL    0	// Use weighting model.
#define CHECK_OUTLIERS         0	// Check faulty measurements in LSQ solution.
#define CALCULATE_IBPL         1	// Calculate IBPL.
#define CALC_ATMSPHCORR	       1	// Calculate iono and tropo corrections instead of using Teseo-II corrections.

#if USE_IMU == 1

#define CALC_INSTALL_MAT		0
#define CONSTANT_GRAVITY    	1	// Assume that the gravity vector is constant throughout a second.
#define CORRECT_ANGLE       	0	// Correct earth rotation term in gyro data.
#define CORRECT_ACC         	0	// Correct Coriolis and centripetal terms in the gravity vector.
#define MECHANIZED_VEL_TYPE 	1	// Select different types of mechanizing the vehicle position (0,1,2)
#define USE_FAST_MECHANIZATION  1   // Use CPU optimized (but ad hoc) method to perform IMU mechanization
#define UPDATE_ATT_MODE     	1	// Select update mode for the attitude matrix.
#define MIN_IMU_SAMPLES			95
#define MAX_IMU_SAMPLES     	105	// Maximum number of samples in IMU matrix for a single epoch.
#define CALC_GNSS_ONLY          1   // Calculate a GNSS-only solution in parallel to the GNSS+IMU solution.
#define USE_ANTISPOOFING		0

#if USE_ANTISPOOFING == 1
// Anti-spoofing
#define AS_WSIZE 			  	5
#define AS_LAMBDA 			  	0.995
#define AS_EFFWIN				(1/(1-AS_LAMBDA))
#define AS_THRESHOLD			0.10
#endif

#else
#define CALC_GNSS_ONLY         0   // Calculate a GNSS-only solution in parallel to the GNSS+IMU solution.
#endif

#if RECEIVER_MODEL == 0

static const double FSTX_ISBIAS[NUM_SYS] = {220, (-0.08*SPEED_OF_LIGHT), (-0.08*SPEED_OF_LIGHT), 0};

// Definition of central frequencies used in Fastrax/STM for the code loop.
#define FSTXFREQ_GPS 	-47122.395
#define FSTXFREQ_GLO 	-47897.569
#define FSTXFREQ_GAL 	-47122.395

static const double drift_values[3] = {-(FSTXFREQ_GPS)*LAMBDA_GPS_L1 - SPEED_OF_LIGHT/1000,
								 	   -(FSTXFREQ_GPS)*LAMBDA_GPS_L1,
								 	   -(FSTXFREQ_GPS)*LAMBDA_GPS_L1 + SPEED_OF_LIGHT/1000};
#else

static const double FSTX_ISBIAS[NUM_SYS] = {0, 0, 0, 0} ;

#endif

// Global settings
#define LEAPSECONDS 18			// Current number of leap seconds.

#define NUM_PARS    (4 + (NUM_SYS-1) + (NUM_FREQ-1)) 	// Size of state vector [x,y,z,clk] (+ isb per system + ifreqbias per freq band)

#if CALCULATE_IBPL == 1
#define NCONF_LVL   7 	// No of confidence level 10-1...10-7
#endif

#if USE_IMU == 0
#define KNSTATES 			(8 + (NUM_SYS-1) + (NUM_FREQ-1) + TWO_RX)
#else
#define KNSTATES 			(17 + (NUM_SYS-1) + (NUM_FREQ-1) + TWO_RX)
#endif

#define VEL_THRESHOLD		0.2     // Velocity threshold used in GNSS odometer.
#define USE_CHI2TEST        1		// Activate measurement rejection.
#define CLK_MAXITER 		5

#define NORMAL_MATRIX_SIZE ( ( (NUM_PARS + 1) * NUM_PARS) / 2)	// Size of covariance matrix in LSQ.
#define COV_MATRIX_SIZE    ( ( (KNSTATES + 1) * KNSTATES) / 2)	// Size of covariance matrix in KF.

#if CALC_GNSS_ONLY == 1
#define GNSS_KSTATES   (8 + (NUM_SYS-1) + (NUM_FREQ-1) + TWO_RX)
#define GNSS_COV_SIZE (((GNSS_KSTATES + 1) * GNSS_KSTATES) / 2)
#endif

#if USE_IMU == 0 || CALC_GNSS_ONLY == 1
#define USE_CONSTANT_HEIGHT 1
#endif

#if USE_WEIGHTING_MODEL == 1
// Model obtained with dynamic data (I)
static const float sigmas_PR[42] = {51.92f, 46.37f, 42.49f, 33.94f, 30.92f, 26.21f, 25.23f,
									22.20f, 21.26f, 19.85f, 17.45f, 15.73f, 14.61f, 13.16f,
									12.68f, 12.04f, 11.26f, 10.52f, 10.24f,  9.82f,  9.50f,
									 9.36f,  9.25f,  9.06f,  8.84f,  8.60f,  8.51f,  8.41f,
									 7.76f,  7.32f,  7.08f,  6.29f,  5.22f,  4.95f,  5.01f,
									 4.66f,  3.74f,  2.68f,  1.56f,  1.56f,  1.56f,  1.56f};

static const float sigmas_DP[42] = {5.12f, 4.32f, 4.46f, 3.49f, 3.25f, 2.62f, 2.68f,
									2.55f, 2.45f, 2.45f, 2.34f, 2.06f, 1.69f, 1.50f,
									1.28f, 1.09f, 0.91f, 0.75f, 0.62f, 0.48f, 0.39f,
									0.29f, 0.20f, 0.15f, 0.14f, 0.12f, 0.12f, 0.11f,
									0.10f, 0.10f, 0.10f, 0.11f, 0.10f, 0.10f, 0.09f,
									0.09f, 0.09f, 0.09f, 0.08f, 0.08f, 0.08f, 0.08f};
#endif

typedef struct
{
	// LSQ position and velocity state vectors.
	double pos   [NUM_PARS];
	double vel   [4+TWO_RX];

	// Scale bias factor in raw pseudo-range measurements.
	double scaled_bias[MAX_CHANNELS_G];

	int    num_iter;	    // Number of iterations required to converge in LSQ.
	int    total_obs;	    // Total number of measurements used.
	int    total_obs_prev;	// Total number of measurements used in previous epoch.

	int    num_LOS[NUM_SYS];    // Number of satellites used per system.
	int	   total_LOS;           // Total number of satellites in view and used.
	int    used_freq[NUM_FREQ]; // Number of measurements per frequency band.

	int    prn_used[MAX_CHANNELS_G];	// PRN used mask.

	char   pos_valid;	// Position validity flag.

	float  prange_residuals;
	float  doppler_residuals;
	int    num_pos_dof;
	int    num_vel_dof;

#if CALCULATE_IBPL == 1
	// IBPL validity flag.
	char ibpl_valid;
	// Horizontal IBPL for different confidence levels.
	float hibpl[NCONF_LVL];
#endif

	// DOP matrices for position and velocity in LSQ.
	float DOPmat_pos[NORMAL_MATRIX_SIZE];
	float DOPmat_vel[NORMAL_MATRIX_SIZE];
	float hdop;
	float vdop;
	float hdop_vel;
	float vdop_vel;

	int num_pos_pars;
	int num_vel_pars;

} lsq_state_t;

typedef struct
{
	float G   [MAX_CHANNELS_G][NUM_PARS];	// Geometry matrix.
	float zpos[MAX_CHANNELS_G];			    // Residuals of position.
	float zvel[MAX_CHANNELS_G];			    // Residuals of velocity.
	float Wpos[MAX_CHANNELS_G];			    // Weighting matrix for position.
	float Wvel[MAX_CHANNELS_G];			    // Weighting matrix for velocity.

} lsq_matrix_t;

#if CALCULATE_IBPL == 1
typedef struct
{
	float eff_dof;
	float r;
	float R[2][2];
	float dop_doppler;
} epoch_kibpl_parameters_t;

typedef struct
{
	float acc_dof;
	float R[2][2];
	float u[2];
} accumulated_kibpl_parameters_t;

// GICSRx KIBPL Data Structure
typedef struct
{
	char   ibpl_epoch_flag;
	double last_ibpl_tow;
	epoch_kibpl_parameters_t prange_epoch;
	epoch_kibpl_parameters_t doppler_epoch;
	float U[2][2];
	float acc_norm_imu;
	float acc_indicator;
	unsigned int num_acc_samples;
} kibpl_variables_t;

typedef struct
{
	float accumulated_r_prange;
	float accumulated_r_doppler;
	float accumulated_r_vel;
	float acc_dof_prange;
	float acc_dof_doppler;
	float acc_dof_vel;

} kipbl_parameters_t;

typedef struct
{
	kibpl_variables_t common;
	accumulated_kibpl_parameters_t prange_parameters;
	accumulated_kibpl_parameters_t doppler_parameters;
} kibpl_t;

#endif

// GICSRx Exchange Data Structure
typedef struct
{
	int    week;		 // GPS Week.
	double tow;          // Time of Week

	double timeOfLastReset;

	// Kalman filter data.
	double ins_state_ecef [KNSTATES];
	double kal_state_ecef [KNSTATES];

#if USE_IMU == 1
	double ins_state_nav  [KNSTATES];
	float quaternion[4];
#endif

	float kal_cov_mat[COV_MATRIX_SIZE];

	// Number of used measurements in Kalman filter.
	int kal_prev_sats;
	int kal_used_sats;

	int num_LOS[NUM_SYS];
	int kal_prn_used[MAX_CHANNELS_G];

	// Propagation time in Kalman filter.
	float prop_time;

	// Kalman filter status.
	char kalman_on;

	// Distance travelled by the vehicle.
	double initial_tow;
	double initial_distance;
	double travelled_distance;

#if CALCULATE_IBPL == 1
	kibpl_t kibpl;
	kibpl_t kibpl_gnss_only;
#endif

#if CALC_GNSS_ONLY == 1
	double gnss_state   [GNSS_KSTATES];
	double gnss_delta   [GNSS_KSTATES];
	float  gnss_cov_mat [GNSS_COV_SIZE];

	char   kalman2_on;
	int    gnss_prev_sats;
	int    gnss_used_sats;

	int    gnss_gps_sats;
	int    gnss_glo_sats;
	int    gnss_gal_sats;

#endif

	int gnss_frequency_counter;

} exchange_t;

typedef struct
{
	float  mask_angle;		// Elevation mask (degrees).
	float  CN0_thr;			// Minimum threshold for carrier to noise ratio in measurements.

	char   use_pilot;
	char   use_data;
	char   use_sys[NUM_SYS];
	char   use_freq[NUM_FREQ];

	ionoModel_t iono_model;

	int    win_PR; 			// Window size for the observation of PR rejection.
	int    tol_PR; 			// Minimum accepted epochs in window.
	float  per_PR; 			// Minimum PRs accepted threshold for pseudoranges in a epoch.
	float  hkvel0;			// Minimum vehicle speed required to initialize Hybrid GNSS+IMU KF.
	int    noSats_min;		// Minimum number of satellites required to initialize GNSS+IMU KF.
	float  chi_test1;		// Chi-square test threshold for Hybrid GNSS+IMU KF.
#if CALC_GNSS_ONLY == 1
	float  chi_test2;		// Chi-square test threshold for GNSS-only KF.
#endif
	int    maxKalx;			// Maximum number of epochs for no updates.

	int    GNSS_rate;

#if USE_WEIGHTING_MODEL == 0
	float  sigma_prange;	// Variance for pseudorange measurements.
	float  sigma_doppler;	// Variance for doppler measurements.
#endif

#if (USE_IMU == 0 || CALC_GNSS_ONLY == 1)
	float  sigma2_vel;		// Variance for process noise in speed.
#endif
	float sigma2_clk;  				// Variance for process noise in clock drift.
	float sigma2_isb[NUM_SYS-1];	// Variance for process noise in isb
	float sigma2_ifreqb[NUM_FREQ-1];// Variance for process noise in interfrequency bias

	float ktimeNoSmooth;
	float kvelStopThr;		// Threshold velocity to consider the vehicle moving: used for weighting ranges.
	float kvelStopFactor;	// Factor used to smooth the effect of the pseudo-range measurements while stopped.

#if USE_IMU == 1
	float sigma2_bacc;		// Variance of accelerometer bias in Gauss-Markov.
	float sigma2_bgyr;		// Variance of accelerometer for random noise in Gauss-Markov.

	float sigma2_wnacc;		// Variance of gyro bias in Gauss-Markov.
	float sigma2_wngyr;		// Variance of gyro for random noise in Gauss-Markov.

	float tau_acc;			// Time constant in Gauss-Markov for accelerometer.
	float tau_gyr;			// Time constant in Gauss-Markov for gyro.

	int   Imat_ready;
	float Imat[3][3];		// Installation matrix.

	int IMU_rate;

	FILE *fp_imu;
#endif

#if CALCULATE_IBPL == 1

	// IBPL: SUM OF PSEUDORANGE AND DOPPLER ERRORS
	char errors_absolute_sum;

	// IBPL: CORRELATIONS
	float gamma_corr_prange;
	float gamma_corr_doppler;

	// IBPL: DOF COMPUTATION
	float dof_tail_depth;
	float max_dof;

	// IBPL: EPOCH DOP FACTOR;
	float prange_dop_factor;
	float doppler_dop_factor;

	// IBPL: ACCELERATION INDICATOR
	float gamma_acc_indicator;
	float acc_indicator_factor_gnss;
	float acc_indicator_factor_imu;
#endif

} kconf_t;

typedef struct
{
	exchange_t exchange;	// Algorithms information.
	kconf_t    kconf;		// Algorithms configuration.

} gicsrx_t;

typedef enum
{
	NOSOL         = 0,
	WLSQ_SOL      = 1,
	KALMAN_ONLY_0 = 2,
	KALMAN_HIMU_0 = 3,
	KALMAN_ONLY_1 = 4,
	KALMAN_HIMU_1 = 5

} quality_t;

typedef struct
{
	int    week;
	double tow;

	char kalman_valid;
	double position[3];
	double velocity[3];
	double clock_bias;
	double clock_drift[1+TWO_RX];
	double isb[NUM_SYS-1];
	double ifreqBias[NUM_FREQ-1];

	int prev_used_sats;
	int curr_used_sats;

	int num_LOS[NUM_SYS];
	int prn_used[MAX_CHANNELS_G];

	float prop_time;

	double t_dist_0;
	double distance_0;
	double distance_t;

	quality_t quality;

	float hdop;
	float vdop;

#if CALCULATE_IBPL == 1
	char kalman_ibpl_valid;
	float hibpl[NCONF_LVL];	// Isotropy-Based HPL
#endif

	double raw_imu[3];

#if CALC_GNSS_ONLY == 1
	char   discr_alarm;
#endif
#if USE_ANTISPOOFING == 1
	char   spoof_alarm;
	double spoof_var;
#endif
} nav_sol_t;

#if USE_IMU == 1
typedef struct
{
	double gps_time;		// GPS time of the sensor measurement.
	double cpu_time;		// CPU time of the sensor measurement.

	float delta_t;			// Elapsed time from the last sensor measurement.

	double  acc_b[3];		// Accelerometer value in current time.
	double  gyr_b[3];		// Gyro value in current time.

} sensor_t;

#endif

typedef enum
{
	KIBPL_OK           = 0,
	KIBPL_NOT_VALID    = 1,
	KIBPL_RESET_FILTER = 2

} kibpl_code_t;

#endif
