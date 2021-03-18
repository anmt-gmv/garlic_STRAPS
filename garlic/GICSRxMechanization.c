#include <math.h>
#include <string.h>

#include "GICSRxMechanization.h"

#include "matrix.h"
#include "algebra.h"

#include "GICSRxObs.h"
#include "GICSRxAux.h"
#include "GICSRxKalman.h"
#include "GICSRxKalman_gnss.h"

#if USE_IMU == 1
void read_sensor_file(double tow, sensor_t imu[MAX_IMU_SAMPLES], int *nsamples, kconf_t *kconf)
{
	int k = 0;

	*nsamples = 0;

	static sensor_t sample = {0};
	static double last_gpstime = 0;

	if(sample.gps_time == 0 || sample.gps_time < tow)
	{
		double acc_b[3], gyr_b[3];
		while(!feof(kconf->fp_imu))
		{
			fscanf(kconf->fp_imu,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &sample.gps_time, &sample.cpu_time,
				   &acc_b[0], &acc_b[1], &acc_b[2], &gyr_b[0], &gyr_b[1], &gyr_b[2]);

			matVecMul_f(sample.acc_b, kconf->Imat, acc_b, 3, 3);
			matVecMul_f(sample.gyr_b, kconf->Imat, gyr_b, 3, 3);

			if(sample.gps_time > (tow - 1.000) && sample.gps_time < (tow + 0.005))
			{
				if(k == 0)
				{
					if(last_gpstime == 0)
					{
						sample.delta_t = 0.01;
					}
					else
					{
						sample.delta_t = sample.gps_time - last_gpstime;
					}
				}
				else
				{
					sample.delta_t = (sample.gps_time - imu[k-1].gps_time);
				}

				memcpy(&imu[k], &sample, sizeof(sensor_t));
				k++;
			}

			if(sample.gps_time >= (tow - 0.005))
			{
				break;
			}

			if(k == MAX_IMU_SAMPLES)
			{
				break;
			}
		}
		*nsamples = k;
	}
	last_gpstime = (k > 0 ? imu[k-1].gps_time : tow);
}

static void propagate_kalimu_deterministic(gicsrx_t *gicsrx, const float acc_b[3], float delta_t, float acc_e [3])
{
	// Calculate different factors that will be used in the function repeatedly.
	float invtau_acc_dt  = 1/gicsrx->kconf.tau_acc*delta_t;
	float invtau_gyr_dt  = 1/gicsrx->kconf.tau_gyr*delta_t;

	float Cbn   [3][3];
	float Cne   [3][3];
	float Cbe   [3][3];

	// Obtain attitude matrix.
	quaternion2mat(gicsrx->exchange.quaternion, Cbn);

	NAVtoECEF_mat(Cne, gicsrx->exchange.ins_state_nav[0], gicsrx->exchange.ins_state_nav[1]);

	// Calculate the rotation matrix from body to ECEF coordinate frame.
	matMul_f(Cbe, Cne, Cbn, 3, 3, 3);

	// Calculate the corrected acceleration in ECEF coordinates.
	matVecMul_f(acc_e, Cbe, acc_b, 3, 3);

	// The transition matrix is expressed as a sparse matrix, since it includes many 0 elements.
	// The use of sparse matrix helps to speed up the mechanization progress, which is crucial
	// in order to satisfy the Teseo-II real-time requirements.
	sparse_matrix_t F_mat;
	reset_sparse_matrix(KNSTATES, &F_mat);

	add_element_to_sparse_matrix(0, 3, delta_t, &F_mat);
	add_element_to_sparse_matrix(1, 4, delta_t, &F_mat);
	add_element_to_sparse_matrix(2, 5, delta_t, &F_mat);

	// NOTE: Some of the following terms could be neglected in some configurations
	add_element_to_sparse_matrix(3, 0, +OMEGA_EARTH*OMEGA_EARTH*delta_t, &F_mat);
	add_element_to_sparse_matrix(4, 1, +OMEGA_EARTH*OMEGA_EARTH*delta_t, &F_mat);
	add_element_to_sparse_matrix(3, 4, +2*OMEGA_EARTH*delta_t, &F_mat);
	add_element_to_sparse_matrix(4, 3, -2*OMEGA_EARTH*delta_t, &F_mat);

	add_element_to_sparse_matrix( 9+KAC, 10+KAC, +OMEGA_EARTH*delta_t, &F_mat);
	add_element_to_sparse_matrix(10+KAC,  9+KAC, -OMEGA_EARTH*delta_t, &F_mat);

	add_element_to_sparse_matrix(12+KAC, 12+KAC, -invtau_acc_dt, &F_mat);
	add_element_to_sparse_matrix(13+KAC, 13+KAC, -invtau_acc_dt, &F_mat);
	add_element_to_sparse_matrix(14+KAC, 14+KAC, -invtau_acc_dt, &F_mat);
	add_element_to_sparse_matrix(15+KAC, 15+KAC, -invtau_gyr_dt, &F_mat);
	add_element_to_sparse_matrix(16+KAC, 16+KAC, -invtau_gyr_dt, &F_mat);
	add_element_to_sparse_matrix(17+KAC, 17+KAC, -invtau_gyr_dt, &F_mat);

	// Clock drift term.
	add_element_to_sparse_matrix( 6,  8+KAC, delta_t, &F_mat);

	add_element_to_sparse_matrix( 3, 10+KAC, -acc_e[2]*delta_t, &F_mat);
	add_element_to_sparse_matrix( 3, 11+KAC, +acc_e[1]*delta_t, &F_mat);
	add_element_to_sparse_matrix( 4,  9+KAC, +acc_e[2]*delta_t, &F_mat);
	add_element_to_sparse_matrix( 4, 11+KAC, -acc_e[0]*delta_t, &F_mat);
	add_element_to_sparse_matrix( 5,  9+KAC, -acc_e[1]*delta_t, &F_mat);
	add_element_to_sparse_matrix( 5, 10+KAC, +acc_e[0]*delta_t, &F_mat);

	add_element_to_sparse_matrix( 3, 12+KAC, +Cbe[0][0]*delta_t, &F_mat);
	add_element_to_sparse_matrix( 3, 13+KAC, +Cbe[0][1]*delta_t, &F_mat);
	add_element_to_sparse_matrix( 3, 14+KAC, +Cbe[0][2]*delta_t, &F_mat);
	add_element_to_sparse_matrix( 4, 12+KAC, +Cbe[1][0]*delta_t, &F_mat);
	add_element_to_sparse_matrix( 4, 13+KAC, +Cbe[1][1]*delta_t, &F_mat);
	add_element_to_sparse_matrix( 4, 14+KAC, +Cbe[1][2]*delta_t, &F_mat);
	add_element_to_sparse_matrix( 5, 12+KAC, +Cbe[2][0]*delta_t, &F_mat);
	add_element_to_sparse_matrix( 5, 13+KAC, +Cbe[2][1]*delta_t, &F_mat);
	add_element_to_sparse_matrix( 5, 14+KAC, +Cbe[2][2]*delta_t, &F_mat);

	add_element_to_sparse_matrix( 9+KAC, 15+KAC, -Cbe[0][0]*delta_t, &F_mat);
	add_element_to_sparse_matrix( 9+KAC, 16+KAC, -Cbe[0][1]*delta_t, &F_mat);
	add_element_to_sparse_matrix( 9+KAC, 17+KAC, -Cbe[0][2]*delta_t, &F_mat);
	add_element_to_sparse_matrix(10+KAC, 15+KAC, -Cbe[1][0]*delta_t, &F_mat);
	add_element_to_sparse_matrix(10+KAC, 16+KAC, -Cbe[1][1]*delta_t, &F_mat);
	add_element_to_sparse_matrix(10+KAC, 17+KAC, -Cbe[1][2]*delta_t, &F_mat);
	add_element_to_sparse_matrix(11+KAC, 15+KAC, -Cbe[2][0]*delta_t, &F_mat);
	add_element_to_sparse_matrix(11+KAC, 16+KAC, -Cbe[2][1]*delta_t, &F_mat);
	add_element_to_sparse_matrix(11+KAC, 17+KAC, -Cbe[2][2]*delta_t, &F_mat);

	// State noise covariance propagation
	// P(k) = PHI*P(k-1)*PHI' + Q(k);
	//
	// Get M*P(k-1)*M'
#if USE_FAST_MECHANIZATION == 1
	update_matrix_sparse_fast(&F_mat, gicsrx->exchange.kal_cov_mat);
#else
	update_matrix_sparse(KNSTATES, &F_mat, gicsrx->exchange.kal_cov_mat);
#endif
}

static void propagate_kalimu_stochastic(const kconf_t * kconf, const float acc_e[3], float * kal_cov_mat, float delta_t)
{
	// Calculate different factors that will be used in the function repeatedly.
	float invtau_acc_dt  = 1/kconf->tau_acc*delta_t;
	float invtau_gyr_dt  = 1/kconf->tau_gyr*delta_t;

	float delta_t2 = delta_t  * delta_t / 2.0;
	float delta_t3 = delta_t * delta_t * delta_t / 3.0;

	// Get P(k) = M*P(k-1)*M' + Q(k)
	// Include IMU process noise model.
	float sigma2_bacc_dt = kconf->sigma2_bacc*delta_t;
	update_sym_element(kal_cov_mat,  3,  3, sigma2_bacc_dt);
	update_sym_element(kal_cov_mat,  4,  4, sigma2_bacc_dt);
	update_sym_element(kal_cov_mat,  5,  5, sigma2_bacc_dt);

	float sigma2_bgyr_dt = kconf->sigma2_bgyr*delta_t;
	update_sym_element(kal_cov_mat,  9+KAC,  9+KAC, sigma2_bgyr_dt);
	update_sym_element(kal_cov_mat, 10+KAC, 10+KAC, sigma2_bgyr_dt);
	update_sym_element(kal_cov_mat, 11+KAC, 11+KAC, sigma2_bgyr_dt);

	float sigma2_wnacc = 2 * kconf->sigma2_wnacc * invtau_acc_dt;
	update_sym_element(kal_cov_mat, 12+KAC, 12+KAC, sigma2_wnacc);
	update_sym_element(kal_cov_mat, 13+KAC, 13+KAC, sigma2_wnacc);
	update_sym_element(kal_cov_mat, 14+KAC, 14+KAC, sigma2_wnacc);

	float sigma2_wngyr = 2 * kconf->sigma2_wngyr * invtau_gyr_dt;
	update_sym_element(kal_cov_mat, 15+KAC, 15+KAC, sigma2_wngyr);
	update_sym_element(kal_cov_mat, 16+KAC, 16+KAC, sigma2_wngyr);
	update_sym_element(kal_cov_mat, 17+KAC, 17+KAC, sigma2_wngyr);

	float sigma2_bacc_dt2 = kconf->sigma2_bacc * delta_t2;
	float sigma2_bacc_dt3 = kconf->sigma2_bacc * delta_t3;
	update_sym_element(kal_cov_mat,  0,  3, sigma2_bacc_dt2);
	update_sym_element(kal_cov_mat,  1,  4, sigma2_bacc_dt2);
	update_sym_element(kal_cov_mat,  2,  5, sigma2_bacc_dt2);
	update_sym_element(kal_cov_mat,  0,  0, sigma2_bacc_dt3);
	update_sym_element(kal_cov_mat,  1,  1, sigma2_bacc_dt3);
	update_sym_element(kal_cov_mat,  2,  2, sigma2_bacc_dt3);

	float sigma2_bgyr_dt2 = kconf->sigma2_bgyr * delta_t2;
	update_sym_element(kal_cov_mat,  3, 10+KAC, -acc_e[2]*sigma2_bgyr_dt2);
	update_sym_element(kal_cov_mat,  3, 11+KAC, +acc_e[1]*sigma2_bgyr_dt2);
	update_sym_element(kal_cov_mat,  4,  9+KAC, +acc_e[2]*sigma2_bgyr_dt2);
	update_sym_element(kal_cov_mat,  4, 11+KAC, -acc_e[0]*sigma2_bgyr_dt2);
	update_sym_element(kal_cov_mat,  5,  9+KAC, -acc_e[1]*sigma2_bgyr_dt2);
	update_sym_element(kal_cov_mat,  5, 10+KAC, +acc_e[0]*sigma2_bgyr_dt2);

	float sigma2_bgyr_dt3 = kconf->sigma2_bgyr * delta_t3;
	update_sym_element(kal_cov_mat,  3,  3, acc_e[2]*acc_e[2]*sigma2_bgyr_dt3);
	update_sym_element(kal_cov_mat,  3,  3, acc_e[1]*acc_e[1]*sigma2_bgyr_dt3);
	update_sym_element(kal_cov_mat,  4,  4, acc_e[2]*acc_e[2]*sigma2_bgyr_dt3);
	update_sym_element(kal_cov_mat,  4,  4, acc_e[0]*acc_e[0]*sigma2_bgyr_dt3);
	update_sym_element(kal_cov_mat,  5,  5, acc_e[1]*acc_e[1]*sigma2_bgyr_dt3);
	update_sym_element(kal_cov_mat,  5,  5, acc_e[0]*acc_e[0]*sigma2_bgyr_dt3);

	// Include clock drift process noise model.
	update_sym_element(kal_cov_mat,  6,  6, kconf->sigma2_clk*delta_t3);
	update_sym_element(kal_cov_mat,  7,  7, kconf->sigma2_isb_g1 * delta_t);
#if KAC == 0
	update_sym_element(kal_cov_mat,  8,  8, kconf->sigma2_clk    * delta_t);
	update_sym_element(kal_cov_mat,  8,  6, kconf->sigma2_clk    * delta_t2);
#else
	update_sym_element(kal_cov_mat,  8,  8, kconf->sigma2_isb_g2 * delta_t);
	update_sym_element(kal_cov_mat,  9,  9, kconf->sigma2_clk    * delta_t);
	update_sym_element(kal_cov_mat,  9,  6, kconf->sigma2_clk    * delta_t2);
#endif
}

/* This function propagates the covariance matrix of the Hybrid Kalman Filter, taking into
 * account the current value of the accelerometer (in body-frame coordinates, corrected with
 * sensor biases). The propagation is performed throughout a delta_t time period.
 */
static void propagate_kalimu(gicsrx_t *gicsrx, const float acc_b[3], float delta_t)
{
	float acc_e [3];
	propagate_kalimu_deterministic(gicsrx, acc_b, delta_t, acc_e);
	propagate_kalimu_stochastic(&gicsrx->kconf, acc_e, gicsrx->exchange.kal_cov_mat, delta_t);
}

/* This function shall be invoked every time a new pair of samples (gyro, accelerometer) is received, in
 * order to mechanize/propagate the state vector (position and velocity) and the KF covariance matrix.
 * It can be configured by the programmer to use a constant gravity vector, which is provided as input,
 * taking advantage of the fact that the variations throughout a second are tiny, or to calculate it
 * each time the function is invoked.
 */
#if CONSTANT_GRAVITY == 1
static void mechanization_step(gicsrx_t *gicsrx, const float rawacc[3], const float rawgyr[3], const float geod_acc[3], const float delta_t, float acc_n[3])
#else
static void mechanization_step(gicsrx_t *gicsrx, const float rawacc[3], const float rawgyr[3], const float delta_t, float acc_n[3])
#endif
{
#if CONSTANT_GRAVITY == 0
	float  geod_acc[3];
#endif

	float  prop_vel[3], Cbn[3][3], acc_b[3], gyr_ib[3];
	double M, N, *ins_state;

	// Avoid transitory periods in the mechanization process, checking that the timestamp is below 1 second limit.
	if (delta_t > 1)
	{
		return;
	}

	// Obtain pointer to the propagated state vector (in navigation coordinate-frame).
	ins_state = &gicsrx->exchange.ins_state_nav[0];

	// Calculate M and N factors.
	compute_M_N(ins_state[0], &M, &N);

#if CONSTANT_GRAVITY == 0
#if CORRECT_ACC == 0
	calculate_geod_acc(&ins_state[0], geod_acc);
#else
	calculate_geod_acc(&ins_state[0], &ins_state[3], M, N, geod_acc);
#endif
#endif

	// Correct accelerometer bias.
	acc_b[0] = rawacc[0] - ins_state[12+KAC];
	acc_b[1] = rawacc[1] - ins_state[13+KAC];
	acc_b[2] = rawacc[2] - ins_state[14+KAC];

	// Get transformation matrix (from body to navigation).
	quaternion2mat(gicsrx->exchange.quaternion, Cbn);

	// Convert accelerometer measurement to navigation coordinate frame, and add gravitational model.
	matVecMul_f(acc_n, Cbn, acc_b, 3, 3);
	matVecSum  (acc_n, acc_n, geod_acc, 3);

	// This is the velocity used to propagate the position.
#if MECHANIZED_VEL_TYPE == 0
	// Mechanization used in ATLANTIDA.
	prop_vel[0]  = +(ins_state[3])/((M + ins_state[2]));
	prop_vel[1]  = +(ins_state[4])/((N + ins_state[2])*cos(ins_state[0]));
	prop_vel[2]  = -(ins_state[5]);
#elif MECHANIZED_VEL_TYPE == 1
	// Mechanization used in GARLIC.
	prop_vel[0]  = +(ins_state[3] + 0.5*delta_t*acc_n[0])/((M + ins_state[2]));
	prop_vel[1]  = +(ins_state[4] + 0.5*delta_t*acc_n[1])/((N + ins_state[2])*cos(ins_state[0]));
	prop_vel[2]  = -(ins_state[5] + 0.5*delta_t*acc_n[2]);
#else
	prop_vel[0]  = +(ins_state[3] + delta_t*acc_n[0])/((M + ins_state[2]));
	prop_vel[1]  = +(ins_state[4] + delta_t*acc_n[1])/((N + ins_state[2])*cos(ins_state[0]));
	prop_vel[2]  = -(ins_state[5] + delta_t*acc_n[2]);
#endif

	// Propagate vehicle position.
	ins_state[0] += delta_t*prop_vel[0];
	ins_state[1] += delta_t*prop_vel[1];
	ins_state[2] += delta_t*prop_vel[2];

	// Propagate vehicle velocity.
	ins_state[3] += delta_t*acc_n[0];
	ins_state[4] += delta_t*acc_n[1];
	ins_state[5] += delta_t*acc_n[2];

	// Correct gyro bias.
	gyr_ib[0] = rawgyr[0] - ins_state[15+KAC];
	gyr_ib[1] = rawgyr[1] - ins_state[16+KAC];
	gyr_ib[2] = rawgyr[2] - ins_state[17+KAC];

	// Correct attitude through quaternion update.
#if CORRECT_ANGLE == 1
	update_quaternion(gicsrx->exchange.quaternion, &ins_state[0], &ins_state[3], gyr_ib, delta_t);
#else
	update_quaternion(gicsrx->exchange.quaternion, gyr_ib, delta_t);
#endif

	/*
	// Update attitude matrix.
	quaternion2mat(gicsrx->exchange.quaternion, Cbn);

	// Calculate attitude angles from quaternions.
	get_attitude_angles_from_matrix(Cbn, &ins_state[9], &ins_state[10], &ins_state[11]);
	 */

#if CLK_STEERING == 0
	// Propagate clock model. Notice that drift step is not taken into account at this point.
	ins_state[6] += delta_t*ins_state[8+KAC]; // To be used only with clock-steering OFF in SRX-10.
#endif

	propagate_kalimu(gicsrx, acc_b, delta_t);
}

static void calculate_distance_metrics(gicsrx_t *gicsrx, nav_sol_t *nav_sol, double raw_imu[3])
{
	static char calc_init_dist = 1;
	if(calc_init_dist == 1)
	{
		double tmp_tow = gicsrx->exchange.tow;
		if(nav_sol->quality >= KALMAN_ONLY_1)
		{
			gicsrx->exchange.initial_distance = calc_initial_distance(&tmp_tow, raw_imu, &gicsrx->exchange, 1);
			gicsrx->exchange.initial_tow = tmp_tow;
			calc_init_dist = 0;
		}
		else
		{
			calc_initial_distance(&tmp_tow, raw_imu, &gicsrx->exchange, 0);
		}
	}
	nav_sol->t_dist_0   = gicsrx->exchange.initial_tow;
	nav_sol->distance_0 = gicsrx->exchange.initial_distance;

	// Used for distance computation.
	char prev_state = (nav_sol->prev_used_sats > 0 && nav_sol->kalman_valid == 1),
			curr_state = (nav_sol->curr_used_sats > 0 && nav_sol->kalman_valid == 1);

	// Update the travelled distance using the Kalman filter velocity.
	update_distance(nav_sol, nav_sol->prop_time, prev_state, curr_state, &gicsrx->exchange.travelled_distance);

	nav_sol->distance_t = gicsrx->exchange.travelled_distance;
}

static void calculate_spoofing_metrics(gicsrx_t *gicsrx, nav_sol_t *nav_sol, double acc_n_1Hz[3])
{
	double acc_bias[3] = {gicsrx->exchange.ins_state_ecef[12+KAC],
                          gicsrx->exchange.ins_state_ecef[13+KAC],
                          gicsrx->exchange.ins_state_ecef[14+KAC]};

	// If biases are too high, the GNSS+IMU KF status is not OK. Force reset.
	if(Vec3Norm(acc_bias) > 1.0) { gicsrx->exchange.kalman_on = 0; }

#if USE_ANTISPOOFING == 1

	double as_velocity[3], acc_n_kal[3];
	static double prop_kalvel[3] = {0, 0, 0};

#if CALC_GNSS_ONLY == 1
	ECEFtoNAV_vel(&gicsrx->exchange.gnss_state[3], as_velocity, gicsrx->exchange.ins_state_nav[0], gicsrx->exchange.ins_state_nav[1]);
#else
	as_velocity[0] = gicsrx->exchange.ins_state_nav[3];
	as_velocity[1] = gicsrx->exchange.ins_state_nav[4];
	as_velocity[2] = gicsrx->exchange.ins_state_nav[5];
#endif

	acc_n_kal[0] = as_velocity[0] - prop_kalvel[0];
	acc_n_kal[1] = as_velocity[1] - prop_kalvel[1];
	acc_n_kal[2] = as_velocity[2] - prop_kalvel[2];

	prop_kalvel[0] = as_velocity[0];
	prop_kalvel[1] = as_velocity[1];
	prop_kalvel[2] = as_velocity[2];

	antispoofing_algorithm(gicsrx, acc_n_1Hz, acc_n_kal, nav_sol);
#endif

#if CALC_GNSS_ONLY == 1
	nav_sol->discr_alarm = check_filter_inconsistencies(gicsrx->exchange.ins_state_ecef, gicsrx->exchange.gnss_state, gicsrx->exchange.gnss_used_sats);
	if(nav_sol->discr_alarm == 1)
	{
		gicsrx->exchange.kalman_on = 0;
	}
#endif
}

#if CONSTANT_GRAVITY == 1
static void mechanize_current_sample(gicsrx_t *gicsrx, gnssdata_t *gnssdata, sensor_t *imu_sample, float geod_acc[3],
		double acc_n_1Hz[3], double raw_imu[3], float *delta_t_accum, int *mech_steps, int left_samples)
#else
static void mechanize_current_sample(gicsrx_t *gicsrx, gnssdata_t *gnssdata, sensor_t *imu_sample,
		double acc_n_1Hz[3], double raw_imu[3], float *delta_t_accum, int *mech_steps, int left_samples)
#endif
{
	static int   ss_index    =  0;
	static float acc_b_ss[3] = {0};
	static float gyr_b_ss[3] = {0};

	matVecSum(acc_b_ss, acc_b_ss, imu_sample->acc_b, 3);
	matVecSum(gyr_b_ss, gyr_b_ss, imu_sample->gyr_b, 3);

	ss_index++;

	if(ss_index == gicsrx->kconf.IMU_rate || (left_samples + ss_index) <= gicsrx->kconf.IMU_rate)
	{
		double delta_t = imu_sample->gps_time - gicsrx->exchange.tow;
		*delta_t_accum += delta_t;

		int k;
		for(k = 0; k < 3; k++)
		{
			acc_b_ss[k] /= ss_index;
			gyr_b_ss[k] /= ss_index;
			raw_imu [k] += acc_b_ss[k];
		}

		float acc_n[3] = {0, 0, 0};
#if CONSTANT_GRAVITY == 1
		mechanization_step(gicsrx, acc_b_ss, gyr_b_ss, geod_acc, delta_t, acc_n);
#else
		mechanization_step(gicsrx, acc_b_ss, gyr_b_ss, delta_t, acc_n);
#endif
		acc_n_1Hz[0] += acc_n[0];
		acc_n_1Hz[1] += acc_n[1];
		acc_n_1Hz[2] += acc_n[2];

		gicsrx->exchange.tow = imu_sample->gps_time;

		// Set to 0 the cumulative variables.
		acc_b_ss[0] = acc_b_ss[1] = acc_b_ss[2] = 0;
		gyr_b_ss[0] = gyr_b_ss[1] = gyr_b_ss[2] = 0;
		ss_index    = 0;

		// Increase the number of samples used for mechanization in current epoch.
		(*mech_steps)++;
	}
}

#if CONSTANT_GRAVITY == 1
static void mechanize_without_samples(gicsrx_t *gicsrx, gnssdata_t *gnssdata, float geod_acc[3], float *delta_t_accum)
#else
static void mechanize_without_samples(gicsrx_t *gicsrx, gnssdata_t *gnssdata, float *delta_t_accum)
#endif
{
	const double delta_t_step = 0.01;

	if(gicsrx->exchange.kalman_on == 1)
	{
		*delta_t_accum = (gnssdata->tow - gicsrx->exchange.tow);

		int nsamples = (int)(*delta_t_accum /delta_t_step);

		double *ins_state = gicsrx->exchange.ins_state_nav;

		float acc_n[3] = {0, 0, 0};
		float acc_b[3] = {ins_state[12+KAC], ins_state[13+KAC], ins_state[14+KAC] - GRAVITY_G};
		float gyr_b[3] = {ins_state[15+KAC], ins_state[16+KAC], ins_state[17+KAC]};

		int k;
		for(k = 0; k < nsamples; k++)
		{
#if CONSTANT_GRAVITY == 1
			mechanization_step(gicsrx, acc_b, gyr_b, geod_acc, delta_t_step, acc_n);
#else
			mechanization_step(gicsrx, acc_b, gyr_b, delta_t_step, acc_n);
#endif
		}
	}
	gicsrx->exchange.tow = gnssdata->tow;
}

/* This function implements a Hybrid Kalman Filter and shall be invoked every epoch to perform properly the
 * propagation and update functions. In the Hybrid KF, the propagation functionality implements a mechanization
 * procedure, where data from inertial sensors (accelerometers and gyros) is incorporated into the navigation
 * system. The samples from IMU shall be read by an external function, since this function takes them as input.
 * Additional features are included, apart from navigation, as for instance consistency check of the filter,
 * distance computation or anti-spoofing algorithm (if activated).
 */
void kalman_hybrid(gicsrx_t *gicsrx, gnssdata_t *gnssdata, lsq_state_t *lsq_state, sensor_t imu[MAX_IMU_SAMPLES], int nsamples, nav_sol_t *nav_sol)
{
	int k;
	int mech_steps = 0;

	float delta_t_accum = 0.00;

	double acc_n_1Hz[3] = {0, 0, 0};
	double raw_imu  [3] = {0, 0, 0};

#if CONSTANT_GRAVITY == 1
	float geod_acc[3];
#if CORRECT_ACC == 0
	calculate_geod_acc(&gicsrx->exchange.ins_state_nav[0], geod_acc);
#else
	double M, N;
	compute_M_N(gicsrx->exchange.ins_state_nav[0], &M, &N);

	calculate_geod_acc(&gicsrx->exchange.ins_state_nav[0], &gicsrx->exchange.ins_state_nav[3], M, N, geod_acc);
#endif
#endif
	for(k = 0; k < nsamples; k++)
	{
#if CONSTANT_GRAVITY == 1
		mechanize_current_sample(gicsrx, gnssdata, &imu[k], geod_acc, acc_n_1Hz, raw_imu, &delta_t_accum, &mech_steps, nsamples - k);
#else
		mechanize_current_sample(gicsrx, gnssdata, &imu[k], acc_n_1Hz, raw_imu, &delta_t_accum, &mech_steps, nsamples - k);
#endif
	}

	// If the buffer of IMU samples is not empty, the average acceleration at 1Hz (in navigation frame coordinates)
	// and the raw value of the accelerometer (in body frame) are calculated by dividing by the number of mechanization
	// steps.
	// If the buffer of IMU samples is empty, it means that the IMU samples are not ready for the current epoch, for two
	// possible reasons: first, there is simply no data; second, a Teseo-II clock reset has been produced and a wrong
	// timestamp has been detected, such that the IMU data related to the current epoch is flagged as unreliable.
	if(nsamples > 0)
	{
		for(k = 0; k < 3; k++)
		{
			acc_n_1Hz[k] /= mech_steps;
			raw_imu  [k] /= mech_steps;
		}
	}
	else
	{
		// In case no samples are available, mechanize the whole period between epochs (typically, 2 seconds) with perfect
		// samples of accelerometer and gyro equivalent to be static (considering also the corresponding sensor biases).
#if CONSTANT_GRAVITY == 1
		mechanize_without_samples(gicsrx, gnssdata, geod_acc, &delta_t_accum);
#else
		mechanize_without_samples(gicsrx, gnssdata, &delta_t_accum);
#endif
	}
	gicsrx->exchange.prop_time = delta_t_accum;

#if CALCULATE_IBPL == 1
	nav_sol->kalman_ibpl_valid = 0;
	gicsrx->exchange.kibpl.common.acc_norm_imu = Vec3Norm(acc_n_1Hz);
#endif
	nav_sol->raw_imu[0] = raw_imu[0];
	nav_sol->raw_imu[1] = raw_imu[1];
	nav_sol->raw_imu[2] = raw_imu[2];

#if CALC_GNSS_ONLY == 1
	kalman_gnss_only(gicsrx, gnssdata, lsq_state->pos_valid, lsq_state, nav_sol);
#endif
	kalman_gnss(gicsrx, gnssdata, lsq_state->pos_valid, lsq_state, nav_sol);

	calculate_distance_metrics (gicsrx, nav_sol, raw_imu);

#if CALC_GNSS_ONLY == 1
	if(gicsrx->exchange.kalman_on == 1 && gicsrx->exchange.kalman2_on == 1)
#else
		if(gicsrx->exchange.kalman_on == 1)
#endif
		{
			calculate_spoofing_metrics (gicsrx, nav_sol, acc_n_1Hz);
		}
}
#endif
