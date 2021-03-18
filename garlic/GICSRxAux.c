#include "matrix.h"
#include "algebra.h"
#include "GICSRxAux.h"

#if USE_IMU == 1

#define DISTANCE_WINDOW  60
#define DISCREPANCY_SIZE 10
#define DISCREPANCY_THR  50

#if USE_ANTISPOOFING == 1

#define AS_CONV_TIME  30
#define AS_NSATS_THR  8
#define AS_ACCDIV_THR 10

/* This function implements an anti-spoofing algorithm, which checks the consistency between the IMU
 * acceleration (properly corrected with the corresponding accelerometer biases, attitude and gravity, and
 * expressed in navigation-frame), and GNSS-based acceleration, calculated as the difference between
 * consecutive velocity estimations in a time lapse of 1-second, since the current implementation of the
 * Kalman Filter does not include acceleration in the state vector.
 * The IMU acceleration can be obtained, for instance, as the averaged value over a second of 100 samples,
 * assuming that its value remain constant throughout this period.
 */
double antispoofing_algorithm(gicsrx_t *gicsrx, double acc_imu[3], double acc_kal[3], nav_sol_t *nav_sol)
{
	// Dynamic variance of the error.
	static double as_variance = 0;

	// This window stores dynamically the accelerations calculated using the KF [][0][], and the IMU [][1][].
	// The purpose is to smooth the noise of the acceleration averaging the AS_WSIZE samples in the array.
	static double antispoofing_window[AS_WSIZE][2][3] = {{{0}}};

	// Shift window and store the new observations.
	int k;
	for(k = 0; k < AS_WSIZE-1; k++)
	{
		antispoofing_window[k][0][0] = antispoofing_window[k+1][0][0];
		antispoofing_window[k][0][1] = antispoofing_window[k+1][0][1];
		antispoofing_window[k][0][2] = antispoofing_window[k+1][0][2];

		antispoofing_window[k][1][0] = antispoofing_window[k+1][1][0];
		antispoofing_window[k][1][1] = antispoofing_window[k+1][1][1];
		antispoofing_window[k][1][2] = antispoofing_window[k+1][1][2];
	}
	antispoofing_window[AS_WSIZE-1][0][0] = acc_kal[0];
	antispoofing_window[AS_WSIZE-1][0][1] = acc_kal[1];
	antispoofing_window[AS_WSIZE-1][0][2] = acc_kal[2];

	antispoofing_window[AS_WSIZE-1][1][0] = acc_imu[0];
	antispoofing_window[AS_WSIZE-1][1][1] = acc_imu[1];
	antispoofing_window[AS_WSIZE-1][1][2] = acc_imu[2];

	// Calculate an averaged estimation of the acceleration over a time period with AS_WSIZE seconds.
	// This average is performed for each acceleration component.
	double acc_kal_f[3] = {0, 0, 0}, acc_imu_f[3] = {0, 0, 0};
	for(k = 0; k < AS_WSIZE; k++)
	{
		matVecSum(acc_kal_f, acc_kal_f, antispoofing_window[k][0], 3);
		matVecSum(acc_imu_f, acc_imu_f, antispoofing_window[k][1], 3);
	}

	// Calculate the modulus of the smoothed accelerations.
	double norm_kal_f   = Vec3Norm(acc_kal_f)/AS_WSIZE;
	double norm_imu_f   = Vec3Norm(acc_imu_f)/AS_WSIZE;

	// The error signal is defined as the difference between the acceleration estimated by the KF and
	// the obtained from IMU. Shall be remarked that the IMU case is also altered by the KF, since
	// attitude and bias is corrected prior to enter the anti-spoofing algorithm.
	double norm_error_f = (norm_kal_f - norm_imu_f);

	// Check that the KF is turned on and the number of satellites used. It is desirable to work with
	// the GNSS-only Kalman Filter, such that we have an independent source for acceleration estimation,
	// instead of using the Hybrid information, which is coupled with the IMU, due to the bias and attitude
	// corrections.
#if CALC_GNSS_ONLY == 1
	int nsats_used = gicsrx->exchange.gnss_used_sats;
#else
	int nsats_used = gicsrx->exchange.kal_used_sats;
#endif

	// If the following conditions are met (convergence, visibility and maximum divergence, to avoid false alarms)
	// the anti-spoofing variance estimation is updated.
	if((gicsrx->exchange.tow - gicsrx->exchange.timeOfLastReset > AS_CONV_TIME) &&
	   (nsats_used >= AS_NSATS_THR) && (fabs(norm_error_f) < AS_ACCDIV_THR))
	{
		as_variance = AS_LAMBDA*as_variance + norm_error_f*norm_error_f;
	}

	// The variance value must be normalized using the effective window size.
	double as_normvar    = as_variance/AS_EFFWIN;
	nav_sol->spoof_var   = as_normvar;

	// Check that current variance is not above the alarm limit.
	nav_sol->spoof_alarm = (as_normvar > AS_THRESHOLD);

	return as_normvar;
}
#endif

#if CALC_GNSS_ONLY == 1
/* This function calculates, when both Hybrid and GNSS-only Kalman Filters are active, the discrepancy
 * between both solutions, in terms of position. When the distance between both is higher than a
 * pre-determined threshold, an alarm flag is returned, so the operator can take specific decisions.
 */
char check_filter_inconsistencies(double hybrid_state[KNSTATES], double gnss_state[GNSS_KSTATES], char noOfSats_GNSS)
{
	char  alarm_signal    = 0;
	int   numOfPrevAlarms = 0;

	// The purpose of this function is to detect inconsistencies between the GNSS+IMU and the GNSS-only KF.
	// According to REQ-GARLIC6-PRF-010 in SRD, the system must detect divergences between both filters higher
	// than 50 meters which hold for more than 10 seconds.

	// This vector contains flags which determine the consistency for an epoch inside the observation window.
	static char consistency_window[DISCREPANCY_SIZE] = {0};

	// Calculate distance between navigation solutions.
	double diff_vector[3] = {hybrid_state[0] - gnss_state[0], hybrid_state[1] - gnss_state[1], hybrid_state[2] - gnss_state[2]};
	double distance       = Vec3Norm(diff_vector);

	// Shift the flags vector and calculate the number of alarms present in the window.
	int k;
	for(k = 0; k < DISCREPANCY_SIZE-1; k++)
	{
		numOfPrevAlarms      += consistency_window[k+1];
		consistency_window[k] = consistency_window[k+1];
	}
	numOfPrevAlarms                       += ((distance > DISCREPANCY_THR) && (noOfSats_GNSS >= 8));
	consistency_window[DISCREPANCY_SIZE-1] = ((distance > DISCREPANCY_THR) && (noOfSats_GNSS >= 8));

	// If all the epochs are set as divergent, set alarm signal to 1.
	if(numOfPrevAlarms == DISCREPANCY_SIZE)
	{
		alarm_signal = 1;
	}
	return alarm_signal;
}
#endif

/* This function calculates the initial distance, using for this purpose a buffer that stores the
 * raw imu accelerations (averaged 100 samples in one second) during a certain period (defined by
 * the length of the distance window). In REQ-GARLIC6-PRF-040, a 60 second period previous to the
 * first fix shall be integrated to obtain the initial distance.
 */
float calc_initial_distance(double *tow_imu, double raw_imu[3], exchange_t *exchange, char get_curr_dist)
{
	// The current implementation of the algorithm is simplified, since both the bias and attitude
	// estimation need to be improved prior to accurately correct gravity vector and biases.
	// Therefore this version considers that the gravity vector is pointed through the z-axis of
	// the sensor coordinate frame.
	const float gravity_geod[3] = {0.00, 0.00, GRAVITY_G};

	static int no_of_stored_seconds = 0;

	// This window stores the raw imu samples (corrected with gravity vector) during the required
	// time period.
	static float tow_window[DISTANCE_WINDOW]    =  {0};
	static float acc_window[DISTANCE_WINDOW][3] = {{0}};

	float distance = 0.0;

	int k;
	float mean_acc[3] = {0,0,0};

	if(no_of_stored_seconds > 0 && fabs(*tow_imu - tow_window[no_of_stored_seconds-1]) > 1.5)
	{
		no_of_stored_seconds = 0;
	}

	if(no_of_stored_seconds < DISTANCE_WINDOW)
	{
		tow_window[no_of_stored_seconds]    = *tow_imu;
		acc_window[no_of_stored_seconds][0] = raw_imu[0] + gravity_geod[0];
		acc_window[no_of_stored_seconds][1] = raw_imu[1] + gravity_geod[1];
		acc_window[no_of_stored_seconds][2] = raw_imu[2] + gravity_geod[2];

		no_of_stored_seconds++;
	}
	else
	{
		// Shift the raw imu accelerations in the observation window, correcting gravity in the current sample.
		// The current version of the algorithm is only capable of estimating more or less accurately the bias
		// component in the direction of the gravity vector. Therefore, the approach used is to calculate the
		// average value of each component (or low-pass value), interpreting it as the component bias. This
		// approach assumes that each acceleration component should be centered around 0, which is not normally
		// true.
		for(k = 0; k < DISTANCE_WINDOW-1; k++)
		{
			tow_window[k] = tow_window[k+1];

			acc_window[k][0] = acc_window[k+1][0];
			acc_window[k][1] = acc_window[k+1][1];
			acc_window[k][2] = acc_window[k+1][2];
		}
		acc_window[DISTANCE_WINDOW-1][0] = raw_imu[0] + gravity_geod[0];
		acc_window[DISTANCE_WINDOW-1][1] = raw_imu[1] + gravity_geod[1];
		acc_window[DISTANCE_WINDOW-1][2] = raw_imu[2] + gravity_geod[2];
	}
	for(k = 0; k < no_of_stored_seconds; k++)
	{
		// Calculate average acceleration.
		matVecSum(mean_acc, mean_acc, acc_window[k], 3);
	}

	// Calculate the average acceleration over the period dividing by the number of samples.
	mean_acc[0] /= no_of_stored_seconds;
	mean_acc[1] /= no_of_stored_seconds;
	mean_acc[2] /= no_of_stored_seconds;

	// If the user requires to get the current distance value, compute it. As the vehicle can be
	// in movement when starting the distance computation, a initial velocity of 0 m/s might not
	// be realistic. Instead of this, the current velocity is set as initial velocity and the
	// integration is performed backwards.
	if(get_curr_dist == 1)
	{
		float Cbn[3][3], Cnb[3][3];

		// The initial velocity is the current velocity vector.
		float vnav[3] = {exchange->ins_state_nav[3], exchange->ins_state_nav[4], exchange->ins_state_nav[5]};
		float vbody[3];

		quaternion2mat(exchange->quaternion, Cbn);
		Transpose(Cnb, Cbn, 3,3);
		matVecMul_f(vbody, Cnb, vnav, 3,3);

		// Integrate distance backwards.
		for(k = no_of_stored_seconds-1; k >= 0; k--)
		{
			distance += Vec3Norm(vbody);

			float acc_body[3] = {(acc_window[k][0] - mean_acc[0]),
					(acc_window[k][1] - mean_acc[1]),
					(acc_window[k][2] - mean_acc[2])};

			vbody[0] -= acc_body[0];
			vbody[1] -= acc_body[1];
			vbody[2] -= acc_body[2];
		}
		*tow_imu = tow_window[0];
	}
	return distance;
}
#endif

/* This function updates the value of the distance travelled by the vehicle, considering for that purpose the
 * Kalman Filter velocity. To avoid divergences, two status flag, referred to the quality of the current and the
 * previous epoch are required. Distance is only computed when both flags are set to 1. If they are not, last
 * valid position is stored to update the distance count when possible.
 */
void update_distance(const nav_sol_t *nav_sol, float delta_t, char prev_state, char curr_state, double *distance)
{
	// Static variables.
	static char   first = 1;
	static double last_valid_ecef[3] = {0};

	// Check whether filter state in the current epoch is OK.
	if(curr_state == 1)
	{
		// Check whether filter state in the previous epoch is OK.
		if(prev_state == 1)
		{
			// Calculate the module of the velocity.
			float kalvel = Vec3Norm(nav_sol->velocity);

			// In order to take into account the current velocity measurement, check whether it is above
			// a noise threshold.
			if(kalvel >= VEL_THRESHOLD)
			{
				*distance += kalvel*delta_t;
			}
		}
		else
		{
			// Check whether the last stored position is valid.
			if(first == 0)
			{
				// Distance variables.
				double dkalpos[3], normdpos;

				// Calculate distance vector between current state position and stored one.
				dkalpos[0] = nav_sol->position[0] - last_valid_ecef[0];
				dkalpos[1] = nav_sol->position[1] - last_valid_ecef[1];
				dkalpos[2] = nav_sol->position[2] - last_valid_ecef[2];

				// Calculate the module of the distance vector and update the distance element.
				normdpos = Vec3Norm(dkalpos);
				*distance += normdpos;
			}
			else
			{
				first = 0;
			}
		}
		if(nav_sol->quality >= KALMAN_ONLY_1)
		{
			// Store the last valid position.
			last_valid_ecef[0] = nav_sol->position[0];
			last_valid_ecef[1] = nav_sol->position[1];
			last_valid_ecef[2] = nav_sol->position[2];
		}
	}
}

