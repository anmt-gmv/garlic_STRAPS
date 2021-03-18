#ifndef _ALGEBRA_H
#define _ALGEBRA_H

#include "GICSRx_defines.h"

void compute_median(float *list, int size, float *median);

#if USE_IMU == 1
/* This function extracts the attitude angles (heading, pitch, roll) from the attitude matrix Cbn */
void get_attitude_angles_from_matrix(float Cbn[3][3], double * angle_1, double * angle_2, double *angle_3);
#endif

/* This function calculates the rotation matrix from ECEF to NED coordinate frame Cen */
void ECEFtoNAV_mat(float Cen[3][3], double latitude, double longitude);

/* This function calculates the rotation matrix from NED to ECEF coordinate frame Cne */
void NAVtoECEF_mat(float Cne[3][3], double latitude, double longitude);

/* Convert a WGS-84 ECEF position into geodetic coordinates (Lat, Lon, h).
 * Latitude and longitude are expressed in radians.
 */
void ECEFtoNAV_pos(const double pos_ECEF[3], double pos_NAV[3]);

/* Convert a geodetic position (Lat, Lon, h) into WGS-84 ECEF coordinates.
 * Latitude and longitude are expressed in radians.
 */
void NAVtoECEF_pos(double pos_ECEF[3], double pos_NAV[3]);

/* This function convert a NED velocity into WGS-84 ECEF velocity
 * Latitude and longitude are expressed in radians.
 */
void NAVtoECEF_vel(double vel_ECEF[3], double vel_LH[3], double latitude, double longitude);

/* This function convert a WGS-84 ECEF velocity into NED velocity
 * Latitude and longitude are expressed in radians.
 */
void ECEFtoNAV_vel(double vel_ECEF[3], double vel_LH[3], double latitude, double longitude);

#if USE_IMU == 1

/* This function calculates the M and N factors in terms of latitude */
void compute_M_N(double latitude, double *M, double *N);

/* This function calculates the gravity vector expressed in navigation coordinate frame with respect to
 * a user position, expressed also in navigation coordinate frame. It is allowed to activate the correction
 * of additional terms, as Coriolis or centripetal forces, for more accuracy.
 */
#if CORRECT_ACC == 0
void calculate_geod_acc(double geod_pos[3], float geod_acc[3]);
#else
void calculate_geod_acc(double geod_pos[3], double geod_vel[3], double M, double N, float geod_acc[3]);
#endif

/* This function convert a quaternion into its associated attitude matrix. */
void quaternion2mat(float q[4], float Cba[3][3]);

/* This function converts the attitude matrix into its associated quaternion vector.
 * Since float data is used, it is suitable that numerical errors can bring a divergence in
 * the procedure. Thus, a classification of the information is performed before calculate it,
 * in order to apply the most accurate and numerically stable procedure.
 */
void mat2quaternion(float Cba[3][3], float q[4]);

#if CORRECT_ANGLE == 1
/* This function calculates the angular velocity of the rotation between inertial and navigation coordinate
 * frame, expressed in navigation coordinate frame.
 */
void calculate_win_n(float win_n[3], double geodpos[3], double geodvel[3]);
#endif

/* This function updates the quaternion vector applying the angular velocity from gyr_corr (gyro measurements
 * in body frame corrected with bias) throughout a time period delta_t. In case that the programmer would want
 * to correct the earth rotation, increasing accuracy, position and velocity in navigation coordinate frame
 * shall also be provided.
 */
#if CORRECT_ANGLE == 1
void update_quaternion(float q[4], double geodpos[3], double geodvel[3], float gyr_corr[3], float delta_t);
#else
void update_quaternion(float q[4], float gyr_corr[3], float delta_t);
#endif

#endif

#endif //_ALGEBRA_H
