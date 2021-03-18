#include "matrix.h"
#include "algebra.h"

#include <stdlib.h>

static int compare_float_function(const void *a, const void *b)
{
	float *x = (float *) a;
	float *y = (float *) b;

	if (*x < *y)
	{
		return -1;
	}
	else if (*x > *y)
	{
		return 1;
	}

	return 0;
}

void compute_median(float *list, int size, float *median)
{
	if (size == 0) {
		*median = 0.0;
		return;
	}
	qsort (list, size, sizeof(float), compare_float_function);

	int k = size / 2;
	if (size % 2 == 0) {
		*median = (list[k-1] + list[k]) / 2.0;
	} else {
		*median = list[k];
	}
}

#if USE_IMU == 1

/* This function allows to generate a continuous evolution of the attitude angles. */
static double get_continuous_angle(double new_angle, double old_angle)
{
	int num_cycles = (int)((new_angle - old_angle) / (2*PI_G) + 0.5);
	return (new_angle - num_cycles*2*PI_G);
}

/* This function extracts the attitude angles (heading, pitch, roll) from the attitude matrix Cbn */
void get_attitude_angles_from_matrix(float Cbn[3][3], double *angle_1, double *angle_2, double *angle_3)
{
	*angle_1 = get_continuous_angle(atan2(Cbn[1][0], Cbn[0][0]), *angle_1);
	*angle_2 = get_continuous_angle(-asin(Cbn[2][0]), *angle_2);
	*angle_3 = get_continuous_angle(atan2(Cbn[2][1], Cbn[2][2]), *angle_3);
}
#endif

/* This function calculates the rotation matrix from ECEF to NED coordinate frame Cen */
void ECEFtoNAV_mat(float Cen[3][3], double latitude, double longitude)
{
	float coslon, sinlon, coslat, sinlat;

	coslon = cos(longitude);
	sinlon = sin(longitude);
	coslat = cos(latitude);
	sinlat = sin(latitude);

	Cen[0][0] = -sinlat*coslon;
	Cen[0][1] = -sinlat*sinlon;
	Cen[0][2] = +coslat;
	Cen[1][0] = -sinlon;
	Cen[1][1] = +coslon;
	Cen[1][2] = +0.0;
	Cen[2][0] = -coslat*coslon;
	Cen[2][1] = -coslat*sinlon;
	Cen[2][2] = -sinlat;
}

/* This function calculates the rotation matrix from NED to ECEF coordinate frame Cne */
void NAVtoECEF_mat(float Cne[3][3], double latitude, double longitude)
{
	float coslon, sinlon, coslat, sinlat;

	coslon = cos(longitude);
	sinlon = sin(longitude);
	coslat = cos(latitude);
	sinlat = sin(latitude);

	Cne[0][0] = -sinlat*coslon;
	Cne[0][1] = -sinlon;
	Cne[0][2] = -coslat*coslon;
	Cne[1][0] = -sinlat*sinlon;
	Cne[1][1] = coslon;
	Cne[1][2] = -coslat*sinlon;
	Cne[2][0] = coslat;
	Cne[2][1] = 0.0;
	Cne[2][2] = -sinlat;
}

/* Convert a WGS-84 ECEF position into geodetic coordinates (Lat, Lon, h).
 * Latitude and longitude are expressed in radians.
 */
void ECEFtoNAV_pos(const double pos_ECEF[3], double pos_NAV[3])
{
	const float EPSCOST = 1e-5;

	double p, theta, N, a, b, sint, cost;

	a = R_EARTH;
	b = R_EARTH*sqrt(1. - E2_G);
	p = sqrt(pos_ECEF[0]*pos_ECEF[0] + pos_ECEF[1]*pos_ECEF[1]);

	theta = atan2(pos_ECEF[2], (p*sqrt(1.-E2_G)));
	sint  = sin(theta);
	cost  = cos(theta);

	pos_NAV[0] = atan2((pos_ECEF[2]+E12_G*b*sint*sint*sint), (p-E2_G*a*cost*cost*cost));
	pos_NAV[1] = atan2(pos_ECEF[1], pos_ECEF[0]);

	sint = sin(pos_NAV[0]);
	cost = cos(pos_NAV[0]);
	N    = a*a/ sqrt(a*a*cost*cost + b*b*sint*sint);

	if(fabs(cost) < EPSCOST)
	{
		pos_NAV[2] = pos_ECEF[2]/sqrt(1.-E2_G) - N;
	}
	else
	{
		pos_NAV[2] = p/cost - N;
	}
}

/* Convert a geodetic position (Lat, Lon, h) into WGS-84 ECEF coordinates.
 * Latitude and longitude are expressed in radians.
 */
void NAVtoECEF_pos(double pos_ECEF[3], double pos_NAV[3])
{
	double N, M;

	double sin_lat = sin(pos_NAV[0]);
	double cos_lat = cos(pos_NAV[0]);

	N = R_EARTH/sqrt(1. - E2_G*sin_lat*sin_lat);
	M = N + pos_NAV[2];

	pos_ECEF[0] = M*cos_lat*cos(pos_NAV[1]);
	pos_ECEF[1] = M*cos_lat*sin(pos_NAV[1]);
	pos_ECEF[2] = (N*(1.-E2_G) + pos_NAV[2])*sin_lat;
}

/* This function convert a NED velocity into WGS-84 ECEF velocity
 * Latitude and longitude are expressed in radians.
 */
void NAVtoECEF_vel(double vel_ECEF[3], double vel_LH[3], double latitude, double longitude)
{
	float Cne[3][3];

	NAVtoECEF_mat(Cne, latitude, longitude);

	vel_ECEF[0] = Cne[0][0]*vel_LH[0] + Cne[0][1]*vel_LH[1] + Cne[0][2]*vel_LH[2];
	vel_ECEF[1] = Cne[1][0]*vel_LH[0] + Cne[1][1]*vel_LH[1] + Cne[1][2]*vel_LH[2];
	vel_ECEF[2] = Cne[2][0]*vel_LH[0] + Cne[2][1]*vel_LH[1] + Cne[2][2]*vel_LH[2];
}

/* This function convert a WGS-84 ECEF velocity into NED velocity
 * Latitude and longitude are expressed in radians.
 */
void ECEFtoNAV_vel(double vel_ECEF[3], double vel_LH[3], double latitude, double longitude)
{
	float Cen[3][3];

	ECEFtoNAV_mat(Cen, latitude, longitude);

	vel_LH[0] = Cen[0][0]*vel_ECEF[0] + Cen[0][1]*vel_ECEF[1] + Cen[0][2]*vel_ECEF[2];
	vel_LH[1] = Cen[1][0]*vel_ECEF[0] + Cen[1][1]*vel_ECEF[1] + Cen[1][2]*vel_ECEF[2];
	vel_LH[2] = Cen[2][0]*vel_ECEF[0] + Cen[2][1]*vel_ECEF[1] + Cen[2][2]*vel_ECEF[2];
}

#if USE_IMU == 1

/* This function calculates the M and N factors in terms of latitude */
void compute_M_N(double latitude, double *M, double *N)
{
	double sin_lat  = sin(latitude);

	double sqrtden = sqrt(1 - E2_G*sin_lat*sin_lat);

	*M = R_EARTH*(1 - E2_G)/(sqrtden*sqrtden*sqrtden);
	*N = R_EARTH/sqrtden;
}

/* This function calculates the gravity vector expressed in navigation coordinate frame with respect to
 * a user position, expressed also in navigation coordinate frame. It is allowed to activate the correction
 * of additional terms, as Coriolis or centripetal forces, for more accuracy.
 */
#if CORRECT_ACC == 0
void calculate_geod_acc(double geod_pos[3], float geod_acc[3])
#else
void calculate_geod_acc(double geod_pos[3], double geod_vel[3], double M, double N, float geod_acc[3])
#endif
{
	const double a = 6378137.0;
	const double b = 6356752.3141;
	const double gamma_a = 9.7803267715;
	const double gamma_b = 9.8321863685;
	const double f = 1/298.257222101;
	const double kM = 3.986005e14;
	const double omega = 7.292115e-5;
	const double m = omega * omega * a * a * b / kM;

	double cos_lat = cos(geod_pos[0]);
	double sin_lat = sin(geod_pos[0]);
	double h = geod_pos[2];

	geod_acc[0] = -8.08e-9 * h * sin(2*geod_pos[0]);
	geod_acc[1] = 0.0;
	geod_acc[2] = ( (a*gamma_a*cos_lat*cos_lat + b*gamma_b*sin_lat*sin_lat) / sqrt(a*a*cos_lat*cos_lat + b*b*sin_lat*sin_lat) ) *
			(1 - (2/a) * (1 + f + m - 2*f*sin_lat*sin_lat) * h + 3*h*h/(a*a));

#if CORRECT_ACC == 1
	double vN = geod_vel[0];
	double vE = geod_vel[1];
	double vD = geod_vel[2];
	double lat_dot = vN / (M + h);
	double lon_dot = vE / ((N + h) * cos_lat);
	geod_acc[0] += - 2*omega*sin_lat*vE + lat_dot*vD - lon_dot*sin_lat*vE;
	geod_acc[1] +=   2*omega*sin_lat*vN + 2*omega*cos_lat*vD + lon_dot*sin_lat*vN + lon_dot*cos_lat*vD;
	geod_acc[2] += - 2*omega*cos_lat*vE - lon_dot*cos_lat*vE - lat_dot*vN;
#endif
}

/* This function convert a quaternion into its associated attitude matrix. */
void quaternion2mat(float q[4], float Cba[3][3])
{
	float q0q0 = q[0]*q[0];
	float q1q1 = q[1]*q[1];
	float q2q2 = q[2]*q[2];
	float q3q3 = q[3]*q[3];

	float q0q1 = q[0]*q[1];
	float q0q2 = q[0]*q[2];
	float q0q3 = q[0]*q[3];
	float q1q2 = q[1]*q[2];
	float q1q3 = q[1]*q[3];
	float q2q3 = q[2]*q[3];

	Cba[0][0] = q0q0 + q1q1 - q2q2 - q3q3;
	Cba[1][1] = q0q0 - q1q1 + q2q2 - q3q3;
	Cba[2][2] = q0q0 - q1q1 - q2q2 + q3q3;

	Cba[0][1] = 2*(q1q2 + q0q3);
	Cba[0][2] = 2*(q1q3 - q0q2);
	Cba[1][0] = 2*(q1q2 - q0q3);
	Cba[1][2] = 2*(q2q3 + q0q1);
	Cba[2][0] = 2*(q1q3 + q0q2);
	Cba[2][1] = 2*(q2q3 - q0q1);
}

/* This function converts the attitude matrix into its associated quaternion vector.
 * Since float data is used, it is suitable that numerical errors can bring a divergence in
 * the procedure. Thus, a classification of the information is performed before calculate it,
 * in order to apply the most accurate and numerically stable procedure.
 */
void mat2quaternion(float Cba[3][3], float q[4])
{
	float tmp_qdata[4], normq;

	// Calculate the different combinations for the quaternion.
	tmp_qdata[0] = 1 + Cba[0][0] + Cba[1][1] + Cba[2][2];
	tmp_qdata[1] = 1 + Cba[0][0] - Cba[1][1] - Cba[2][2];
	tmp_qdata[2] = 1 - Cba[0][0] + Cba[1][1] - Cba[2][2];
	tmp_qdata[3] = 1 - Cba[0][0] - Cba[1][1] + Cba[2][2];

	// Update the quaternion using the most suitable value. Typically, most cases will match
	// with the first condition of the if().
	if(tmp_qdata[0] > tmp_qdata[1] && tmp_qdata[0] > tmp_qdata[2] && tmp_qdata[0] > tmp_qdata[3])
	{
		q[0] = 0.5*sqrt(tmp_qdata[0]);
		q[1] = 0.25/q[0]*(Cba[1][2] - Cba[2][1]);
		q[2] = 0.25/q[0]*(Cba[2][0] - Cba[0][2]);
		q[3] = 0.25/q[0]*(Cba[0][1] - Cba[1][0]);
	}
	else if(tmp_qdata[1] > tmp_qdata[0] && tmp_qdata[1] > tmp_qdata[2] && tmp_qdata[1] > tmp_qdata[3])
	{
		q[1] = 0.5*sqrt(tmp_qdata[1]);
		q[0] = 0.25/q[1]*(Cba[1][2] - Cba[2][1]);
		q[2] = 0.25/q[1]*(Cba[0][1] + Cba[1][0]);
		q[3] = 0.25/q[1]*(Cba[2][0] + Cba[0][2]);
	}
	else if(tmp_qdata[2] > tmp_qdata[0] && tmp_qdata[2] > tmp_qdata[1] && tmp_qdata[2] > tmp_qdata[3])
	{
		q[2] = 0.5*sqrt(tmp_qdata[2]);
		q[0] = 0.25/q[2]*(Cba[2][0] - Cba[0][2]);
		q[1] = 0.25/q[2]*(Cba[1][0] + Cba[0][1]);
		q[3] = 0.25/q[2]*(Cba[1][2] + Cba[2][1]);
	}
	else
	{
		q[3] = 0.5*sqrt(tmp_qdata[3]);
		q[0] = 0.25/q[3]*(Cba[0][1] - Cba[1][0]);
		q[1] = 0.25/q[3]*(Cba[2][0] + Cba[0][2]);
		q[2] = 0.25/q[3]*(Cba[1][2] + Cba[2][1]);
	}

	// Calculate norm of the quaternion.
	normq = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);

	// Normalize quaternion to assure that it has norm with value 1.
	q[0] = q[0]/normq;
	q[1] = q[1]/normq;
	q[2] = q[2]/normq;
	q[3] = q[3]/normq;
}

#if CORRECT_ANGLE == 1
/* This function calculates the angular velocity of the rotation between inertial and navigation coordinate
 * frame, expressed in navigation coordinate frame.
 */
void calculate_win_n(float win_n[3], double geodpos[3], double geodvel[3])
{
	// win_n = ((lambda_dot + OMEGA_EARTH)*cos(phi), -phi_dot, -((lambda_dot + OMEGA_EARTH)*sin(phi)))^t

	double M, N;
	compute_M_N(geodpos[0], &M, &N);

	double sinphi  = sin(geodpos[0]);
	double cosphi  = cos(geodpos[0]);

	double phi_dot = geodvel[0] /((M + geodpos[2]));		// phi_dot = vn/((M + h))
	double lam_dot = geodvel[1] /((N + geodpos[2])*cosphi); // lam_dot = ve/((N + h)*cos(phi))

	win_n[0] = +(lam_dot + OMEGA_EARTH)*cosphi;
	win_n[1] = -phi_dot;
	win_n[2] = -(lam_dot + OMEGA_EARTH)*sinphi;
}
#endif

/* This function updates the quaternion vector applying the angular velocity from gyr_corr (gyro measurements
 * in body frame corrected with bias) throughout a time period delta_t. In case that the programmer would want
 * to correct the earth rotation, increasing accuracy, position and velocity in navigation coordinate frame
 * shall also be provided.
 */
#if CORRECT_ANGLE == 1
void update_quaternion(float q[4], double geodpos[3], double geodvel[3], float gyr_corr[3], float delta_t)
#else
void update_quaternion(float q[4], float gyr_corr[3], float delta_t)
#endif
{
	float angle_nb[3], norm_angle, angleinv, coshalf, sinhalf, beta_k[3];
	float q_tmp[4], norm_q, PHI_mat[4][4];

#if CORRECT_ANGLE == 1

	// The activation of this mode corrects the earth rotation term in gyro measurements.
	float Cnb[3][3], Cbn[3][3], win_n[3], win_b[3];

	// Calculate transformation matrix (from body to navigation).
	quaternion2mat(q, Cbn);

	// Calculate inverse transformation matrix (from navigation to body).
	Transpose(Cnb, Cbn, 3, 3);

	// Rotation angle from base to ECEF frame in base coordinates.
	// Equation w_ab(b) = w_ib(b) - C_ab*w_ia(a), being a = n, and
	// w_in(n) = ((lamdot+we)*coslat, -phidot, -(lamdot+we)*sinlat)^t.
	//
	// Calculate w_in(n).
	calculate_win_n(win_n, geodpos, geodvel);
	//
	// Calculate w_in(b) = Cnb*w_in(n).
	matVecMul_f(win_b, Cnb, win_n, 3, 3);

	angle_nb[0] = (gyr_corr[0] - win_b[0])*delta_t;
	angle_nb[1] = (gyr_corr[1] - win_b[1])*delta_t;
	angle_nb[2] = (gyr_corr[2] - win_b[2])*delta_t;
#else
	// The effect of win_b can be considered negligible in a first approach.
	angle_nb[0] = gyr_corr[0]*delta_t;
	angle_nb[1] = gyr_corr[1]*delta_t;
	angle_nb[2] = gyr_corr[2]*delta_t;
#endif

	if(Vec3Norm(angle_nb) > 1e-5)
	{
		norm_angle = sqrt(angle_nb[0]*angle_nb[0] + angle_nb[1]*angle_nb[1] + angle_nb[2]*angle_nb[2]);
		angleinv   = 1/norm_angle;
		coshalf    = cos(0.5*norm_angle);
		sinhalf    = sin(0.5*norm_angle);

		beta_k[0] = angleinv*sinhalf*angle_nb[0];
		beta_k[1] = angleinv*sinhalf*angle_nb[1];
		beta_k[2] = angleinv*sinhalf*angle_nb[2];

		// Transition matrix.
		PHI_mat[0][0] = coshalf;
		PHI_mat[1][1] = coshalf;
		PHI_mat[2][2] = coshalf;
		PHI_mat[3][3] = coshalf;

		PHI_mat[0][1] = +beta_k[0];
		PHI_mat[0][2] = +beta_k[1];
		PHI_mat[0][3] = +beta_k[2];
		PHI_mat[1][0] = -beta_k[0];
		PHI_mat[1][2] = +beta_k[2];
		PHI_mat[1][3] = -beta_k[1];
		PHI_mat[2][0] = -beta_k[1];
		PHI_mat[2][1] = -beta_k[2];
		PHI_mat[2][3] = +beta_k[0];
		PHI_mat[3][0] = -beta_k[2];
		PHI_mat[3][1] = +beta_k[1];
		PHI_mat[3][2] = -beta_k[0];

		// Calculate the new quaternion value.
		matVecMul_f(q_tmp, PHI_mat, q, 4, 4);

		// Keep the norm of the quaternion equal to 1.
		norm_q = sqrt(q_tmp[0]*q_tmp[0] + q_tmp[1]*q_tmp[1] + q_tmp[2]*q_tmp[2] + q_tmp[3]*q_tmp[3]);

		// Normalize quaternion.
		q[0] = q_tmp[0]/norm_q;
		q[1] = q_tmp[1]/norm_q;
		q[2] = q_tmp[2]/norm_q;
		q[3] = q_tmp[3]/norm_q;
	}
}

#endif
