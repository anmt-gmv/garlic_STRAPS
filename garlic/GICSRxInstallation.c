#include "GICSRxInstallation.h"

#if CALC_INSTALL_MAT == 1

#include <math.h>
#include "algebra.h"

#define vec3norm(A) (sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]))

#define vec3dot(A,B) (A[0]*B[0] + A[1]*B[1] + A[2]*B[2])

#define MIN_STAT_SECS			10
#define MAX_STAT_SECS			30
#define MAX_STAT_AGE			200

#define MIN_TURN_SECS			5
#define MAX_TURN_SECS			30

#define MAX_GYR_AVG_STILL 		0.02
#define MAX_GYR_STD_STILL		0.015
#define MAX_ACC_STD_STILL		0.30

#define MIN_GYR_VERT_COS		0.50
#define MAX_GYR_AVG_DSTILL 		0.002
#define MAX_ACC_AVG_DSTILL		0.02

#define MIN_GYR_AVG_TURN		0.30
#define MAX_GYR_AVG_TURN		0.80
#define MAX_GYR_STD_TURN		0.05
#define MAX_ACC_STD_TURN		1

#define MAXCOUNT_Z_FAULT		2
#define MAX_DIFF_ZLZR			0.10
#define MAXZ_INNOVATION			0.035

#define MINCOUNT_X_VALID		15
#define MAXCOUNT_X_VALID		100

#define MAX_Z_GYR_AVG_LINEAR 	0.015
#define MAX_Z_GYR_STD_LINEAR 	0.01
#define MIN_XY_ACC_AVG_LINEAR	1
#define MAX_Z_ACC_AVG_LINEAR	0.3
#define XY_OVER_Z_FACTOR		15
#define XZ_TOLERANCE			3
#define MIN_ACC_FORWARD_COSINE	0.5

// Used in calculate_static_biases
typedef struct
{
	int data_ready;
	int data_age;
	int counter;

	float sacc_bias[3];
	float sgyr_bias[3];

	float acc_data[MAX_STAT_SECS][3];
	float gyr_data[MAX_STAT_SECS][3];

} still_data_t;

// Used in calculate_z_axis
typedef struct
{
	int   axis_ready;

	float z_tw;
	float x_tw;

	float z_axis[3];
	float x_axis[3];

	int   left_counter;
	int   right_counter;
	int   count_z_fault;

	float zl[MAX_TURN_SECS][3];
	float zr[MAX_TURN_SECS][3];
	float xl[MAX_TURN_SECS][3];
	float xr[MAX_TURN_SECS][3];

} z_axis_data_t;

// Used in calculate_x_axis
typedef struct
{
	int   axis_ready;
	float x_axis[3];

	int   block_ready;
	int   counter;
	float block[MAXCOUNT_X_VALID][3];

} x_axis_data_t;

static still_data_t   still;
static z_axis_data_t  zaxdata;
static x_axis_data_t  xaxdata;
static output_data_t  outdata;

static void vec3cross(float v1[3], float v2[3], float vres[3])
{
	vres[0] = v1[1]*v2[2] - v1[2]*v2[1];
	vres[1] = v1[2]*v2[0] - v1[0]*v2[2];
	vres[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

static void calculate_mean_vector(float *values, int nsamples, float *result)
{
	*result = 0;

	int k;
	for(k = 0; k < nsamples; k++)
	{
		*result += values[k];
	}
	*result /= nsamples;
}

static void calculate_mean_matrix(float values[][3], int nsamples, float result[3])
{
	int j, k;
	for(k = 0; k < 3; k++)
	{
		result[k] = 0;
		for(j = 0; j < nsamples; j++)
		{
			result[k] += values[j][k];
		}
		result[k] /= nsamples;
	}
}

static void calculate_median_matrix(float values[][3], int nsamples, float result[3])
{
	int j, k;
	float tmp_list[MAXCOUNT_X_VALID];

	for(k = 0; k < 3; k++)
	{
		for(j = 0; j < nsamples; j++)
		{
			tmp_list[j] = values[j][k];
		}
		compute_median(tmp_list, nsamples, &result[k]);
	}
}

static void calculate_std_vector(float *values, int nsamples, float *result)
{
	*result = 0;

	float valmean;
	calculate_mean_vector(values, nsamples, &valmean);

	int k;
	for(k = 0; k < nsamples; k++)
	{
		*result += (values[k] - valmean)*(values[k] - valmean);
	}
	*result /= (nsamples-1);

	*result = sqrt(fabs(*result));
}

static void calculate_std_matrix(float values[][3], int nsamples, float result[3])
{
	float valmean[3];
	calculate_mean_matrix(values, nsamples, valmean);

	int j, k;
	for(k = 0; k < 3; k++)
	{
		result[k] = 0;
		for(j = 0; j < nsamples; j++)
		{
			result[k] += (values[j][k] - valmean[k])*(values[j][k] - valmean[k]);
		}
		result[k] /= (nsamples-1);

		result[k] = sqrt(fabs(result[k]));
	}
}

static void normalize_vector(float vector[3])
{
	float norm_value = vec3norm(vector);

	vector[0] = vector[0] / norm_value;
	vector[1] = vector[1] / norm_value;
	vector[2] = vector[2] / norm_value;
}

static void normalize_matrix(float vector[][3], int nsamples)
{
	int k;
	for(k = 0; k < nsamples; k++)
	{
		normalize_vector(vector[k]);
	}
}

static void calculate_static_biases(float rawmeas_acc[][3], float rawmeas_gyr[][3], int nsamples,
		float acc_mean[3], float gyr_mean[3], float acc_std[3], float gyr_std[3])
{
	if(vec3norm(gyr_mean) < MAX_GYR_AVG_STILL && vec3norm(gyr_std) < MAX_GYR_STD_STILL && vec3norm(acc_std) < MAX_ACC_STD_STILL)
	{
		if(still.counter < MAX_STAT_SECS-1) { still.counter++; }

		if(still.counter >= MIN_STAT_SECS)
		{
			float mean_acc_data[3], mean_gyr_data[3];

			calculate_mean_matrix (still.acc_data, still.counter, mean_acc_data);
			calculate_mean_matrix (still.gyr_data, still.counter, mean_gyr_data);

			if(fabs(vec3norm(acc_mean) - vec3norm(mean_acc_data)) < MAX_ACC_AVG_DSTILL &&
					fabs(vec3norm(gyr_mean) - vec3norm(mean_gyr_data)) < MAX_GYR_AVG_DSTILL)
			{
				still.data_ready = 1;
				still.data_age   = 0;
			}
			else
			{
				still.counter = -1;
			}
		}
		if(still.counter >= 0)
		{
			if(still.counter < MAX_STAT_SECS-1)
			{
				int k;
				for(k = 0; k < 3; k++)
				{
					still.acc_data[still.counter][k] = acc_mean[k];
					still.gyr_data[still.counter][k] = gyr_mean[k];
				}
			}
			else
			{
				int j, k;
				for(k = 0; k < 3; k++)
				{
					for(j = 0; j < MAX_STAT_SECS-1; j++)
					{
						still.acc_data[j][k] = still.acc_data[j+1][k];
						still.gyr_data[j][k] = still.gyr_data[j+1][k];
					}
					still.acc_data[MAX_STAT_SECS-1][k] = acc_mean[k];
					still.gyr_data[MAX_STAT_SECS-1][k] = gyr_mean[k];
				}
			}
			calculate_mean_matrix (still.acc_data, still.counter+1, still.sacc_bias);
			calculate_mean_matrix (still.gyr_data, still.counter+1, still.sgyr_bias);
		}
	}
	else
	{
		still.counter = -1;
	}

	if(still.data_ready == 1)
	{
		still.data_age++;
	}

	if(still.data_age > MAX_STAT_AGE)
	{
		still.data_ready = 0;
		still.data_age   = 0;
	}

	outdata.static_ready = still.data_ready;

	int k;
	for(k = 0; k < 3; k++)
	{
		outdata.static_acc[k] = still.sacc_bias[k];
		outdata.static_gyr[k] = still.sgyr_bias[k];
	}
}

static void calculate_z_axis(int turnInProgress, float acc_data[3], float gyr_data[3])
{
	int turnData_av = (turnInProgress > 0) ? 1 : 0;

	if(still.data_ready == 1)
	{
		if(turnInProgress == 1 && zaxdata.left_counter < MAX_TURN_SECS)
		{
			float tmp_vector[MAX_TURN_SECS][3], xvcross[3], xl_axis_std[3], zl_axis_std[3];

			zaxdata.zl[zaxdata.left_counter][0] = -gyr_data[0];
			zaxdata.zl[zaxdata.left_counter][1] = -gyr_data[1];
			zaxdata.zl[zaxdata.left_counter][2] = -gyr_data[2];

			memcpy(tmp_vector, zaxdata.zl, sizeof(tmp_vector));

			normalize_matrix     (tmp_vector, zaxdata.left_counter);
			calculate_std_matrix (tmp_vector, zaxdata.left_counter, zl_axis_std);

			vec3cross(acc_data, gyr_data, xvcross);

			zaxdata.xl[zaxdata.left_counter][0] = xvcross[0];
			zaxdata.xl[zaxdata.left_counter][1] = xvcross[1];
			zaxdata.xl[zaxdata.left_counter][2] = xvcross[2];

			memcpy(tmp_vector, zaxdata.xl, sizeof(tmp_vector));

			normalize_matrix     (tmp_vector, zaxdata.left_counter);
			calculate_std_matrix (tmp_vector, zaxdata.left_counter, xl_axis_std);

			zaxdata.left_counter++;

			if(vec3norm(zl_axis_std) > MAX_GYR_STD_TURN)
			{
				zaxdata.left_counter = 0;
				turnData_av          = 0;
			}
		}
		if(turnInProgress == 2 && zaxdata.right_counter < MAX_TURN_SECS)
		{
			float tmp_vector[MAX_TURN_SECS][3], xvcross[3], xr_axis_std[3], zr_axis_std[3];

			zaxdata.zr[zaxdata.right_counter][0] = +gyr_data[0];
			zaxdata.zr[zaxdata.right_counter][1] = +gyr_data[1];
			zaxdata.zr[zaxdata.right_counter][2] = +gyr_data[2];

			memcpy(tmp_vector, zaxdata.zr, sizeof(tmp_vector));

			normalize_matrix     (tmp_vector, zaxdata.right_counter);
			calculate_std_matrix (tmp_vector, zaxdata.right_counter, zr_axis_std);

			vec3cross(acc_data, gyr_data, xvcross);

			zaxdata.xr[zaxdata.right_counter][0] = xvcross[0];
			zaxdata.xr[zaxdata.right_counter][1] = xvcross[1];
			zaxdata.xr[zaxdata.right_counter][2] = xvcross[2];

			memcpy(tmp_vector, zaxdata.xr, sizeof(tmp_vector));

			normalize_matrix     (tmp_vector, zaxdata.right_counter);
			calculate_std_matrix (tmp_vector, zaxdata.right_counter, xr_axis_std);

			zaxdata.right_counter++;

			if(vec3norm(zr_axis_std) > MAX_GYR_STD_TURN)
			{
				zaxdata.right_counter = 0;
				turnData_av           = 0;
			}
		}
	}

	if(turnData_av == 1 && zaxdata.left_counter >= MIN_TURN_SECS && zaxdata.right_counter >= MIN_TURN_SECS)
	{
		float tmp_vector[MAX_TURN_SECS][3];
		float zl_axis_mean[3], zr_axis_mean[3], xl_axis_mean[3], xr_axis_mean[3];

		memcpy(tmp_vector, zaxdata.zl, sizeof(tmp_vector));
		normalize_matrix      (tmp_vector,   zaxdata.left_counter);
		calculate_mean_matrix (tmp_vector,   zaxdata.left_counter, zl_axis_mean);
		normalize_vector      (zl_axis_mean);

		memcpy(tmp_vector, zaxdata.zr, sizeof(tmp_vector));
		normalize_matrix      (tmp_vector,   zaxdata.right_counter);
		calculate_mean_matrix (tmp_vector,   zaxdata.right_counter, zr_axis_mean);
		normalize_vector      (zr_axis_mean);

		memcpy(tmp_vector, zaxdata.xl, sizeof(tmp_vector));
		normalize_matrix      (tmp_vector,   zaxdata.left_counter);
		calculate_mean_matrix (tmp_vector,   zaxdata.left_counter, xl_axis_mean);
		normalize_vector      (xl_axis_mean);

		memcpy(tmp_vector, zaxdata.xr, sizeof(tmp_vector));
		normalize_matrix      (tmp_vector,   zaxdata.right_counter);
		calculate_mean_matrix (tmp_vector,   zaxdata.right_counter, xr_axis_mean);
		normalize_vector      (xr_axis_mean);

		float diff_zl_zr[3] = {zl_axis_mean[0] - zr_axis_mean[0],
				zl_axis_mean[1] - zr_axis_mean[1],
				zl_axis_mean[2] - zr_axis_mean[2]};

		if(vec3norm(diff_zl_zr) < MAX_DIFF_ZLZR)
		{
			float z_lw, z_rw, z_tw_aux, z_axis_aux[3];
			float x_lw, x_rw, x_tw_aux, x_axis_aux[3];

			z_lw = 1/(vec3dot(zl_axis_mean, zl_axis_mean));
			z_rw = 1/(vec3dot(zr_axis_mean, zr_axis_mean));

			z_axis_aux[0] = z_lw*zl_axis_mean[0] + z_rw*zr_axis_mean[0];
			z_axis_aux[1] = z_lw*zl_axis_mean[1] + z_rw*zr_axis_mean[1];
			z_axis_aux[2] = z_lw*zl_axis_mean[2] + z_rw*zr_axis_mean[2];

			x_lw = 1/(vec3dot(xl_axis_mean, xl_axis_mean));
			x_rw = 1/(vec3dot(xr_axis_mean, xr_axis_mean));

			x_axis_aux[0] = x_lw*xl_axis_mean[0] + x_rw*xr_axis_mean[0];
			x_axis_aux[1] = x_lw*xl_axis_mean[1] + x_rw*xr_axis_mean[1];
			x_axis_aux[2] = x_lw*xl_axis_mean[2] + x_rw*xr_axis_mean[2];

			z_tw_aux = z_lw + z_rw;
			x_tw_aux = x_lw + x_rw;

			normalize_vector (z_axis_aux);
			normalize_vector (x_axis_aux);

			if(zaxdata.axis_ready == 0)
			{
				memcpy(zaxdata.z_axis, z_axis_aux, 3*sizeof(float));
				memcpy(zaxdata.x_axis, x_axis_aux, 3*sizeof(float));

				zaxdata.z_tw = z_tw_aux;
				zaxdata.x_tw = x_tw_aux;

				zaxdata.axis_ready    = 1;
				zaxdata.count_z_fault = 0;
			}
			else
			{
				float diffz_innovation[3] = {z_axis_aux[0] - zaxdata.z_axis[0],
						z_axis_aux[1] - zaxdata.z_axis[1],
						z_axis_aux[2] - zaxdata.z_axis[2]};

				if(vec3norm(diffz_innovation) < MAXZ_INNOVATION)
				{
					zaxdata.z_tw = zaxdata.z_tw + z_tw_aux;
					zaxdata.x_tw = zaxdata.x_tw + x_tw_aux;

					int k;
					for(k = 0; k < 3; k++)
					{
						zaxdata.z_axis[k] = z_tw_aux*z_axis_aux[k] + zaxdata.z_tw*zaxdata.z_axis[k];
						zaxdata.x_axis[k] = x_tw_aux*x_axis_aux[k] + zaxdata.x_tw*zaxdata.x_axis[k];
					}
					normalize_vector(zaxdata.z_axis);
					normalize_vector(zaxdata.x_axis);

					zaxdata.count_z_fault = 0;
				}
				else
				{
					zaxdata.count_z_fault++;
				}

				if(zaxdata.count_z_fault >= MAXCOUNT_Z_FAULT)
				{
					zaxdata.axis_ready    = 0;
					zaxdata.left_counter  = 0;
					zaxdata.right_counter = 0;
				}
			}
		}
		else
		{
			zaxdata.axis_ready    = 0;
			zaxdata.left_counter  = 0;
			zaxdata.right_counter = 0;
		}
	}

	outdata.z_axis_ready = zaxdata.axis_ready;

	int k;
	for(k = 0; k < 3; k++)
	{
		outdata.z_axis   [k] = zaxdata.z_axis[k];
		outdata.x_axis_t [k] = zaxdata.x_axis[k];
	}
}

static void calculate_x_axis(float acc_data[3], float gyr_data[3], float std_zgyr,
		int reset_x_comp)
{
	float x_axis_local[3], z_axis_local[3];

	memcpy(z_axis_local, outdata.z_axis, 3*sizeof(float));

	if(reset_x_comp == 1)
	{
		xaxdata.counter     = 0;
		xaxdata.axis_ready  = 0;
		xaxdata.block_ready = 0;
	}

	if(still.data_ready == 1 && xaxdata.counter < MAXCOUNT_X_VALID)
	{
		if(outdata.z_axis_ready == 0)
		{
			z_axis_local[0] = still.sacc_bias[0];
			z_axis_local[1] = still.sacc_bias[1];
			z_axis_local[2] = still.sacc_bias[2];

			normalize_vector(z_axis_local);
		}

		float z_omega  = vec3dot(gyr_data, z_axis_local);
		float z_accel  = vec3dot(acc_data, z_axis_local);

		float xy_accel[3];
		xy_accel[0] = acc_data[0] - z_accel*z_axis_local[0];
		xy_accel[1] = acc_data[1] - z_accel*z_axis_local[1];
		xy_accel[2] = acc_data[2] - z_accel*z_axis_local[2];

		if(std_zgyr < MAX_Z_GYR_STD_LINEAR && fabs(z_omega) < MAX_Z_GYR_AVG_LINEAR && vec3norm(xy_accel) > MIN_XY_ACC_AVG_LINEAR &&
				fabs(z_accel) < MAX_Z_ACC_AVG_LINEAR && vec3norm(xy_accel) > XY_OVER_Z_FACTOR*fabs(z_accel))
		{
			if(outdata.z_axis_ready == 0)
			{
				if(xaxdata.block_ready == 1)
				{
					calculate_mean_matrix (xaxdata.block, xaxdata.counter, x_axis_local);
					normalize_vector      (x_axis_local);
				}
				else
				{
					memcpy(x_axis_local, xy_accel, 3*sizeof(float));
					normalize_vector(x_axis_local);
				}
			}
			else
			{
				float x_block_mean[3];
				calculate_mean_matrix (xaxdata.block, xaxdata.counter, x_block_mean);

				memcpy(x_axis_local, outdata.x_axis_t, 3*sizeof(float));
				if(vec3dot(x_block_mean, x_axis_local) < 0)
				{
					int j, k;
					for(j = 0; j < xaxdata.counter; j++)
					{
						for(k = 0; k < 3; k++)
						{
							xaxdata.block[j][k] = -xaxdata.block[j][k];
						}
					}
				}
			}
			float acc_sign = vec3dot(x_axis_local, xy_accel)/vec3norm(xy_accel);

			if(fabs(acc_sign) > +MIN_ACC_FORWARD_COSINE)
			{
				int xy_sign = (acc_sign < 0 ? -1 : +1);

				xaxdata.block[xaxdata.counter][0] = xy_sign*xy_accel[0];
				xaxdata.block[xaxdata.counter][1] = xy_sign*xy_accel[1];
				xaxdata.block[xaxdata.counter][2] = xy_sign*xy_accel[2];

				xaxdata.counter++;
				xaxdata.block_ready = 1;
			}

			if(xaxdata.counter >= MINCOUNT_X_VALID)
			{
				float x_block_tmp[MAXCOUNT_X_VALID][3];
				memcpy(x_block_tmp, xaxdata.block, sizeof(x_block_tmp));

				normalize_matrix        (x_block_tmp, xaxdata.counter);
				//calculate_mean_matrix (x_block_tmp, xaxdata.counter, xaxdata.x_axis);
				calculate_median_matrix (x_block_tmp, xaxdata.counter, xaxdata.x_axis);
				normalize_vector        (xaxdata.x_axis);

				xaxdata.axis_ready = 1;
			}
		}
	}
	outdata.x_axis_ready = xaxdata.axis_ready;

	outdata.x_axis[0] = xaxdata.x_axis[0];
	outdata.x_axis[1] = xaxdata.x_axis[1];
	outdata.x_axis[2] = xaxdata.x_axis[2];
}

static void calculate_zgyr_std(float rawmeas_gyr[][3], int nsamples, float *std_zgyr)
{
	float z_axis_local [3];
	float omega_z_vec  [MAX_IMU_SAMPLES];

	float *z_axis_tmp = (outdata.z_axis_ready == 0 ? still.sacc_bias : outdata.z_axis);

	memcpy(z_axis_local, z_axis_tmp, 3*sizeof(float));
	normalize_vector(z_axis_local);

	int k;
	for(k = 0; k < nsamples; k++)
	{
		omega_z_vec[k] = vec3dot(rawmeas_gyr[k], z_axis_local);
	}
	calculate_std_vector(omega_z_vec, nsamples, std_zgyr);
}

void initialize_installation()
{
	memset(&still,   0, sizeof(still_data_t));
	memset(&zaxdata, 0, sizeof(z_axis_data_t));
	memset(&xaxdata, 0, sizeof(x_axis_data_t));
	memset(&outdata, 0, sizeof(output_data_t));

	still.counter = -1;
}

static void calculate_y_axis(int *reset_x_comp, int *imat_ready, float imat[3][3]/*, float att_angles[3]*/)
{
	*reset_x_comp = 0;

	if(outdata.x_axis_ready == 1 && outdata.z_axis_ready == 1)
	{
		if(vec3dot(outdata.x_axis_t, outdata.x_axis) < 0)
		{
			outdata.x_axis[0] = -outdata.x_axis[0];
			outdata.x_axis[1] = -outdata.x_axis[1];
			outdata.x_axis[2] = -outdata.x_axis[2];
		}

		if(fabs(asin(vec3dot(outdata.x_axis, outdata.z_axis))*RADTODEG < XZ_TOLERANCE))
		{
			float xzproj = vec3dot(outdata.x_axis, outdata.z_axis);

			outdata.x_axis[0] = outdata.x_axis[0] - xzproj*outdata.z_axis[0];
			outdata.x_axis[1] = outdata.x_axis[1] - xzproj*outdata.z_axis[1];
			outdata.x_axis[2] = outdata.x_axis[2] - xzproj*outdata.z_axis[2];

			normalize_vector(outdata.x_axis);

			vec3cross(outdata.z_axis, outdata.x_axis, outdata.y_axis);

			outdata.y_axis_ready = 1;
		}
		else
		{
			*reset_x_comp = 1;
		}

		if(outdata.y_axis_ready == 1)
		{
			*imat_ready = 1;

			imat[0][0]  = outdata.x_axis[0];
			imat[0][1]  = outdata.x_axis[1];
			imat[0][2]  = outdata.x_axis[2];

			imat[1][0]  = outdata.y_axis[0];
			imat[1][1]  = outdata.y_axis[1];
			imat[1][2]  = outdata.y_axis[2];

			imat[2][0]  = outdata.z_axis[0];
			imat[2][1]  = outdata.z_axis[1];
			imat[2][2]  = outdata.z_axis[2];

			float tmpvec_1[3] = {imat[1][0], -imat[0][0], 0};
			float tmpvec_2[3] = {imat[0][0],  imat[1][0], imat[2][0]};
			float tmpvec_3[3] = {imat[0][2],  imat[1][2], imat[2][2]};
			float tmpvec_4[3];

			vec3cross(tmpvec_1, tmpvec_2, tmpvec_4);

			float p1 = vec3dot(tmpvec_1, tmpvec_3);
			float p2 = vec3dot(tmpvec_4, tmpvec_3);

			outdata.att_angles[0] = atan2(imat[1][0], imat[0][0])*RADTODEG;
			outdata.att_angles[1] = -asin(imat[2][0])*RADTODEG;
			outdata.att_angles[2] = atan2(p1, p2)*RADTODEG;

			fprintf(stderr, "installation matrix calculated: %.4f,%.4f,%.4f\n",outdata.att_angles[0],outdata.att_angles[1],outdata.att_angles[2]);
		}
	}
}

static void correct_static_biases(float rawmeas_acc[][3], float rawmeas_gyr[][3], int nsamples)
{
	// Correct static biases.
	int j, k;
	for(k = 0; k < 3; k++)
	{
		for(j = 0; j < nsamples; j++)
		{
			rawmeas_acc[j][k] -= still.sacc_bias[k];
			rawmeas_gyr[j][k] -= still.sgyr_bias[k];
		}
	}
}

void reset_axis_status()
{
	memset(&still,   0, sizeof(still_data_t));
	memset(&zaxdata, 0, sizeof(z_axis_data_t));
	memset(&xaxdata, 0, sizeof(x_axis_data_t));
	memset(&outdata, 0, sizeof(output_data_t));

	still.counter = -1;
}

void set_axis_status(float imat[3][3], float att_angles[3])
{
	outdata.static_ready = 1;

	outdata.static_acc[0] = 0;
	outdata.static_acc[1] = 0;
	outdata.static_acc[2] = -GRAVITY;

	outdata.static_gyr[0] = 0;
	outdata.static_gyr[1] = 0;
	outdata.static_gyr[2] = 0;

	outdata.x_axis_ready  = 1;
	outdata.y_axis_ready  = 1;
	outdata.z_axis_ready  = 1;

	outdata.x_axis[0] = outdata.x_axis_t[0] = imat[0][0];
	outdata.x_axis[1] = outdata.x_axis_t[1] = imat[0][1];
	outdata.x_axis[2] = outdata.x_axis_t[2] = imat[0][2];

	outdata.y_axis[0] = imat[1][0];
	outdata.y_axis[1] = imat[1][1];
	outdata.y_axis[2] = imat[1][2];

	outdata.z_axis[0] = imat[2][0];
	outdata.z_axis[1] = imat[2][1];
	outdata.z_axis[2] = imat[2][2];

	outdata.att_angles[0] = att_angles[0];
	outdata.att_angles[1] = att_angles[1];
	outdata.att_angles[2] = att_angles[2];
}

void get_axis_status(output_data_t *outcpy)
{
	memcpy(outcpy, &outdata, sizeof(output_data_t));
}

void calculate_installation_matrix(float rawmeas_acc[][3], float rawmeas_gyr[][3], int nsamples,
		int *imat_ready, float imat[3][3])
{
	static int reset_x_comp = 0;

	float acc_mean[3], gyr_mean[3], acc_std[3], gyr_std[3];

	calculate_mean_matrix (rawmeas_acc, nsamples, acc_mean);
	calculate_mean_matrix (rawmeas_gyr, nsamples, gyr_mean);
	calculate_std_matrix  (rawmeas_acc, nsamples, acc_std);
	calculate_std_matrix  (rawmeas_gyr, nsamples, gyr_std);

	calculate_static_biases(rawmeas_acc, rawmeas_gyr, nsamples, acc_mean, gyr_mean, acc_std, gyr_std);

	int turnInProgress = 0;
	if(vec3norm(gyr_mean) > MIN_GYR_AVG_TURN && vec3norm(gyr_mean) < MAX_GYR_AVG_TURN &&
			vec3norm(gyr_std)  < MAX_GYR_STD_TURN && vec3norm(acc_std)  < MAX_ACC_STD_TURN)
	{
		if(fabs(vec3dot(acc_mean, gyr_mean)) > +MIN_GYR_VERT_COS*vec3norm(gyr_mean)*vec3norm(acc_mean))
		{
			turnInProgress = vec3dot(acc_mean, gyr_mean) > 0 ? 1 : 2;
		}
	}

	correct_static_biases (rawmeas_acc, rawmeas_gyr, nsamples);

	calculate_mean_matrix (rawmeas_acc, nsamples, acc_mean);
	calculate_mean_matrix (rawmeas_gyr, nsamples, gyr_mean);

	calculate_z_axis(turnInProgress, acc_mean, gyr_mean);

	float std_zgyr = 0.0;
	calculate_zgyr_std(rawmeas_gyr, nsamples, &std_zgyr);

	calculate_x_axis(acc_mean, gyr_mean, std_zgyr, reset_x_comp);

	calculate_y_axis(&reset_x_comp, imat_ready, imat);
}

#endif
