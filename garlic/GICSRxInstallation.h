#include "GICSRx_defines.h"

#if CALC_INSTALL_MAT == 1

#ifndef GICSRXINSTALLATION_H_
#define GICSRXINSTALLATION_H_

typedef struct
{
	int   static_ready;

	float static_acc [3];
	float static_gyr [3];

	int   x_axis_ready;
	int   y_axis_ready;
	int   z_axis_ready;

	float x_axis_t   [3];

	float x_axis     [3];
	float y_axis     [3];
	float z_axis     [3];

	float att_angles [3];

} output_data_t;

void reset_axis_status();

void set_axis_status(float imat[3][3], float att_angles[3]);

void get_axis_status(output_data_t *outcpy);

void initialize_installation();

void calculate_installation_matrix(float rawmeas_acc[][3], float rawmeas_gyr[][3], int nsamples,
		int *imat_ready, float imat[3][3]);
#endif /* GICSRXINSTALLATION_H_ */

#endif
