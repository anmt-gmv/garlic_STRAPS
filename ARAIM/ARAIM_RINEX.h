/*
 * ARAIM_RINEX.h
 *
 *  Created on: 22 ago. 2019
 *      Author: anmt
 */

#ifndef ARAIM_RINEX_H_
#define ARAIM_RINEX_H_

 /*---------------------------------------------------------------------------------------------------------------
										INCLUDE LIBRARIES
 ---------------------------------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include "dbtypedef.h"

 /*---------------------------------------------------------------------------------------------------------------
										STRUCTS DEFINITIONS
 ---------------------------------------------------------------------------------------------------------------*/

/* Struct with ARAIM configuration parameters */
typedef struct
{
	double Psat_GAL, Psat_GPS, Pconst_GAL, Pconst_GPS, bnom_GAL, bnom_GPS;
	double UERE_usr_GAL_deg[19], UERE_usr_GAL[19];
	double SISE_GAL, ura_GAL, URE_GPS, ura_GPS;
	double TOL_PL, Psat_THRES, Phmi_VERT, Phmi_HOR, Pfa_VERT, Pfa_HOR;
	double KACC, KFF, ACCthres, ACCffthres, Pemt, EMTthres, HAL, VAL, std_v_req, std_h_req;
	double Pfa_CHI2;

} configuration;


/*---------------------------------------------------------------------------------------------------------------
										FUNCTIONS DECLARATIONS
---------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------Functions for hard debugging (folder "hard_debugging") ------------------------------------*/

void disp_matrix( int n_rows, int n_cols, double A[n_rows][n_cols] );

void disp_matrix_float( int n_rows, int n_cols, float A[n_rows][n_cols] );

void disp_vect(int size, double a[size]);

void disp_vect_int(int size, int a[size]);

void disp_matrix_int( int n_rows, int n_cols, int A[n_rows][n_cols] );


/*-----------------------------------Load ARAIM configuration parameters---------------------------------------*/

void load_ARAIM_parameters( configuration* ARAIM );


/*------------------------------  Functions to extract data from RINEX structures --------------------------------------*/

void from_rinex_structures( lsq_state_t *lsq_state, lsq_matrix_t *lsq_mat, gnssdata_t *gnssdata, configuration *ARAIM,
		int Nsat, int Nx, int* lsq_num_LOS, int* PRN_used, double* G, double* w_ura, double* w_ure, double* c_ure, double* pr_res );


/*---------------------------------------  Functions specific for fault mode solutions (folder "fault_mode_solution" ) ---------------------------------------------*/

void ECEF_to_ENU( int Nsat, int Nx, int PRN[Nsat], double G_ECEF[Nsat][Nx], double usrpos[3], double* G_ENU );

/*------------------------------------------------ Matrix/vector operations (folder "matrix_vect")----------------------------------------------------*/

void rank1up( int Nx, int Nsat, double C[Nx][Nx], double G[Nsat][Nx], double wi, int idx_sat, double* C_red );

void rank1up_col( int Nx, double C[Nx][Nx], int idx_col, double* C_up );

void matrix_mul( int P, int Q, int R, double A[P][Q], double B[Q][R], double* C1);

void matrix_diff( int n_rows, int n_cols, double A[n_rows][n_cols], double B[n_rows][n_cols], double* C );

void matrix_sum( int n_rows, int n_cols, double A[n_rows][n_cols], double B[n_rows][n_cols], double* C );

void transpose(int n_rows, int n_cols, double A[n_rows][n_cols], double* At);

void matrix_diff( int n_rows, int n_cols, double A[n_rows][n_cols], double B[n_rows][n_cols], double* C );

int inv( int Nrows, int Ncols, double A[Nrows][Ncols], double *p );

void minor_matrix( int Nrows, int Ncols, double A[Nrows][Ncols], int idx_r, int idx_c, double *p );

double det( int Nrows, int Ncols, double A[Nrows][Ncols] );

void cofactor_matrix( int Nrows, int Ncols, double A[Nrows][Ncols], double *p );

void init_array_double( int nrows, int ncols, double* A_init );

void float_to_double( int init_rows, int init_cols,  int n_rows, int n_cols, float A[init_rows][init_cols], double* A_db );

void copy_array_double( int n_rows, int n_cols, double A[n_rows][n_cols], double* copy );

void copy_vect_double( int n, double A[n], double* copy );

void drop_vector_el_int( int Nsat , int pr_res[Nsat], int idx_sat, int* pr_res_fde );

void copy_vect_int( int n, int A[n], int* copy );

void drop_matrix_row( int n_rows, int n_cols, double A[n_rows][n_cols], int row_to_drop, double* A_drop );

void drop_matrix_col( int n_rows, int n_cols, double A[n_rows][n_cols], int col_to_drop, double* A_drop );

double sum_elements_vect( int size, double a[size] );

void drop_vector_el( int Nsat , double pr_res[Nsat], int idx_sat, double* pr_res_fde );

// new functions
int check_new_sat_in( gnssdata_t *gnssdata, int nsat, int *nsat_previous,
		int sat_iw[MAX_CHANNELS_G], int sat_iw_previous[MAX_CHANNELS_G], int sat_new[MAX_CHANNELS_G] );


#endif /* ARAIM_RINEX_H_ */
