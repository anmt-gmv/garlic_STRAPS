/////////////////////////////////////////////////////////////////////////////////
// Project:                     GICSRX
// Purpose:                     Functions for IBPL computation.
// File:                        GICSRxIBPL.c
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
// DESCRIPTION: this file includes the definitions of the least-squares function,
//              atmospheric correction functions and IBPL integrity functions.
//
//
/////////////////////////////////////////////////////////////////////////////////
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"
#include "algebra.h"
#include "GICSRxIBPL.h"
#include "GICSRxObs.h"
#include "GICSRxWNavSol.h"
#include "GICSRxKalman.h"

static double PLtable_2D[96][7] = {{9.949874,      99.995  ,    999.9995 , 9999.99995  ,     100000   ,   1000000  ,   10000000}   ,
        {3,             9.949874 ,          31.606961  ,        99.995,      316.226185     ,         999.9995   ,         3162.277502},
        {1.908295,      4.532587 ,   9.949874        ,    21.521126     ,     46.405115    ,          99.995  ,    215.441148},
        {1.470469,      3      ,     5.533785  ,   9.949874     ,       17.754655    ,          31.606961     ,     56.22524},
        {1.229588,      2.304251 ,   3.853431   ,         6.229825   ,         9.949874  ,            15.817353     ,     25.098951},
        {1.074446,      1.908295  ,  3   ,        4.532587   ,          6.739131         ,     9.949874       ,     14.643888},
        {0.964727,      1.651543 ,   2.489349  ,          3.590955     ,       5.082023    ,          7.127043       ,     9.949874},
        {0.882201,      1.470469 ,   2.150212   ,         3     ,     4.096681          ,    5.533785     ,       7.431967},
        {0.817374172,     1.335125238  ,   1.908294745 ,    2.596658781  ,            3.451883059  ,   4.532587219 ,    5.910849062},
        {0.764783102 ,    1.229587911  ,   1.726578033,     2.304251168  ,   3  ,         3.853431189,     4.91109604},
        {0.7210486  ,   1.144608973   ,       1.584674015 ,    2.082474305  ,            2.666703641,     3.365778869 ,    4.211671192},
        {0.683958528 ,    1.074446225  ,   1.470468517    , 1.908294745    ,          2.410999936 ,    3   ,        3.698377033},
        {0.65199898  ,      1.015341135,     1.376323408  ,   1.767661275  ,            2.208623117 ,    2.716022165 ,    3.307229417},
        {0.624095741 ,    0.964726764  ,   1.297187648    , 1.65154283   ,           2.044376355  ,   2.489348656  ,   3},
        {0.599463419 ,    0.920787596  ,   1.229587911    , 1.553881873  ,            1.908294745 ,    2.304251168 ,    2.752627651},
        {0.577513145 ,    0.882201457  ,   1.171056662    , 1.470468517  ,            1.793589985 ,    2.150212374  ,   2.549302276},
        {0.557793813 ,    0.847981251  ,   1.119791902    , 1.398287966       ,       1.695486691  ,   2.019954962 ,    2.379253305},
        {0.539953392  ,   0.817374172  ,   1.074446225    , 1.335125238       ,       1.610532106  ,   1.908294745 ,    2.234914429},
        {0.523712694  ,   0.789795378  ,   1.033991335    , 1.279316575       ,       1.536170006  ,   1.811444837  ,   2.110827985},
		{0.50884714   ,     0.764783102,     0.997628345  ,   1.229587911     ,         1.470468517 ,    1.726578033,     2.002965885},
        {0.495173792  ,   0.741967506  ,   0.964726764    , 1.184946912       ,       1.411940967  ,   1.65154283   ,    1.908294745},
        {0.482541956  ,   0.7210486    ,      0.934782019 ,    1.144608973    ,          1.35942483 ,      1.584674015 ,    1.824489321} ,
        {0.470826241  ,   0.701780268  ,   0.907385255    , 1.107945374       ,       1.311998006 ,    1.524663219  ,   1.749738834},
        {0.459921361 ,    0.683958528  ,   0.882201457   ,  1.074446225       ,       1.268919704 ,    1.470468517  ,   1.682613102},
        {0.449738185 ,    0.667412744  ,   0.858953333   ,  1.043693504       ,       1.229587911  ,   1.421250056  ,   1.621968109},
        {0.440200683 ,    0.65199898   ,    0.837409267  ,   1.015341135      ,        1.193508281 ,    1.376323408 ,    1.5668782},
        {0.431243565 ,    0.637594889  ,   0.817374172   ,  0.989100037       ,       1.160271024  ,   1.335125238  ,   1.516586559},
        {0.422810401 ,    0.624095741  ,   0.798682482   ,  0.964726764       ,       1.129533499  ,   1.297187648 ,    1.470468517},
        {0.414852139 ,    0.611411315  ,   0.781192696   ,  0.942014767       ,       1.101006944  ,   1.262118773 ,    1.428003985},
        {0.407325915  ,   0.599463419  ,   0.764783102   ,  0.920787596     ,         1.074446225  ,   1.229587911  ,   1.388756481},
        {0.40019408   ,     0.58818392 ,      0.749348388,     0.900893553,              1.049641849,     1.199314012 ,    1.352356992},
        {0.393423417  ,   0.577513145  ,   0.734796928   ,  0.882201457  ,            1.026413672,     1.171056662  ,   1.318491424} ,
        {0.38698449   ,     0.56739858 ,      0.7210486  ,        0.864597247 ,             1.004605894,     1.144608973  ,   1.286890742},
        {0.380851113 ,    0.557793813  ,   0.708033      ,      0.847981251   ,           0.984083049,     1.119791902 ,    1.257323165},
        {0.374999898 ,    0.54865765   ,    0.695687997  ,   0.832265952      ,        0.964726764 ,    1.096449702 ,    1.229587911},
        {0.369409889 ,    0.539953392  ,   0.683958528   ,  0.817374172       ,       0.946433121 ,    1.074446225 ,    1.203510177},
        {0.364062241 ,    0.531648221  ,   0.672795605  ,   0.803237566       ,       0.929110495 ,    1.053661914 ,    1.178937049},
        {0.358939956 ,    0.523712694  ,   0.662155486  ,   0.789795378       ,       0.912677769 ,    1.033991335 ,    1.155734169},
        {0.354027658 ,    0.516120309  ,   0.65199898   ,    0.7769934        ,      0.897062856 ,    1.015341135  ,   1.133782988},
        {0.3493114   ,  0.50884714     ,   0.642290857  ,   0.764783102   ,           0.882201457  ,   0.997628345  ,   1.112978499},
        {0.344778496 ,    0.501871523 ,    0.632999351  ,   0.753120894       ,       0.86803603 ,      0.980778963 ,    1.093227342},
        {0.340417380  ,   0.495173792 ,    0.624095741  ,   0.741967506       ,       0.854514914 ,    0.964726764  ,   1.074446225},
        {0.336217479  ,   0.488736047 ,    0.615553983  ,   0.731287457       ,       0.841591584 ,    0.949412294  ,   1.056560594},
        {0.332169111  ,   0.482541956 ,    0.607350401  ,   0.721048600       ,       0.829224020 ,    0.934782019  ,   1.039503506},
        {0.328263386  ,   0.476576583 ,    0.599463419  ,   0.711221734       ,       0.817374172 ,    0.920787596  ,   1.023214675},
        {0.324492128  ,   0.470826241 ,    0.591873329  ,   0.701780268   ,           0.806007490  ,   0.907385255 ,    1.007639661},
        {0.320847802  ,   0.465278356 ,    0.584562089  ,   0.692699934   ,           0.795092528  ,   0.894535261 ,    0.992729166},
        {0.317323449   ,  0.459921361 ,    0.577513145  ,   0.683958528   ,           0.784600598  ,   0.882201457  ,   0.978438440},
        {0.313912634   ,  0.454744588 ,    0.570711280  ,   0.675535697   ,           0.774505468  ,   0.870350863  ,   0.964726764},
        {0.310609395 ,    0.449738185 ,    0.564142481  ,   0.667412744   ,           0.764783102  ,   0.858953333 ,    0.951556997},
        {0.307408197 ,    0.444893033 ,    0.557793813  ,   0.659572461   ,           0.755411431  ,   0.847981251   ,  0.938895191},
        {0.304303896 ,    0.440200683 ,    0.551653321 ,    0.651998980   ,           0.746370153 ,    0.837409267  ,   0.926710252},
        {0.301291705 ,    0.435653292 ,    0.545709936 ,    0.644677644   ,           0.737640556 ,    0.827214067  ,   0.914973640},
        {0.298367160   ,  0.431243565  ,   0.539953392 ,    0.637594889   ,           0.729205365 ,    0.817374172  ,   0.903659111},
        {0.295526095  ,   0.426964712  ,   0.534374150 ,    0.630738144   ,           0.721048600   ,  0.807869754  ,   0.892742487},
        {0.292764614,     0.422810401  ,   0.528963339 ,    0.624095741   ,           0.713155459,     0.798682482,     0.882201457},
        {0.290079074 ,    0.418774719  ,   0.523712694 ,    0.617656832   ,           0.705512209,     0.789795378,     0.872015392},
        {0.287466057  ,   0.414852139  ,   0.518614503 ,    0.611411315   ,           0.698106086,     0.781192696,     0.862165194},
        {0.284922360   ,  0.411037488  ,   0.513661565 ,    0.605349774   ,           0.690925215,     0.772859806,     0.852633151},
        {0.282444972 ,    0.407325915  ,   0.508847140 ,    0.599463419   ,           0.683958528,     0.764783102,     0.843402813},
        {0.280031064  ,   0.403712871  ,   0.504164920 ,    0.593744032   ,           0.677195697,     0.756949904,     0.834458877},
        {0.277677972 ,    0.400194080  ,   0.499608988 ,    0.588183920   ,           0.670627073,     0.749348388,     0.825787094},
        {0.275383185 ,    0.396765523  ,   0.495173792 ,    0.582775878   ,           0.664243628,     0.741967506,     0.817374172},
        {0.273144336 ,    0.393423417  ,   0.490854114 ,    0.577513145   ,           0.658036906,     0.734796928,     0.809207699},
        {0.270959190 ,    0.390164193  ,   0.486645042 ,    0.572389369   ,           0.651998980,     0.727826982,     0.801276069},
        {0.268825635  ,   0.386984490  ,   0.482541956 ,    0.567398580   ,           0.646122406,     0.721048600,     0.793568418},
        {0.266741673  ,   0.383881130,     0.478540497 ,    0.562535159   ,           0.640400187,     0.714453270,     0.786074564} ,
        {0.264705413  ,   0.380851113,     0.474636553 ,    0.557793813   ,           0.634825742,     0.708033000,     0.778784954} ,
        {0.262715064  ,   0.377891598,     0.470826241 ,    0.553169547   ,           0.629392875,     0.701780268,     0.771690612} ,
        {0.260768926  ,   0.374999898,     0.467105890 ,    0.548657650   ,           0.624095741,     0.695687997,     0.764783102} ,
        {0.258865387  ,   0.372173467,     0.463472026 ,    0.544253668   ,           0.618928829,     0.689749514,     0.758054481} ,
        {0.257002916  ,   0.369409889,     0.459921361 ,    0.539953392   ,           0.613886932,     0.683958528,     0.751497267} ,
        {0.255180059  ,   0.366706873,     0.456450779 ,    0.535752834   ,           0.608965127,     0.678309096,     0.745104403} ,
        {0.253395430  ,   0.364062241,     0.453057322 ,    0.531648221   ,           0.604158759,     0.672795605,     0.738869228} ,
        {0.251647713  ,   0.361473924,     0.449738185 ,    0.527635972   ,           0.599463419,     0.667412744,     0.732785450} ,
        {0.249935654  ,   0.358939956,     0.446490701 ,    0.523712694   ,           0.594874928,     0.662155486,     0.726847118} ,
        {0.248258055  ,   0.356458461,     0.443312338 ,    0.519875161   ,           0.590389323,     0.657019070,     0.721048600} ,
        {0.246613777  ,   0.354027658,     0.440200683 ,    0.516120309   ,           0.586002844,     0.651998980,     0.715384561} ,
        {0.245001732   ,  0.351645845,     0.437153444 ,    0.512445226   ,           0.581711916,     0.647090932,     0.709849947} ,
        {0.243420881   ,  0.349311400 ,    0.434168432 ,    0.508847140   ,           0.577513145,     0.642290857,     0.704439959} ,
        {0.241870230   ,  0.347022777 ,    0.431243565 ,    0.505323411   ,           0.573403299,     0.637594889,     0.699150045} ,
        {0.240348830  ,   0.344778496 ,    0.428376855 ,    0.501871523   ,           0.569379304,     0.632999351,     0.693975877} ,
        {0.238855774  ,   0.342577147 ,    0.425566404 ,    0.498489080   ,           0.565438231,     0.628500746,     0.688913342} ,
        {0.237390193  ,   0.340417380 ,    0.422810401 ,    0.495173792   ,           0.561577288,     0.624095741,     0.683958528} ,
        {0.235951254  ,   0.338297902 ,    0.420107114 ,    0.491923477   ,           0.557793813,     0.619781163,     0.679107707} ,
        {0.234538160  ,   0.336217479 ,    0.417454888 ,    0.488736047   ,           0.554085264,     0.615553983,     0.674357332} ,
        {0.233150147   ,  0.334174927 ,    0.414852139 ,    0.485609509   ,           0.550449215,     0.611411315,     0.669704017} ,
        {0.231786481   ,  0.332169111 ,    0.412297352 ,    0.482541956   ,           0.546883349,     0.607350401,     0.665144537} ,
        {0.230446460   ,  0.330198945 ,    0.409789074 ,    0.479531563   ,           0.543385447,     0.603368608,     0.660675810} ,
        {0.229129407   ,  0.328263386 ,    0.407325915 ,    0.476576583   ,           0.539953392,     0.599463419,     0.656294895} ,
        {0.227834675   ,  0.326361434 ,    0.404906541 ,    0.473675344   ,           0.536585152,     0.595632427,     0.651998980} ,
        {0.226561641 ,    0.324492128 ,    0.402529671 ,    0.470826241  ,            0.533278786,     0.591873329,     0.647785378} ,
        {0.225309704 ,    0.322654546  ,   0.400194080 ,    0.468027737  ,            0.530032432,     0.588183920,     0.643651518} ,
        {0.224078289 ,    0.320847802  ,   0.397898588 ,    0.465278356  ,            0.526844305,     0.584562089,     0.639594939},
        {0.222866841 ,    0.319071042  ,   0.395642063 ,    0.462576685  ,            0.523712694,     0.581005810,     0.635613283},
        {0.221674827 ,    0.317323449  ,   0.393423417 ,    0.459921361  ,            0.520635957,     0.577513145,     0.631704293}};

static void compute_nav_mat(const double *position, float Rmat[3][3])
{
	double GeoPos[3];

	ECEFtoNAV_pos(position, GeoPos);
	ECEFtoNAV_mat(Rmat, GeoPos[0], GeoPos[1]);
}

void compute_dop(const double *position, int num_pars, const float *cov_mat, float *hdop, float *vdop)
{
	// Definition of local variables.
	int i, j;

	float d2x, d2y, d2z, dxy;
	float Rmat[3][3];
	float Rmat_X[3][KNSTATES];
	float Rmat_X_DOP[3][KNSTATES];

	compute_nav_mat(position, Rmat);

	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
		{
			Rmat_X[i][j] = Rmat[i][j];
		}

		for(; j < KNSTATES; j++) {
			Rmat_X[i][j] = 0.0;
		}
	}

	// Calculate R'*inv(H'*H)*R.
	for (i = 0; i < 3; i++)
	{
		matSymVecMul_f(Rmat_X_DOP[i], cov_mat, Rmat_X[i], num_pars);
	}

	// Extract horizontal variance and covariance values from the
	// state vector error covariance (weighted DoP) matrix.
	DotProduct(Rmat_X[0], Rmat_X_DOP[0], num_pars, d2x);
	DotProduct(Rmat_X[1], Rmat_X_DOP[1], num_pars, d2y);
	DotProduct(Rmat_X[2], Rmat_X_DOP[2], num_pars, d2z);
	DotProduct(Rmat_X[0], Rmat_X_DOP[1], num_pars, dxy);

	// Compute normalized (by the dimension) horizontal and vertical DoP
	*hdop = sqrt((d2x+d2y)/2.0 + sqrt(((d2x-d2y)*(d2x-d2y))/4 + dxy*dxy));
	*vdop = sqrt(d2z);
}

#if CALCULATE_IBPL == 1

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

static float update_qhwh_trace(unsigned int size, const float * h_vec, float sigma_meas, const float * PH)
{
	float dot_prod;
	DotProduct(h_vec, PH, size, dot_prod);
	dot_prod /= (sigma_meas*sigma_meas);

	return dot_prod;

} // END of function update_qhwh_trace

static float update_qhwh_trace_pp(const float * h_vec, float sigma_meas, const float * PH)
{
	float dot_prod =
			PH[0] * h_vec[0] +
			PH[1] * h_vec[1] +
			PH[2] * h_vec[2];
	dot_prod /= (sigma_meas*sigma_meas);

	return dot_prod;

} // END of function update_qhwh_trace

static float update_qhwh_trace_pv(const float * h_vec, float sigma_meas, const float * PH)
{
	float dot_prod =
			PH[0] * h_vec[3] +
			PH[1] * h_vec[4] +
			PH[2] * h_vec[5];
	dot_prod /= (sigma_meas*sigma_meas);

	return dot_prod;
} // END of function update_qhwh_trace

static float update_qhwh_trace_vp(const float * h_vec, float sigma_meas, const float * PH)
{
	float dot_prod =
			PH[3] * h_vec[0] +
			PH[4] * h_vec[1] +
			PH[5] * h_vec[2];
	dot_prod /= (sigma_meas*sigma_meas);

	return dot_prod;
} // END of function update_qhwh_trace

static float update_qhwh_trace_vv(const float * h_vec, float sigma_meas, const float * PH)
{
	float dot_prod =
			PH[3] * h_vec[3] +
			PH[4] * h_vec[4] +
			PH[5] * h_vec[5];
	dot_prod /= (sigma_meas*sigma_meas);

	return dot_prod;

} // END of function update_qhwh_trace

static float update_qhwhq_trace_pp(float sigma_meas, const float * PH)
{
	float dot_prod =
			PH[0] * PH[0] +
			PH[1] * PH[1] +
			PH[2] * PH[2];
	dot_prod /= (sigma_meas*sigma_meas);

	return dot_prod;

} // END of function update_qhwhq_trace

static float update_qhwhq_trace_pv(float sigma_meas, const float * PH)
{
	float dot_prod =
			PH[0] * PH[3] +
			PH[1] * PH[4] +
			PH[2] * PH[5];
	dot_prod /= (sigma_meas*sigma_meas);

	return dot_prod;

} // END of function update_qhwhq_trace

static float update_qhwhq_trace_vv(float sigma_meas, const float * PH)
{
	float dot_prod =
			PH[3] * PH[3] +
			PH[4] * PH[4] +
			PH[5] * PH[5];
	dot_prod /= (sigma_meas*sigma_meas);

	return dot_prod;

} // END of function update_qhwhq_trace

// This function calculates the Horizontal IBPL for a given confidence levels.
static char compute_single_ibpl(float conf_level, float num_dof, float indicator, float * HIBPL)
{
	*HIBPL = 1000.0;

	// Only calculate IBPL if there are enough measurements.
	if(num_dof > 0)
	{
		float aux = pow(conf_level, - 2 / num_dof) - 1;
		*HIBPL = sqrt(num_dof)*indicator*sqrt(aux);
		return 1;
	} else {
		return 0;
	}

	return 0;
}

static void compute_vec_2D(float Rmat[3][3], const float * hvec, float * hvec_2D)
{
	int k;
	float hR_p[3] = {0.0, 0.0, 0.0};
	float hR_v[3] = {0.0, 0.0, 0.0};
	for (k = 0; k < 3; k++)
	{
		hR_p[0] += hvec[k] * Rmat[0][k];
		hR_p[1] += hvec[k] * Rmat[1][k];
		hR_v[0] += hvec[k+3] * Rmat[0][k];
		hR_v[1] += hvec[k+3] * Rmat[1][k];
	}

	for (k = 0; k < 3; k++)
	{
		hvec_2D[k] = hR_p[0] * Rmat[0][k] + hR_p[1] * Rmat[1][k];
		hvec_2D[k+3] = hR_v[0] * Rmat[0][k] + hR_v[1] * Rmat[1][k];
		hvec_2D[k+6] = 0;
	}

}

static void update_dof(float max_dof, float r1, float r2, float ref_sigma, float acc_r, float first_dof, float second_dof, float * upd_dof)
{
	float min_n = MIN(first_dof, second_dof);

	if (min_n <= 0.5) {
		*upd_dof = 1.0;
		return;
	}

	float alpha = first_dof / 2;
	float beta = second_dof / 2;
	float A = r1 * r1;
	float B = r2 * r2;
	float C = ref_sigma * ref_sigma;

	float scale_factor = C + A + B;
	if (scale_factor < 1.0e-9) {
		*upd_dof = 1.0;
		return;
	}

	float r2_power = pow(B / scale_factor, beta);
	float argA = (alpha + beta) * beta * A / (alpha * scale_factor);
	float expA = exp(argA);

	float r1_power = pow(A / scale_factor, alpha);
	float argB = (alpha + beta) * alpha * B / (beta * scale_factor);
	float expB = exp(argB);

	double prob = r2_power * expA + r1_power * expB;

	double aux = 1 + C / (second_dof*acc_r*acc_r);
	*upd_dof = - 2 * log(prob) / log(aux);
	aux = 1 + C / (*upd_dof*acc_r*acc_r);
	*upd_dof = - 2 * log(prob) / log(aux);

	if (*upd_dof < min_n) {
		*upd_dof = min_n;
	}

	if (*upd_dof > max_dof) {
		*upd_dof = max_dof;
	}

	if (*upd_dof < 4) {
		*upd_dof = 4;
	}

}

static void fit_distribution(float dof_tail_depth, float max_dof, float R1[2][2], float dof_1, float R2[2][2], float dof_2, float Dm[2][2], float R[2][2], float * dof)
{
	int i, j;
	for (i = 0; i < 2; i++) {
		for (j = 0; j < 2; j++) {
			R[i][j] = R1[i][j] + R2[i][j] + Dm[i][j] + Dm[j][i];
		}
	}

	float t1 = R1[0][0] + R1[1][1];
	float t2 = R2[0][0] + R2[1][1];
	float t  = R[0][0]  + R[1][1];

	float A = sqrt(dof_1 * t1);
	float B = sqrt(dof_2 * t2);

	float ref_sigma = sqrt(A*A+B*B) * dof_tail_depth;
	update_dof(max_dof, A, B, ref_sigma, sqrt(t), dof_1, dof_2, dof);

}

static void add_ibpl_lines(
		char flag,
		char doppler_flag,
		double tow,
		float dof_tail_depth,
		float max_dof,
		float gamma,
		float U[2][2],
		const epoch_kibpl_parameters_t * epoch_parameters,
		const accumulated_kibpl_parameters_t * acc_parameters,
		accumulated_kibpl_parameters_t * out_parameters)
{
	// Compute correlation matrix Dm
	float Dm[2][2];
	float mu_A = sqrt(epoch_parameters->R[0][0]);
	float mu_B = sqrt(epoch_parameters->R[1][1]);
	float matrix[2][2];
	matrix[0][0] = epoch_parameters->r * mu_A * acc_parameters->u[0];
	matrix[0][1] = epoch_parameters->r * mu_A * acc_parameters->u[1];
	matrix[1][0] = epoch_parameters->r * mu_B * acc_parameters->u[0];
	matrix[1][1] = epoch_parameters->r * mu_B * acc_parameters->u[1];

	float Ut[2][2];
	Transpose(Ut, U, 2, 2);
	matMul_f(Dm, matrix, Ut, 2, 2, 2);


	float UR[2][2];
	float URUt[2][2];
	matMul_f(UR, U, acc_parameters->R, 2, 2, 2);
	matMul_f(URUt, UR, Ut, 2, 2, 2);

	float r2 = epoch_parameters->r * epoch_parameters->r;
	float R1[2][2];
	int i, j;
	for (i = 0; i < 2; i++) {
		for (j = 0; j < 2; j++) {
			R1[i][j] = r2 * epoch_parameters->R[i][j];
		}
	}

	fit_distribution(dof_tail_depth, max_dof, R1, epoch_parameters->eff_dof, URUt, acc_parameters->acc_dof, Dm, out_parameters->R, &out_parameters->acc_dof);

	// UPDATE CORRELATIONS
	out_parameters->u[0] = gamma * (epoch_parameters->r * sqrt(epoch_parameters->R[0][0]) +	acc_parameters->u[0] * Ut[0][0] + acc_parameters->u[1] * Ut[1][0]);
	out_parameters->u[1] = gamma * (epoch_parameters->r * sqrt(epoch_parameters->R[1][1]) +	acc_parameters->u[0] * Ut[0][1] + acc_parameters->u[1] * Ut[1][1]);

	if (flag)
	{
		FILE * fid = 0;
		if (doppler_flag == 0) {
			fid = fopen("/tmp/kibpl_dof_prange.txt", "a+");
		} else if (doppler_flag == 1) {
			fid = fopen("/tmp/kibpl_dof_doppler.txt", "a+");
		}

		if (fid != 0) {

			fprintf(fid, "%.4f  %.4f     %.4f  %.4f  %.4f  %.4f     %.4f  %.4f  %.4f  %.4f    %.4f  %.4f  %.4f  %.4f   %.4f  %.4f  %.4f   %.9f  %.9f  %.9f  %.9f\n",
					tow, epoch_parameters->r,
					epoch_parameters->R[0][0], epoch_parameters->R[0][1], epoch_parameters->R[1][0], epoch_parameters->R[1][1],
					U[0][0], U[0][1], U[1][0], U[1][1],
					acc_parameters->R[0][0], acc_parameters->R[0][1], acc_parameters->R[1][0], acc_parameters->R[1][1],
					epoch_parameters->eff_dof, acc_parameters->acc_dof, out_parameters->acc_dof,
					acc_parameters->u[0], acc_parameters->u[1], out_parameters->u[0], out_parameters->u[1]);
			fclose(fid);

		}
	}
}

static void compute_kalman_prange_indicators(double tow, float prange_dop_factor, float gamma_prange, const double * state, const double * delta_state, const float * kal_cov_mat, float Id_KH[2][2], kibpl_variables_t * kibpl_state, gnssdata_t *gnssdata)
{
	float prange_trace = 0.0;
	float prange_qhwh_trace_pp = 0.0;
	float prange_qhwh_trace_pv = 0.0;
	float prange_qhwh_trace_vp = 0.0;
	float prange_qhwh_trace_vv = 0.0;
	float prange_qhwhq_trace_pp = 0.0;
	float prange_qhwhq_trace_pv = 0.0;
	float prange_qhwhq_trace_vv = 0.0;

	int k;

	float Rmat[3][3];
	compute_nav_mat(state, Rmat);

	int num_valid_obs = 0;
	float residuals_sq = 0.0;
	for(k = 0; k < gnssdata->noOfChannelsAv; k++)
	{
		float h_vec[KNSTATES];
		float residual;
		compute_prange_residual(state, delta_state, &gnssdata->OBS[k], gnssdata->tow, tow, &residual, h_vec, KNSTATES);

		gnssdata->OBS[k].prange_residual = residual;
		if (gnssdata->OBS[k].prange_status == OBS_VALID) {
			residuals_sq += (residual * residual)  / (gnssdata->OBS[k].prange_sigma  * gnssdata->OBS[k].prange_sigma);

			// Calculate PH = P*h
			float PH [9];
			matSymVecMul_f(PH, kal_cov_mat, h_vec, 9);

			float hvec_2D[9];
			compute_vec_2D(Rmat, h_vec, hvec_2D);

			float PH_2D [9];
			compute_vec_2D(Rmat, PH, PH_2D);

			prange_trace += update_qhwh_trace(9, h_vec, gnssdata->OBS[k].prange_sigma, PH);;

			prange_qhwh_trace_pp += update_qhwh_trace_pp(hvec_2D, gnssdata->OBS[k].prange_sigma, PH_2D);
			prange_qhwh_trace_pv += update_qhwh_trace_pv(hvec_2D, gnssdata->OBS[k].prange_sigma, PH_2D);
			prange_qhwh_trace_vp += update_qhwh_trace_vp(hvec_2D, gnssdata->OBS[k].prange_sigma, PH_2D);
			prange_qhwh_trace_vv += update_qhwh_trace_vv(hvec_2D, gnssdata->OBS[k].prange_sigma, PH_2D);

			prange_qhwhq_trace_pp += update_qhwhq_trace_pp(gnssdata->OBS[k].prange_sigma, PH_2D);
			prange_qhwhq_trace_pv += update_qhwhq_trace_pv(gnssdata->OBS[k].prange_sigma, PH_2D);
			prange_qhwhq_trace_vv += update_qhwhq_trace_vv(gnssdata->OBS[k].prange_sigma, PH_2D);

			num_valid_obs++;
		}
	}

	if (prange_qhwh_trace_pp < 0) {
		kibpl_state->ibpl_epoch_flag = 0;
		return;
	}

	Id_KH[0][0] -= prange_qhwh_trace_pp / 2.0;
	Id_KH[0][1] -= prange_qhwh_trace_pv / 2.0;
	Id_KH[1][0] -= prange_qhwh_trace_vp / 2.0;
	Id_KH[1][1] -= prange_qhwh_trace_vv / 2.0;

	kibpl_state->prange_epoch.R[0][0] = prange_qhwhq_trace_pp / 2.0;
	kibpl_state->prange_epoch.R[0][1] = prange_qhwhq_trace_pv / 2.0;
	kibpl_state->prange_epoch.R[1][0] = prange_qhwhq_trace_pv / 2.0;
	kibpl_state->prange_epoch.R[1][1] = prange_qhwhq_trace_vv / 2.0;

	// Number of degrees of freedom and r
	float d_eff = prange_trace * (1 + gamma_prange);

	kibpl_state->prange_epoch.eff_dof = num_valid_obs - d_eff;
	if (kibpl_state->prange_epoch.eff_dof < 0) {
		kibpl_state->prange_epoch.eff_dof = 0;
	}

	kibpl_state->prange_epoch.r = 0.0;
	if (kibpl_state->prange_epoch.eff_dof > 0.0) {
		kibpl_state->prange_epoch.r = sqrt(residuals_sq / kibpl_state->prange_epoch.eff_dof);
	}

	kibpl_state->prange_epoch.r *= prange_dop_factor;

} // END of function update

static void compute_kalman_doppler_indicators(float doppler_dop_factor, float gamma_doppler, const double * state, const double * delta_state, const float * kal_cov_mat, float Id_KH[2][2], kibpl_variables_t * kibpl_state, gnssdata_t *gnssdata)
{
	float doppler_trace = 0.0;
	float doppler_qhwh_trace_pp = 0.0;
	float doppler_qhwh_trace_pv = 0.0;
	float doppler_qhwh_trace_vp = 0.0;
	float doppler_qhwh_trace_vv = 0.0;
	float doppler_qhwhq_trace_pp = 0.0;
	float doppler_qhwhq_trace_pv = 0.0;
	float doppler_qhwhq_trace_vv = 0.0;

	int k;

	float Rmat[3][3];
	compute_nav_mat(state, Rmat);

	int num_valid_obs = 0;
	float residuals_sq = 0.0;
	for(k = 0; k < gnssdata->noOfChannelsAv; k++)
	{
		float h_vec[KNSTATES];
		float residual;
		compute_doppler_residual(state, delta_state, &gnssdata->OBS[k], &residual, h_vec, KNSTATES);

		gnssdata->OBS[k].doppler_residual = residual;
		if (gnssdata->OBS[k].doppler_status == OBS_VALID) {
			residuals_sq += (residual * residual)  / (gnssdata->OBS[k].doppler_sigma  * gnssdata->OBS[k].doppler_sigma);

			// Calculate PH = P*h
			float PH [9];
			matSymVecMul_f(PH, kal_cov_mat, h_vec, 9);

			float hvec_2D[9];
			compute_vec_2D(Rmat, h_vec, hvec_2D);

			float PH_2D [9];
			compute_vec_2D(Rmat, PH, PH_2D);

			doppler_trace += update_qhwh_trace(9, h_vec, gnssdata->OBS[k].doppler_sigma, PH);

			doppler_qhwh_trace_pp += update_qhwh_trace_pp(hvec_2D, gnssdata->OBS[k].doppler_sigma, PH_2D);
			doppler_qhwh_trace_pv += update_qhwh_trace_pv(hvec_2D, gnssdata->OBS[k].doppler_sigma, PH_2D);
			doppler_qhwh_trace_vp += update_qhwh_trace_vp(hvec_2D, gnssdata->OBS[k].doppler_sigma, PH_2D);
			doppler_qhwh_trace_vv += update_qhwh_trace_vv(hvec_2D, gnssdata->OBS[k].doppler_sigma, PH_2D);;

			doppler_qhwhq_trace_pp += update_qhwhq_trace_pp(gnssdata->OBS[k].doppler_sigma, PH_2D);
			doppler_qhwhq_trace_pv += update_qhwhq_trace_pv(gnssdata->OBS[k].doppler_sigma, PH_2D);
			doppler_qhwhq_trace_vv += update_qhwhq_trace_vv(gnssdata->OBS[k].doppler_sigma, PH_2D);

			num_valid_obs++;
		}
	}

	if (num_valid_obs < 4) {
		kibpl_state->ibpl_epoch_flag = 0;
	}

	if (doppler_qhwh_trace_vv < 0) {
		kibpl_state->ibpl_epoch_flag = 0;
		return;
	}

	Id_KH[0][0] -= doppler_qhwh_trace_pp / 2.0;
	Id_KH[0][1] -= doppler_qhwh_trace_pv / 2.0;
	Id_KH[1][0] -= doppler_qhwh_trace_vp / 2.0;
	Id_KH[1][1] -= doppler_qhwh_trace_vv / 2.0;

	kibpl_state->doppler_epoch.R[0][0] = doppler_qhwhq_trace_pp / 2.0;
	kibpl_state->doppler_epoch.R[0][1] = doppler_qhwhq_trace_pv / 2.0;
	kibpl_state->doppler_epoch.R[1][0] = doppler_qhwhq_trace_pv / 2.0;
	kibpl_state->doppler_epoch.R[1][1] = doppler_qhwhq_trace_vv / 2.0;


	// Number of degrees of freedom and r
	float d_eff = doppler_trace * (1 + gamma_doppler);

	const float phi = 0.5;
	kibpl_state->doppler_epoch.eff_dof *= phi;

	float r2 = residuals_sq + kibpl_state->doppler_epoch.eff_dof * kibpl_state->doppler_epoch.r * kibpl_state->doppler_epoch.r;

	kibpl_state->doppler_epoch.eff_dof = (num_valid_obs - d_eff) + kibpl_state->doppler_epoch.eff_dof;
	if (kibpl_state->doppler_epoch.eff_dof < 0) {
		kibpl_state->doppler_epoch.eff_dof = 0;
	}

	kibpl_state->doppler_epoch.r = 0.0;
	if (kibpl_state->doppler_epoch.eff_dof > 0.0) {
		kibpl_state->doppler_epoch.r = sqrt(r2 / kibpl_state->doppler_epoch.eff_dof);
	}

	kibpl_state->doppler_epoch.r *= doppler_dop_factor;

	kibpl_state->doppler_epoch.dop_doppler = 1.0;
	float total_lambda_c = sqrt(kibpl_state->doppler_epoch.R[0][0]);
	if (doppler_qhwh_trace_pv > 0.0) {
		kibpl_state->doppler_epoch.dop_doppler = total_lambda_c / (doppler_qhwh_trace_pv / 2.0);
	}

} // END of function update

static void correct_kalman_doppler_indicators(float tow, const double * delta_state, float gamma_acc_indicator, float acc_indicator_factor, kibpl_variables_t * kibpl_state)
{
	// ACCELERATION INDICATOR
	float acc_indicator = 0.0;
	if (kibpl_state->acc_norm_imu > 0.0) {
		acc_indicator = kibpl_state->acc_norm_imu;
	} else {
		float delta_vel_norm_sq = delta_state[3] * delta_state[3] + delta_state[4] * delta_state[4] + delta_state[5] * delta_state[5];
		acc_indicator = sqrt(delta_vel_norm_sq);
	}

	kibpl_state->acc_indicator =
			sqrt(gamma_acc_indicator * gamma_acc_indicator * kibpl_state->acc_indicator * kibpl_state->acc_indicator + acc_indicator * acc_indicator);
	float reduced_indicator_sq = kibpl_state->acc_indicator * kibpl_state->acc_indicator;
	if (reduced_indicator_sq < 0.0) {
		reduced_indicator_sq = 0.0;
	}

	kibpl_state->doppler_epoch.r = sqrt(kibpl_state->doppler_epoch.r * kibpl_state->doppler_epoch.r + acc_indicator_factor * reduced_indicator_sq / (kibpl_state->doppler_epoch.dop_doppler * kibpl_state->doppler_epoch.dop_doppler));

} // END of function update

void initialize_kibpl_indicators(const kconf_t * conf, const lsq_state_t * lsq_state, double tow, kibpl_t * kibpl_state)
{
	kibpl_state->common.last_ibpl_tow = tow;
	kibpl_state->common.U[0][0] = 1;
	kibpl_state->common.U[0][1] = 0;
	kibpl_state->common.U[1][0] = 0;
	kibpl_state->common.U[1][1] = 1;

	// PSEUDORANGE
	kibpl_state->common.prange_epoch.r = lsq_state->prange_residuals;
	kibpl_state->common.prange_epoch.R[0][0] = lsq_state->hdop * lsq_state->hdop;
	kibpl_state->common.prange_epoch.R[0][1] = 0;
	kibpl_state->common.prange_epoch.R[1][0] = 0;
	kibpl_state->common.prange_epoch.R[1][1] = 0;
	kibpl_state->common.prange_epoch.eff_dof = lsq_state->num_pos_dof;

	kibpl_state->prange_parameters.acc_dof = kibpl_state->common.prange_epoch.eff_dof;
	float r2 = kibpl_state->common.prange_epoch.r * kibpl_state->common.prange_epoch.r;
	kibpl_state->prange_parameters.R[0][0] = r2 * kibpl_state->common.prange_epoch.R[0][0];
	kibpl_state->prange_parameters.R[0][1] = r2 * kibpl_state->common.prange_epoch.R[0][1];
	kibpl_state->prange_parameters.R[1][0] = r2 * kibpl_state->common.prange_epoch.R[1][0];
	kibpl_state->prange_parameters.R[1][1] = r2 * kibpl_state->common.prange_epoch.R[1][1];
	kibpl_state->prange_parameters.u[0] = conf->gamma_corr_prange * kibpl_state->common.prange_epoch.r * sqrt(kibpl_state->common.prange_epoch.R[0][0]);
	kibpl_state->prange_parameters.u[1] = conf->gamma_corr_prange * kibpl_state->common.prange_epoch.r * sqrt(kibpl_state->common.prange_epoch.R[1][1]);

	// DOPPLER
	kibpl_state->common.doppler_epoch.r = lsq_state->doppler_residuals;
	kibpl_state->common.doppler_epoch.R[0][0] = 0;
	kibpl_state->common.doppler_epoch.R[0][1] = 0;
	kibpl_state->common.doppler_epoch.R[1][0] = 0;
	kibpl_state->common.doppler_epoch.R[1][1] = lsq_state->hdop_vel * lsq_state->hdop_vel;
	kibpl_state->common.doppler_epoch.eff_dof = lsq_state->num_vel_dof;

	kibpl_state->doppler_parameters.acc_dof = kibpl_state->common.doppler_epoch.eff_dof;
	r2 = kibpl_state->common.doppler_epoch.r * kibpl_state->common.doppler_epoch.r;
	kibpl_state->doppler_parameters.R[0][0] = r2 * kibpl_state->common.doppler_epoch.R[0][0];
	kibpl_state->doppler_parameters.R[0][1] = r2 * kibpl_state->common.doppler_epoch.R[0][1];
	kibpl_state->doppler_parameters.R[1][0] = r2 * kibpl_state->common.doppler_epoch.R[1][0];
	kibpl_state->doppler_parameters.R[1][1] = r2 * kibpl_state->common.doppler_epoch.R[1][1];
	kibpl_state->doppler_parameters.u[0] = conf->gamma_corr_prange * kibpl_state->common.doppler_epoch.r * sqrt(kibpl_state->common.doppler_epoch.R[0][0]);
	kibpl_state->doppler_parameters.u[1] = conf->gamma_corr_prange * kibpl_state->common.doppler_epoch.r * sqrt(kibpl_state->common.doppler_epoch.R[1][1]);

	// ACCELERATION INDICATORS
	kibpl_state->common.acc_indicator = 0.0;
	kibpl_state->common.acc_norm_imu = 0.0;

}

void compute_kibpl_indicators(
		double tow,
		unsigned int state_size,
		const kconf_t * kconf,
		float acc_indicator_factor,
		const double * state,
		const double * delta_state,
		const float  * kal_cov_mat,
		kibpl_variables_t * kibpl_state,
		gnssdata_t * gnssdata)
{
	float prop_time = tow - kibpl_state->last_ibpl_tow;
	float F[2][2] = {{1,0},{0,1}};
	F[0][1] = prop_time;
	float Id_KH[2][2] = {{1,0},{0,1}};

	compute_kalman_prange_indicators(tow, kconf->prange_dop_factor, kconf->gamma_corr_prange, state, delta_state, kal_cov_mat, Id_KH, kibpl_state, gnssdata);
	compute_kalman_doppler_indicators(kconf->doppler_dop_factor, kconf->gamma_corr_doppler, state, delta_state, kal_cov_mat, Id_KH, kibpl_state, gnssdata);
	correct_kalman_doppler_indicators(tow, delta_state, kconf->gamma_acc_indicator, acc_indicator_factor, kibpl_state);

	matMul_f(kibpl_state->U, Id_KH, F, 2, 2, 2);

}

static char update_kibpl(char flag, const kconf_t * conf, double tow, kibpl_t * kibpl_state, float * hibpl)
{
	// UPDATE PSEUDORANGE INDICATORS
	add_ibpl_lines(
			flag, 0, tow,
			conf->dof_tail_depth, conf->max_dof, conf->gamma_corr_prange,
			kibpl_state->common.U,
			&kibpl_state->common.prange_epoch,
			&kibpl_state->prange_parameters,
			&kibpl_state->prange_parameters);

	// UPDATE DOPPLER INDICATORS
	add_ibpl_lines(
			flag, 1, tow,
			conf->dof_tail_depth, conf->max_dof, conf->gamma_corr_doppler,
			kibpl_state->common.U,
			&kibpl_state->common.doppler_epoch,
			&kibpl_state->doppler_parameters,
			&kibpl_state->doppler_parameters);

	// COMPUTE GLOBAL INDICATOR
	float Dm[2][2] = {{0, 0}, {0, 0}};
	float globalR[2][2] = {{0, 0}, {0, 0}};
	float global_dof = 0;
	if (conf->errors_absolute_sum == 0) {
		fit_distribution(
				conf->dof_tail_depth, conf->max_dof,
				kibpl_state->prange_parameters.R, kibpl_state->prange_parameters.acc_dof,
				kibpl_state->doppler_parameters.R, kibpl_state->doppler_parameters.acc_dof,
				Dm,
				globalR, &global_dof);
	}

	// COMPUTE IBPL
	int k;
	for (k = 0; k < NCONF_LVL; k++)
	{
		if (conf->errors_absolute_sum) {
			float hibpl_prange, hibpl_doppler;
			float indicator_prange = sqrt(kibpl_state->prange_parameters.R[0][0]);
			float indicator_doppler = sqrt(kibpl_state->doppler_parameters.R[0][0]);
			compute_single_ibpl(k < 7 ? pow(10, -(k+1)) : 0.05, kibpl_state->prange_parameters.acc_dof, indicator_prange, &hibpl_prange);
			compute_single_ibpl(k < 7 ? pow(10, -(k+1)) : 0.05, kibpl_state->doppler_parameters.acc_dof, indicator_doppler, &hibpl_doppler);

			hibpl[k] = hibpl_prange + hibpl_doppler;
		} else {
			compute_single_ibpl(k < 7 ? pow(10, -(k+1)) : 0.05, global_dof, sqrt(globalR[0][0]), &hibpl[k]);
		}

		if(isnan(hibpl[k]))
		{
			hibpl[k] = 0.0;
			return 0;
		}
	}

	return 1;
}

static kibpl_code_t check_ibpl_inputs(int kal_used_sats, const kibpl_t * kibpl_state)
{
	if (kal_used_sats < 4) {
		return KIBPL_NOT_VALID;
	}

	if ((kibpl_state->common.prange_epoch.eff_dof < 1) || (kibpl_state->common.doppler_epoch.eff_dof < 1)) {
		return KIBPL_NOT_VALID;
	}

	if (kibpl_state->common.ibpl_epoch_flag == 0) {
		return KIBPL_NOT_VALID;
	}

	return KIBPL_OK;
}

char GICSRx_LSQ_IBPL(lsq_state_t *lsq_state)
{

	int k;
	for (k = 0; k < NCONF_LVL; k++)
	{
		//compute_single_ibpl(pow(10, -(k+1)), lsq_state->num_pos_dof, lsq_state->prange_residuals * lsq_state->hdop, &lsq_state->hibpl[k]);
		lsq_state->hibpl[k] = PLtable_2D[lsq_state->num_pos_dof-1][k] * lsq_state->prange_residuals * lsq_state->hdop * sqrt(lsq_state->num_pos_dof);
		if(isnan(lsq_state->hibpl[k]))
		{
			return 0;
		}
	}

	return 1;
}

kibpl_code_t GICSRx_KF_IBPL(char flag, const kconf_t * conf, double tow, int kal_used_sats, kibpl_t * kibpl_state, nav_sol_t * nav_sol)
{
	// CHECK INPUTS
	kibpl_code_t kibpl_return_code = check_ibpl_inputs(kal_used_sats, kibpl_state);
	if (kibpl_return_code != KIBPL_OK) {
		return kibpl_return_code;
	}

	// COMPUTE KIBPL AND UPDATE PARAMETERS
	if (update_kibpl(flag, conf, tow, kibpl_state, nav_sol->hibpl) == 0) {
		return KIBPL_RESET_FILTER;
	}

	return KIBPL_OK;
}

#endif
