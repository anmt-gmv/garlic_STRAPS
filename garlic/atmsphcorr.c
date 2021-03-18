#include "atmsphcorr.h"

#include "algebra.h"
#include "calendar.h"
#include "NeQuick.h"

#include <math.h>
#include <limits.h>
#include <string.h>

#if CALC_ATMSPHCORR == 1

static float NeQuick_table_corr[MAX_N_SAT_G][3];

///////////////////////////////////////////////////////////////////////////////
//
// Function : Neill_tropo_model
//
// Purpose  : It computes tropospheric delay correction based on Neill model
//
// Args I   :
//            float  latitude    	user latitude         (radians)
//            float  height      	user height           (meters)
//			  float  elevation  	satellite elevation   (radians)
//            float  day         	day within year
// Args I/O :
// Args O   :
//            float  *delay_tropo  	tropospheric delay    (meters)
//			  float  *sigma_tropo   delay variance        (meters)
// Returns  : void
// Depends  :
// Calls    :
// Comments : Neill model for the computation of tropospheric delay.
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//          |         |          |      |                                    |
//          |   1.0   | 99/01/18 | GNSS |            First version           |
//			|   2.0   | 15/05/18 | CMVV |           Second version           |
//
///////////////////////////////////////////////////////////////////////////////
static void Neill_tropo_model(float latitude, float height, float elevation, float day, float *delay_tropo, float *sigma_tropo)
{
	const float lat[5]  = {15., 30., 45., 60., 75.};
	const float p0 [5]  = {1013.25, 1017.25, 1015.75, 1011.75, 1013.00};
	const float pd [5]  = {0., -3.75, -2.25, -1.75, -0.50};
	const float t0 [5]  = {299.65, 294.15, 283.15, 272.15, 263.65};
	const float td [5]  = {0.00, 7.00, 11.00, 15.00, 14.50};
	const float w0 [5]  = {26.31, 21.79, 11.66, 6.78, 4.11};
	const float wd [5]  = {0.00, 8.85, 7.24, 5.36, 3.39};
	const float e0 [5]  = {6.30e-3, 6.05e-3, 5.58e-3, 5.39e-3, 4.53e-3};
	const float ed [5]  = {0.00e-3, 0.25e-3, 0.32e-3, 0.81e-3, 0.62e-3};
	const float v0 [5]  = {2.77, 3.15, 2.57, 1.81, 1.55};
	const float vd [5]  = {0.00, 0.33, 0.46, 0.74, 0.30};

	float dmin, zhyd, zwet, dhyd, dwet;
	float press0, tempe0, watva0, telrt0, wvlrt0;
	float pressd, temped, watvad, telrtd, wvlrtd;
	float press, tempe, watva, telrt, wvlrt;
	float melev, dummy;
	int ind;

	double abslat = fabs(latitude*RADTODEG);

	if(abslat > lat[0] && abslat < lat[4])
	{
		if(abslat > lat[0] && abslat <= lat[1])
		{
			ind = 0;
		}
		else if(abslat <= lat[2])
		{
			ind = 1;
		}
		else if(abslat <= lat[3])
		{
			ind = 2;
		}
		else
		{
			ind = 3;
		}

		dummy = (abslat - lat[ind]) / (lat[ind + 1] - lat[ind]);
		press0 = p0[ind] + (p0[ind + 1] - p0[ind]) * dummy;
		tempe0 = t0[ind] + (t0[ind + 1] - t0[ind]) * dummy;
		watva0 = w0[ind] + (w0[ind + 1] - w0[ind]) * dummy;
		telrt0 = e0[ind] + (e0[ind + 1] - e0[ind]) * dummy;
		wvlrt0 = v0[ind] + (v0[ind + 1] - v0[ind]) * dummy;

		pressd = pd[ind] + (pd[ind + 1] - pd[ind]) * dummy;
		temped = td[ind] + (td[ind + 1] - td[ind]) * dummy;
		watvad = wd[ind] + (wd[ind + 1] - wd[ind]) * dummy;
		telrtd = ed[ind] + (ed[ind + 1] - ed[ind]) * dummy;
		wvlrtd = vd[ind] + (vd[ind + 1] - vd[ind]) * dummy;
	}
	else if(abslat <= lat[0])
	{
		press0 = p0[0];
		tempe0 = t0[0];
		watva0 = w0[0];
		telrt0 = e0[0];
		wvlrt0 = v0[0];

		pressd = pd[0];
		temped = td[0];
		watvad = wd[0];
		telrtd = ed[0];
		wvlrtd = vd[0];
	}
	else
	{
		press0 = p0[4];
		tempe0 = t0[4];
		watva0 = w0[4];
		telrt0 = e0[4];
		wvlrt0 = v0[4];

		pressd = pd[4];
		temped = td[4];
		watvad = wd[4];
		telrtd = ed[4];
		wvlrtd = vd[4];
	}

	// Compute the Dmin value depending on the latitude
	dmin = (latitude >= 0 ? 28. : 211.); // North : South.

	// Compute A-3 equation
	dummy = cos(2*PI_G*(day - dmin) / 365.25);
	press = press0 - pressd * dummy;
	tempe = tempe0 - temped * dummy;
	watva = watva0 - watvad * dummy;
	telrt = telrt0 - telrtd * dummy;
	wvlrt = wvlrt0 - wvlrtd * dummy;

	// Compute A-6 equation
	zhyd = (1.e-6 * 77.604 * 287.054 * press) / (9.784);

	// Compute A-7 equation
	zwet = (1.e-6 * 382000. * 287.054 * watva)/ ((9.784 * (1. + wvlrt) - 287.054 * telrt) * tempe);

	if (height < -100.) { height = -100.; }

	dummy = 1. - telrt * height / tempe;
	if(dummy <= 0.)
	{
		dhyd = 0.;
		dwet = 0.;
	}
	else
	{
		dhyd = zhyd * pow(dummy, 9.80665 / (287.054 * telrt));
		dwet = zwet	* pow(dummy, 9.80665 * (wvlrt + 1) / (287.054 * telrt) - 1.);
	}

	// CMVV: Do not apply minimum elevation correction, to match with RTKLIB.
	/*if (elevdeg < 5.)
	{
		elevdeg = 5.;
	}
	*/

	// Compute mapping function
	dummy = sin(elevation);
	melev = 1.001 / sqrt(0.002001 + dummy * dummy);

	// Compute A-11 equation
	*delay_tropo = (dhyd + dwet) * melev;
	*sigma_tropo = 0.12 * melev;
}

//////////////////////////////////////////////////////////////////////////////
//
// Function : get_sigma_uire
// Purpose  : Computes sigma_uire (Klobuchar model standard deviation)
//
// Args I   : double     delay_iono    Ionospheric delay (m)
//            double     PHIm          Geomagnetic latitude of the Earth
//                                     projection at IIP point (degrees)
//            double     elevation     Satellite elevation angle (radians)
// Args I/O :
// Args O   :
// Returns  : double                  Sigma_uire (m)
// Depends  :
// Calls    :
// Comments : Model from MOPS change 3 assumed.
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//          |         |          |      |                                    |
//          |   1.0   | 99/01/18 | GNSS |            First version           |
//          |   2.0   | 15/05/18 | CMVV |           Second version           |
//////////////////////////////////////////////////////////////////////////////
static float get_sigma_uire(float delay_iono, float PHIm, float elevation)
{
	const float I_H_IONO = 350000.0;

	float sigma_iono;

	float Fpp = 1/sqrt(1 - pow((R_EARTH*cos(elevation))/(R_EARTH + I_H_IONO), 2));

	float tv = 6.0;
	if(fabs(PHIm) <= 20)
	{
		tv = 9.0;
	}
	else if(fabs(PHIm) <= 55)
	{
		tv = 4.5;
	}

	sigma_iono = (delay_iono/5.0 > Fpp*tv ? delay_iono/5.0 : Fpp*tv);

	return sigma_iono;
}

////////////////////////////////////////////////////////////////////////////////
//
// Function : Klobuchar_iono_model
// Purpose  : Compute ionospheric corrections using Klobuchar model if available
//
// Args I   : double    userlat        User latitude                  (radians)
//            double    userlon        User longitude                 (radians)
//            double    elevation      Satellite elevation angle      (radians)
//            double    azimuth    	   Satellite azimuth   angle      (radians)
//            double    GPStime    	   GPS time within the week       (seconds)
//
//            ionoutc_gps_t *iono_gps  Structure containing Klobuchar parameters
//
// Args I/O :
// Args O   :
//            double    *delay_iono    Ionospheric error               (meters)
//            double    *sigma_iono    Sigma UIRE                      (meters)
// Returns  : char     valid
// Depends  :
// Calls    : get_sigma_uire
// Comments :
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE                 |
//          |         |          |      |                                      |
//          |   1.0   | 99/01/18 | GNSS |            First version             |
//          |   1.1   | 06/08/02 | EJGG |          Sign in x*4 term            |
//          |   2.0   | 15/05/18 | CMVV |             New release              |
////////////////////////////////////////////////////////////////////////////////
static char Klobuchar_iono_model(float userlat, float userlon, float elevation, float azimuth, float GPStime,
				ionoutc_gps_t *iono_gps, double gamma, float *delay_iono, float *sigma_iono)
{
	const float atmsph_epsilon = 1e-7;

	// Returns false if the structure for the ionospheric almanac is not ready.
	if (iono_gps->vflg == 0 || elevation < 0 || elevation > PI_G/2 || fabs(azimuth) > 2*PI_G)
	{
		return 0;
	}

	// Get elevation in semi-circles.
	float elevation_sc = elevation/PI_G;

	// Sub-ionospheric latitude
	float chi  = 0.0137/(elevation_sc + 0.11) - 0.022;
	float PHIi = userlat/PI_G + chi*cos(azimuth);

	if (PHIi > 0.416)
	{
		PHIi = 0.416;
	}
	else if (PHIi < -0.416)
	{
		PHIi = -0.416;
	}

	float cosPHIi = cos(PI_G*PHIi);

	// Sub-ionospheric longitude
	float lambdai = userlon/PI_G;
	if(fabs(cosPHIi) >= atmsph_epsilon)
	{
		lambdai += chi*sin(azimuth)/cosPHIi;
	}

	// Geo-magnetic latitude of sub-ionospheric location
	float PHIm = PHIi + 0.064*cos(PI_G*(lambdai - 1.617));

	float t = fmod(GPStime, DAYSECS) + lambdai*DAYSECS/2;
	if(t >= DAYSECS)
	{
		t -= DAYSECS;
	}
	else if (t < 0)
	{
		t += DAYSECS;
	}

	// Compute F, PER, AMP, x and delay_iono

	float Ft = 1.0 + 16*(0.53 - elevation_sc)*(0.53 - elevation_sc)*(0.53 - elevation_sc);

	float PER = iono_gps->beta[0] + PHIm * (iono_gps->beta[1] + PHIm*(iono_gps->beta[2] + PHIm * iono_gps->beta[3]));
	if(PER < 72000.0)
	{
		PER = 72000.0;
	}

	float AMP = iono_gps->alpha[0] + PHIm*(iono_gps->alpha[1] + PHIm*(iono_gps->alpha[2] + PHIm * iono_gps->alpha[3]));
	if(AMP < 0)
	{
		AMP = 0.0;
	}

	float x = 2*PI_G*(t - 50400.0)/PER;

	float tiono;
	if(fabs(x) < 1.57)
	{
		tiono = Ft*(5e-9 + AMP*(1 - x*x*(0.5 - (x*x/24))));
	}
	else
	{
		tiono = Ft*5e-9;
	}

	*delay_iono = SPEED_OF_LIGHT*tiono*gamma;

	// Compute standard deviation associated to the computed delay, taking into
	// account that the geo-magnetic latitude is expressed in semi-circles.
	*sigma_iono = get_sigma_uire(*delay_iono, PHIm*180, elevation);

	return 1;
}

static IntegrateData_st stIntegrateData;

void initialize_NeQuick_iono_model(char *mapping_path, ionoutc_gal_t *iono_gal)
{
	char bError1 = NeqReadMODIPFiles(mapping_path, &stIntegrateData.pstNeQuick.pstModip);
	char bError2 = NeqReadCCIRFiles (mapping_path, &stIntegrateData.pstNeQuick.pstCCIR);

	if(bError1 == 1 || bError2 == 1)
	{
//		fprintf(stderr, "Error reading NeQuick files\n");
		iono_gal->vflg[0] = 0;
	}

	stIntegrateData.pstNeQuick.siMaxRecurse    = 5;
	stIntegrateData.pstNeQuick.pdKronrodTol[0] = 0.1;
	stIntegrateData.pstNeQuick.pdKronrodTol[1] = 1.0;

	int k;
	for(k = 0; k < MAX_N_SAT_G; k++)
	{
		NeQuick_table_corr[k][0] = -100;
	}
}

static char NeQuick_iono_model(int month, double UTCtime, double gpos_user[3], double gpos_ssat[3], float elevation,
						ionoutc_gal_t *iono_gal, double signal_frequency, float *delay_iono, float *sigma_iono)
{
	if(iono_gal->vflg[0] == 0)
	{
		return 0;
	}

	NeQuick_status valid;

	stIntegrateData.pstNeQuick.siMonth = month;
	stIntegrateData.pstNeQuick.dUT     = UTCtime;

	stIntegrateData.pstNeQuick.pdGssPosLLH[0]  = gpos_user[0]*RADTODEG;
	stIntegrateData.pstNeQuick.pdGssPosLLH[1]  = gpos_user[1]*RADTODEG;
	stIntegrateData.pstNeQuick.pdGssPosLLH[2]  = gpos_user[2]/1000;

	stIntegrateData.pstNeQuick.pdSatPosLLH[0]  = gpos_ssat[0]*RADTODEG;
	stIntegrateData.pstNeQuick.pdSatPosLLH[1]  = gpos_ssat[1]*RADTODEG;
	stIntegrateData.pstNeQuick.pdSatPosLLH[2]  = gpos_ssat[2]/1000;

	stIntegrateData.pstNeQuick.siNumCoeff = 3;
	stIntegrateData.pstNeQuick.pdCoeff[0] = iono_gal->ai[0];
	stIntegrateData.pstNeQuick.pdCoeff[1] = iono_gal->ai[1];
	stIntegrateData.pstNeQuick.pdCoeff[2] = iono_gal->ai[2];

	double pdSTEC;
	valid = NeQuick(&stIntegrateData, &pdSTEC);

	*delay_iono = 40.3*pdSTEC/signal_frequency/signal_frequency;

	// The computation of the variance for the ionospheric delay does not correspond to the NeQuick
	// algorithm. The method used is described in MOPS as an approach for the Klobuchar model.
	double magnetic_dip = tan(stIntegrateData.dModipRx*DEGTORAD)*sqrt(cos(gpos_user[0]));

	double PHIm = atan(tan(magnetic_dip)/2)*RADTODEG;

	*sigma_iono = get_sigma_uire(*delay_iono, PHIm, elevation);

	return (valid == E_OK);
}

char calculate_atmsphCorr(int month, double day, double tow, double geod_pos[3], ionoutc_t *iono_utc, obsdata_t *obsdata, ionoModel_t iono_model)
{
	char update_eph = 0;

	float delay_tropo = 0.00, sigma_tropo = 0.00;
	float delay_iono  = 0.00, sigma_iono  = 0.00;


	// Compute tropospheric correction.
	Neill_tropo_model(geod_pos[0], geod_pos[2], obsdata->elev*DEGTORAD, day, &delay_tropo, &sigma_tropo);

	// Calculate ionospheric delay, using Klobuchar or NeQuick model.
	double signal_frequency = SPEED_OF_LIGHT/obsdata->lambda;

	if((iono_model == ION_AUTO || iono_model == ION_NEQU || (iono_model == ION_FREE && obsdata->iono_free == 0)) && iono_utc->gal.vflg[0] == 1)
	{
		// UTC time expressed in hours (approximated value, it does not use Galileo clock parameters).
		double UTCtime = fmod((tow - iono_utc->gal.dtlsf), DAYSECS)/3600;

		double gpos_ssat[3];
		ECEFtoNAV_pos(obsdata->sat.pos, gpos_ssat);

		if(fabs(NeQuick_table_corr[obsdata->PRN-1][0] - tow) < 60)
		{
			delay_iono = NeQuick_table_corr[obsdata->PRN-1][1];
		}
		else
		{
			NeQuick_iono_model(month, UTCtime, geod_pos, gpos_ssat, obsdata->elev*DEGTORAD, &iono_utc->gal,
						   	   signal_frequency, &delay_iono, &sigma_iono);

			NeQuick_table_corr[obsdata->PRN-1][0] = tow;
			NeQuick_table_corr[obsdata->PRN-1][1] = delay_iono;

			update_eph = 1;
		}
		NeQuick_table_corr[obsdata->PRN-1][2] = tow;
	}
	else if ((iono_model == ION_AUTO || iono_model == ION_KLOB || (iono_model == ION_FREE && obsdata->iono_free == 0)) && iono_utc->gps.vflg == 1)
	{
		double gamma = (BANDFREQL1*BANDFREQL1)/(signal_frequency*signal_frequency);
		Klobuchar_iono_model(geod_pos[0], geod_pos[1], obsdata->elev*DEGTORAD, obsdata->azim*DEGTORAD,
							tow, &iono_utc->gps, gamma, &delay_iono, &sigma_iono);
	}

	obsdata->atmsphCorr = delay_iono + delay_tropo;

	return update_eph;
}

void update_NeQuick_table(gnssdata_t *gnssdata, ionoutc_gal_t *iono_utc_gal, double geod_pos[3], int month)
{
	float delay_iono = 0;
	float sigma_iono = 0;

	int k, max_dist, max_index;
	for(k = 0, max_dist = 0, max_index = -1; k < gnssdata->noOfChannelsAv; k++)
	{
		if (gnssdata->OBS[k].prange_status == OBS_VALID) {
			int time_dist1 = ((int)(gnssdata->tow - NeQuick_table_corr[gnssdata->OBS[k].PRN-1][0]) + SECONDSONWEEK) % SECONDSONWEEK;
			int time_dist2 = ((int)(gnssdata->tow - NeQuick_table_corr[gnssdata->OBS[k].PRN-1][2]) + SECONDSONWEEK) % SECONDSONWEEK;
			if(time_dist1 > max_dist && time_dist2 < 5)
			{
				max_dist  = time_dist1;
				max_index = k;
			}
		}
	}

	// Get pointer to the last updated ephemeris.
	obsdata_t *obsdata = &gnssdata->OBS[max_index];

	if(max_dist >= 30)
	{
		// UTC time expressed in hours (approximated value, it does not use Galileo clock parameters).
		double UTCtime = fmod((gnssdata->tow - iono_utc_gal->dtlsf), DAYSECS)/3600;

		double signal_frequency = SPEED_OF_LIGHT/obsdata->lambda;

		double gpos_ssat[3];
		ECEFtoNAV_pos(obsdata->sat.pos, gpos_ssat);

		NeQuick_iono_model(month, UTCtime, geod_pos, gpos_ssat, obsdata->elev*DEGTORAD, iono_utc_gal,
				signal_frequency, &delay_iono, &sigma_iono);

		NeQuick_table_corr[obsdata->PRN-1][0] = gnssdata->tow;
		NeQuick_table_corr[obsdata->PRN-1][1] = delay_iono;
		NeQuick_table_corr[obsdata->PRN-1][2] = gnssdata->tow;
	}
}

#endif
