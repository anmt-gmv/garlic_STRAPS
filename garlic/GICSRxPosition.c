///////////////////////////////////////////////////////////////////////////////
//
// File     : GICSRxPosition.c
// Purpose  : Useful subroutines for getting user position with GPS and GLONASS
//
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE                 |
//          |         |          |      |                                      |
//          |   1.0   | 12/11/29 | cmvv |            First version             |
//
///////////////////////////////////////////////////////////////////////////////

#include "GICSRxPosition.h"

#include <stdio.h>
#include <string.h>

#include "matrix.h"
#include "algebra.h"
#include "calendar.h"

#define SATEPSILON        1.e-7
#define RELATIV_CONST     -4.442807633e-10
#define RELATIV_CONST_GAL -4.442807309e-10

/* 1msec, used to calculate position derivative and compute satellite velocity. */
#define DELTA_T_DIFF   0.001

/* Constants used for SP3 */
#define DAYHOURS             (24L)
#define MAXSP3INTERVALS      (8*DAYHOURS+1)
#define MININTERP			 (8)
#define BAD_CLK_VALUE_SP3    (999999.) /* Bad or absent clock values are to be set at 999999.*/
#define BAD_POS_VALUE_SP3    (0.)      /* Bad or absent positional values are to be set at 0.000000 */
#define SP3_LINE_LENGTH_SI   (90)
#define MAX_SATS_INSP3_SI    (85)

static double compute_tow_difference(double tow1, double tow2)
{
	double dt = tow1 - tow2;

	if(dt > HALFWEEK)
	{
		dt -= SECONDSONWEEK;
	}
	else if(dt < -HALFWEEK)
	{
		dt += SECONDSONWEEK;
	}
	return dt;
}

///////////////////////////////////////////////////////////////////////////////
//
// Function : GICSRxGPSPos
// Purpose  : Computes Position, Velocity in WGS84 and
//            eccentric anomaly for a GPS satellite.
//
// Args I   :  ephem     (ephgps_t*)   GPS Navigation message structure
//             TxTime    (double)      Desired time for output data
//	      	   RxTime    (double)      User reception time
//			   coef_TGD  (double)      Flag for TGD application
// Args I/O : pephsat    (ephemsat_t*) Satellite state vector in WGS84
// Args O   :
// Returns  : char                   1 if everything is OK
// Depends  :
// Calls    :
// Comments : GPS ICD-200 models for satellite position computation
//            using the broadcast navigation message.
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//          |         |          |      |                                    |
//          |   1.0   | 97/04/29 | GNSS |            First version           |
//	   		|   1.1   | 08/01/21 | SSKK | Position calculated to user        |
//	    	|         |          |      | reception time 		     		 |
//
///////////////////////////////////////////////////////////////////////////////
char GICSRxGPSPos(ephgps_t *ephem, double RxTime, double TxTime, ephemsat_t *pephsat, double coef_TGD)
{
	double SemiAxis, MeanMotion, MeanAnom, Old_EccAnom, trueAnom, EccAnom;
	double ArgLat, TwoArgLat, Radius;
	double Delta_ik, Delta_rk, Delta_Uk, Xk, Yk;
	double AscNode, e1, e2;
	double Delta_T, ik, a, b, GoodAnom;
	double cosT,sinT,cosAr,sinAr,cosAs,sinAs,cosi,sini;
	int count;

	char converged = 0;

	double aux[2];
	double travelTime;
	double Rsv0[3], Rsv1[3]; // for velocity computation

	const int    MAX_KEPLER_ITER = 100;
    const double KEPLER_EPSILON  = 1.e-12;

	// Mean Motion:
	// ------------
	SemiAxis  = ephem->dAroot;
	SemiAxis *= SemiAxis;

	if(SemiAxis < 10000.)
	{
		return 0;
	}

	MeanMotion = ephem->dDn + sqrt(MU_EARTH/(SemiAxis*SemiAxis*SemiAxis));

	// Compute lapsed time since toe when required for end of week cross-over
	Delta_T = compute_tow_difference(TxTime, ephem->dtoe);

	// Compute travel time.
	travelTime = compute_tow_difference(RxTime, TxTime);

	// Mean Anomaly:
	// -------------
	MeanAnom = ephem->dMo + MeanMotion*Delta_T;

	// Solves Kepler for eccentric anomaly:
	// ------------------------------------
	e1        = ephem->dEcc;
	e2        = e1*e1;
	converged = 0;
	count     = 0;
	EccAnom   = MeanAnom;

	while(converged == 0 && count++ < MAX_KEPLER_ITER)
	{
		Old_EccAnom = EccAnom;
		EccAnom     = MeanAnom + e1*sin(EccAnom);

		if(fabs(EccAnom - Old_EccAnom) < KEPLER_EPSILON)
		{
			converged = 1;
		}
	}

	// Obtains true anomaly:
	a = sqrt(1-e2)* sin(EccAnom);
	b = cos(EccAnom)-e1;
	if(fabs(a) < SATEPSILON && fabs(b) < SATEPSILON)
	{
		trueAnom = 0.;
	} else {
		trueAnom = atan2(a,b);
	}
	GoodAnom = EccAnom;
	EccAnom  = acos((e1+cos(trueAnom))/(1+e1*cos(trueAnom)));

	// Argument of Latitude:
	// ---------------------
	ArgLat = trueAnom + ephem->d_w;

	// Corrections:
	// ------------
	TwoArgLat = 2.*ArgLat;
	cosT      = cos(TwoArgLat);
	sinT      = sin(TwoArgLat);
	Delta_Uk  = ephem->dCuc*cosT + ephem->dCus*sinT;
	Delta_rk  = ephem->dCrc*cosT + ephem->dCrs*sinT;
	Delta_ik  = ephem->dCic*cosT + ephem->dCis*sinT;

	// Corrected radius, argument of latitude,
	// inclination and ascending node:
	// ---------------------------------------
	Radius = SemiAxis*(1 -e1*cos(EccAnom)) + Delta_rk;

	ArgLat = ArgLat + Delta_Uk;

	ik = ephem->d_io + Delta_ik + ephem->dIDOT*Delta_T;

	AscNode = ephem->dOMEGAo + Delta_T*(ephem->dOMEGADOT - OMEGA_EARTH) - OMEGA_EARTH*ephem->dtoe;

	// Position in Orbital Plane:
	// --------------------------
	cosAr = cos(ArgLat);
	sinAr = sin(ArgLat);
	Xk 	  = Radius*cosAr;
	Yk    = Radius*sinAr;

	// WGS-84 Position coordinates:
	// ----------------------------
	cosAs  = cos(AscNode);
	sinAs  = sin(AscNode);
	cosi   = cos(ik);
	sini   = sin(ik);

	// Satellite position at Tx time expressed in WGS-84 axes as they are at Tx time
	Rsv0[0] = Xk*cosAs - Yk*sinAs*cosi;
	Rsv0[1] = Xk*sinAs + Yk*cosAs*cosi;
	Rsv0[2] = Yk*sini;

	// Satellite position at Tx time expressed in WGS-84 axes as they are at Rx time
	pephsat->pos[0] = Rsv0[0] + OMEGA_EARTH*travelTime*Rsv0[1];
	pephsat->pos[1] = Rsv0[1] - OMEGA_EARTH*travelTime*Rsv0[0];
	pephsat->pos[2] = Rsv0[2];

	// Now the satellite velocity is computed using numerical differentiation between
	// two close satellite positions.

	/////////////////////////////////////////////////////////////////////////////////
	// This method to calculate satellite velocity is faster in computational terms
	// but is more unprecise. If satellite propagation is activated, this method is
	// not recommended.
	/////////////////////////////////////////////////////////////////////////////////
/*
 	//double Vxk, Vyk, derpar, rdot, rLatdot;

	// Calculate satellite velocity using keplerian orbit parameters.
	if(fabs(1.0-e2) < SATEPSILON) {
		derpar = 0.0;
	} else {
		derpar = sqrt(MU_EARTH/SemiAxis/(1.-e2));
	}
	rdot 	= derpar*e1*sin(trueAnom);
	rLatdot = derpar*(1 + e1*cos(trueAnom));
	Vxk 	= rdot*cosAr - rLatdot*sinAr;
	Vyk 	= rdot*sinAr + rLatdot*cosAr;

	Vsv[0] = Vxk*cosAs - Vyk*sinAs*cosi;
	Vsv[1] = Vxk*sinAs + Vyk*cosAs*cosi;
	Vsv[2] = Vyk*sini;

	// CMVV rotate sat velocity due to user reception time
	aux[0] = Vsv[0] + OMEGA_EARTH*TransitTime*Vsv[1];
	aux[1] = Vsv[1] - OMEGA_EARTH*TransitTime*Vsv[0];

	Vsv[0] = aux[0];
	Vsv[1] = aux[1];
*/
	/////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////

	// The velocity calculated with keplerian orbital parameters is ignored. It is
	// considered to be much more precise the computation of the satellite velocity
	// using numerical derivative. To do so, we calculate the satellite position
	// At = 1msec after.

	TxTime += DELTA_T_DIFF;

	// Compute lapsed time since toe when required for end of week cross-over
	Delta_T = compute_tow_difference(TxTime, ephem->dtoe);

	// Mean Anomaly:
	// -------------
	MeanAnom = ephem->dMo + MeanMotion*Delta_T;

	// Solves Kepler for eccentric anomaly:
	// ------------------------------------
	e1 = ephem->dEcc;
	e2 = e1*e1;
	EccAnom = MeanAnom;
	converged = 0;
	count     = 0;

	while(converged == 0 && count < MAX_KEPLER_ITER)
	{
		Old_EccAnom = EccAnom;
		EccAnom = MeanAnom + e1*sin(EccAnom);
		count++;
		if(fabs(EccAnom - Old_EccAnom) < KEPLER_EPSILON) { converged = 1; }
	}

	// Obtains true anomaly:
	a = sqrt(1-e2)* sin(EccAnom);
	b = cos(EccAnom)-e1;
	if(fabs(a)<SATEPSILON && fabs(b)<SATEPSILON) {
		trueAnom = 0.;
	} else {
		trueAnom = atan2(a,b);
	}
	GoodAnom = EccAnom;
	EccAnom = acos((e1+cos(trueAnom))/(1+e1*cos(trueAnom)));

	// Argument of Latitude:
	// ---------------------
	ArgLat = trueAnom + ephem->d_w;

	// Corrections:
	// ------------
	TwoArgLat = 2.*ArgLat;
	cosT = cos(TwoArgLat);
	sinT = sin(TwoArgLat);
	Delta_Uk = ephem->dCuc*cosT + ephem->dCus*sinT;
	Delta_rk = ephem->dCrc*cosT + ephem->dCrs*sinT;
	Delta_ik = ephem->dCic*cosT + ephem->dCis*sinT;

	// Corrected radius, argument of latitude,
	// inclination and ascending node:
	// ---------------------------------------
	Radius = SemiAxis*(1 -e1*cos(EccAnom)) + Delta_rk;

	ArgLat = ArgLat + Delta_Uk;

	ik = ephem->d_io + Delta_ik + ephem->dIDOT*Delta_T;

	AscNode = ephem->dOMEGAo + Delta_T*(ephem->dOMEGADOT - OMEGA_EARTH) - OMEGA_EARTH*ephem->dtoe;

	// Position in Orbital Plane:
	// --------------------------
	cosAr = cos(ArgLat);
	sinAr = sin(ArgLat);
	Xk 	  = Radius*cosAr;
	Yk    = Radius*sinAr;

	// WGS-84 Position coordinates:
	// ----------------------------
	cosAs  = cos(AscNode);
	sinAs  = sin(AscNode);
	cosi   = cos(ik);
	sini   = sin(ik);

	Rsv1[0] = Xk*cosAs - Yk*sinAs*cosi;
	Rsv1[1] = Xk*sinAs + Yk*cosAs*cosi;
	Rsv1[2] = Yk*sini;

	// Satellite inertial velocity at Tx time expressed
	// in WGS-84 axes as they are at Tx time.
	pephsat->vel[0] = (Rsv1[0] - Rsv0[0])/DELTA_T_DIFF - OMEGA_EARTH*Rsv0[1];
	pephsat->vel[1] = (Rsv1[1] - Rsv0[1])/DELTA_T_DIFF + OMEGA_EARTH*Rsv0[0];
	pephsat->vel[2] = (Rsv1[2] - Rsv0[2])/DELTA_T_DIFF;

	// Satellite inertial velocity at Tx time expressed
    // in WGS-84 axes as they are at Rx time.
	aux[0] = pephsat->vel[0] + OMEGA_EARTH*travelTime*pephsat->vel[1];
	aux[1] = pephsat->vel[1] - OMEGA_EARTH*travelTime*pephsat->vel[0];

	pephsat->vel[0] = aux[0];
	pephsat->vel[1] = aux[1];

	/////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////

	// Recompute delta_t with respect to time of clock.
	Delta_T = compute_tow_difference(TxTime, ephem->dTOC);

	// Calculate bias and drift values for the GPS satellite's clock.
	//
	// Satellite Clock Bias.
	pephsat->clk_bias = SPEED_OF_LIGHT*(ephem->dAfo + ephem->dAf1*Delta_T +
						ephem->dAf2*Delta_T*Delta_T - coef_TGD*ephem->dTGD);

	pephsat->clk_bias += RELATIV_CONST*ephem->dEcc*ephem->dAroot*sin(GoodAnom)*SPEED_OF_LIGHT;

	// Satellite Clock Drift.
	pephsat->clk_drift = SPEED_OF_LIGHT*(ephem->dAf1 + 2*ephem->dAf2*Delta_T);

	// Store Clock Drift Rate.
	pephsat->clk_driftrate = ephem->dAf2;

	return 1;

} // END of function GICSRxGPSPos

///////////////////////////////////////////////////////////////////////////////
//
// Function : GICSRxGalPos
// Purpose  : Computes Position, Velocity in WGS84 and
//            eccentric anomaly for a Galileo satellite.
//
// Args I   :  ephem     (ephgal_t*)   Galileo Navigation message structure
//             TxTime    (double)      Desired time for output data
//	      	   RxTime    (double)      User reception time
//			   bgdCorr   (double)      BGD correction to be directly applied
// Args I/O : pephsat    (ephemsat_t*) Satellite state vector in WGS84
// Args O   :
// Returns  : char                   1 if everything is OK
// Depends  :
// Calls    :
// Comments : GAL ICD models for satellite position computation
//            using the broadcast navigation message.
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//          |         |          |      |                                    |
//          |   1.0   | 15/04/14 | CMVV |            First version           |
//
///////////////////////////////////////////////////////////////////////////////
char GICSRxGalPos(ephgal_t *ephem, double RxTime, double TxTime, ephemsat_t *pephsat, double bgdCorr)
{
	double SemiAxis, MeanMotion, MeanAnom, Old_EccAnom, trueAnom, EccAnom;
	double ArgLat, TwoArgLat, Radius;
	double Delta_ik, Delta_rk, Delta_Uk, Xk, Yk;
	double AscNode, e1, e2;
	double Delta_T, ik, a, b, GoodAnom;
	double cosT,sinT,cosAr,sinAr,cosAs,sinAs,cosi,sini;
	int count;

	char converged = 0;

	double aux[2];
	double travelTime;
	double Rsv0[3], Rsv1[3]; // for velocity computation

	const int    MAX_KEPLER_ITER = 100;
    const double KEPLER_EPSILON  = 1.e-12;

	// Mean Motion:
	// ------------
	SemiAxis  = ephem->dAroot;
	SemiAxis *= SemiAxis;

	if(SemiAxis < 10000.)
	{
		return 0;
	}

	MeanMotion = ephem->dDn + sqrt(MU_EARTH_GAL/(SemiAxis*SemiAxis*SemiAxis));

	// Compute lapsed time since toe when required for end of week cross-over
	Delta_T = compute_tow_difference(TxTime, ephem->dtoe);

	// Compute travel time.
	travelTime = compute_tow_difference(RxTime, TxTime);

	// Mean Anomaly:
	// -------------
	MeanAnom = ephem->dMo + MeanMotion*Delta_T;

	// Solves Kepler for eccentric anomaly:
	// ------------------------------------
	e1        = ephem->dEcc;
	e2        = e1*e1;
	converged = 0;
	count     = 0;
	EccAnom   = MeanAnom;

	while(converged == 0 && count++ < MAX_KEPLER_ITER)
	{
		Old_EccAnom = EccAnom;
		EccAnom     = MeanAnom + e1*sin(EccAnom);

		if(fabs(EccAnom - Old_EccAnom) < KEPLER_EPSILON)
		{
			converged = 1;
		}
	}

	// Obtains true anomaly:
	a = sqrt(1-e2)* sin(EccAnom);
	b = cos(EccAnom)-e1;
	if(fabs(a) < SATEPSILON && fabs(b) < SATEPSILON)
	{
		trueAnom = 0.;
	} else {
		trueAnom = atan2(a,b);
	}
	GoodAnom = EccAnom;
	EccAnom  = acos((e1+cos(trueAnom))/(1+e1*cos(trueAnom)));

	// Argument of Latitude:
	// ---------------------
	ArgLat = trueAnom + ephem->d_w;

	// Corrections:
	// ------------
	TwoArgLat = 2.*ArgLat;
	cosT      = cos(TwoArgLat);
	sinT      = sin(TwoArgLat);
	Delta_Uk  = ephem->dCuc*cosT + ephem->dCus*sinT;
	Delta_rk  = ephem->dCrc*cosT + ephem->dCrs*sinT;
	Delta_ik  = ephem->dCic*cosT + ephem->dCis*sinT;

	// Corrected radius, argument of latitude,
	// inclination and ascending node:
	// ---------------------------------------
	Radius = SemiAxis*(1 - e1*cos(EccAnom)) + Delta_rk;

	ArgLat = ArgLat + Delta_Uk;

	ik = ephem->d_io + Delta_ik + ephem->dIDOT*Delta_T;

	AscNode = ephem->dOMEGAo + Delta_T*(ephem->dOMEGADOT - OMEGA_EARTH) - OMEGA_EARTH*ephem->dtoe;

	// Position in Orbital Plane:
	// --------------------------
	cosAr = cos(ArgLat);
	sinAr = sin(ArgLat);
	Xk 	  = Radius*cosAr;
	Yk    = Radius*sinAr;

	// WGS-84 Position coordinates:
	// ----------------------------
	cosAs  = cos(AscNode);
	sinAs  = sin(AscNode);
	cosi   = cos(ik);
	sini   = sin(ik);

	// Satellite position at Tx time expressed in WGS-84 axes as they are at Tx time
	Rsv0[0] = Xk*cosAs - Yk*sinAs*cosi;
	Rsv0[1] = Xk*sinAs + Yk*cosAs*cosi;
	Rsv0[2] = Yk*sini;

	// Satellite position at Tx time expressed in WGS-84 axes as they are at Rx time
	pephsat->pos[0] = Rsv0[0] + OMEGA_EARTH*travelTime*Rsv0[1];
	pephsat->pos[1] = Rsv0[1] - OMEGA_EARTH*travelTime*Rsv0[0];
	pephsat->pos[2] = Rsv0[2];

	// Now the satellite velocity is computed using numerical differentiation between
	// two close satellite positions.
	TxTime += DELTA_T_DIFF;

	// Compute lapsed time since toe when required for end of week cross-over
	Delta_T = compute_tow_difference(TxTime, ephem->dtoe);

	// Mean Anomaly:
	// -------------
	MeanAnom = ephem->dMo + MeanMotion*Delta_T;

	// Solves Kepler for eccentric anomaly:
	// ------------------------------------
	e1 = ephem->dEcc;
	e2 = e1*e1;
	EccAnom = MeanAnom;
	converged = 0;
	count     = 0;

	while(converged == 0 && count < MAX_KEPLER_ITER)
	{
		Old_EccAnom = EccAnom;
		EccAnom = MeanAnom + e1*sin(EccAnom);
		count++;
		if(fabs(EccAnom - Old_EccAnom) < KEPLER_EPSILON) { converged = 1; }
	}

	// Obtains true anomaly:
	a = sqrt(1-e2)* sin(EccAnom);
	b = cos(EccAnom)-e1;
	if(fabs(a) < SATEPSILON && fabs(b) < SATEPSILON) {
		trueAnom = 0.;
	} else {
		trueAnom = atan2(a,b);
	}
	GoodAnom = EccAnom;
	EccAnom = acos((e1+cos(trueAnom))/(1+e1*cos(trueAnom)));

	// Argument of Latitude:
	// ---------------------
	ArgLat = trueAnom + ephem->d_w;

	// Corrections:
	// ------------
	TwoArgLat = 2.*ArgLat;
	cosT = cos(TwoArgLat);
	sinT = sin(TwoArgLat);
	Delta_Uk = ephem->dCuc*cosT + ephem->dCus*sinT;
	Delta_rk = ephem->dCrc*cosT + ephem->dCrs*sinT;
	Delta_ik = ephem->dCic*cosT + ephem->dCis*sinT;

	// Corrected radius, argument of latitude,
	// inclination and ascending node:
	// ---------------------------------------
	Radius = SemiAxis*(1 -e1*cos(EccAnom)) + Delta_rk;

	ArgLat = ArgLat + Delta_Uk;

	ik = ephem->d_io + Delta_ik + ephem->dIDOT*Delta_T;

	AscNode = ephem->dOMEGAo + Delta_T*(ephem->dOMEGADOT - OMEGA_EARTH) - OMEGA_EARTH * ephem->dtoe;

	// Position in Orbital Plane:
	// --------------------------
	cosAr = cos(ArgLat);
	sinAr = sin(ArgLat);
	Xk 	  = Radius*cosAr;
	Yk    = Radius*sinAr;

	// WGS-84 Position coordinates:
	// ----------------------------
	cosAs  = cos(AscNode);
	sinAs  = sin(AscNode);
	cosi   = cos(ik);
	sini   = sin(ik);

	Rsv1[0] = Xk*cosAs - Yk*sinAs*cosi;
	Rsv1[1] = Xk*sinAs + Yk*cosAs*cosi;
	Rsv1[2] = Yk*sini;

	// Satellite inertial velocity at Tx time expressed
	// in WGS-84 axes as they are at Tx time.
	pephsat->vel[0] = (Rsv1[0] - Rsv0[0])/DELTA_T_DIFF - OMEGA_EARTH*Rsv0[1];
	pephsat->vel[1] = (Rsv1[1] - Rsv0[1])/DELTA_T_DIFF + OMEGA_EARTH*Rsv0[0];
	pephsat->vel[2] = (Rsv1[2] - Rsv0[2])/DELTA_T_DIFF;

	// Satellite inertial velocity at Tx time expressed
    // in WGS-84 axes as they are at Rx time.
	aux[0] = pephsat->vel[0] + OMEGA_EARTH*travelTime*pephsat->vel[1];
	aux[1] = pephsat->vel[1] - OMEGA_EARTH*travelTime*pephsat->vel[0];

	pephsat->vel[0] = aux[0];
	pephsat->vel[1] = aux[1];

	/////////////////////////////////////////////////////////////////////////////////

	// Recompute delta_t with respect to time of clock.
	Delta_T = compute_tow_difference(TxTime, ephem->dTOC);

	// Calculate bias and drift values for the Galileo satellite's clock.
	//
	// Satellite Clock Bias.
	pephsat->clk_bias = SPEED_OF_LIGHT*(ephem->dAfo + ephem->dAf1*Delta_T +
						ephem->dAf2*Delta_T*Delta_T - bgdCorr);

	pephsat->clk_bias += RELATIV_CONST_GAL*ephem->dEcc*ephem->dAroot*sin(GoodAnom)*SPEED_OF_LIGHT;

	// Satellite Clock Drift.
	pephsat->clk_drift = SPEED_OF_LIGHT*(ephem->dAf1 + 2*ephem->dAf2*Delta_T);

	// Store Clock Drift Rate.
	pephsat->clk_driftrate = ephem->dAf2;

	return 1;

} // END of function GICSRxGalPos

///////////////////////////////////////////////////////////////////////////////
//
// Function : GICSRxBEIPos
// Purpose  : Computes Position, Velocity in WGS84 and
//            eccentric anomaly for a BEI satellite.
//
// Args I   :  ephem     (ephgps_t*)   GPS Navigation message structure
//             TxTime    (double)      Desired time for output data
//	      	   RxTime    (double)      User reception time
//			   coef_TGD  (double)      Flag for TGD application
// Args I/O : pephsat    (ephemsat_t*) Satellite state vector in WGS84
// Args O   :
// Returns  : char                   1 if everything is OK
// Depends  :
// Calls    :
// Comments : BEI ICD models for satellite position computation
//            using the broadcast navigation message.
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//          |         |          |      |                                    |
//          |   1.0   | 97/04/29 | GNSS |            First version           |
//	   		|   1.1   | 08/01/21 | SSKK | Position calculated to user        |
//	    	|         |          |      | reception time 		     		 |
//
///////////////////////////////////////////////////////////////////////////////
char GICSRxBEIPos(ephbei_t *ephem, double RxTime, double TxTime, ephemsat_t *pephsat, double coef_TGD)
{
	double SemiAxis, MeanMotion, MeanAnom, Old_EccAnom, trueAnom, EccAnom;
	double ArgLat, TwoArgLat, Radius;
	double Delta_ik, Delta_rk, Delta_Uk, Xk, Yk;
	double AscNode, e1, e2;
	double Delta_T, ik, a, b, GoodAnom;
	double cosT,sinT,cosAr,sinAr,cosAs,sinAs,cosi,sini;
	int count;

	char converged = 0;

	double aux[2];
	double travelTime;
	double Rsv0[3], Rsv1[3]; // for velocity computation

	const int    MAX_KEPLER_ITER = 100;
    const double KEPLER_EPSILON  = 1.e-12;

	// Mean Motion:
	// ------------
	SemiAxis  = ephem->dAroot;
	SemiAxis *= SemiAxis;

	if(SemiAxis < 10000.)
	{
		return 0;
	}

	MeanMotion = ephem->dDn + sqrt(MU_EARTH_BEI/(SemiAxis*SemiAxis*SemiAxis));

	// Compute lapsed time since toe when required for end of week cross-over
	Delta_T = compute_tow_difference(TxTime, ephem->dtoe);

	// Compute travel time.
	travelTime = compute_tow_difference(RxTime, TxTime);

	// Mean Anomaly:
	// -------------
	MeanAnom = ephem->dMo + MeanMotion*Delta_T;

	// Solves Kepler for eccentric anomaly:
	// ------------------------------------
	e1        = ephem->dEcc;
	e2        = e1*e1;
	converged = 0;
	count     = 0;
	EccAnom   = MeanAnom;

	while(converged == 0 && count++ < MAX_KEPLER_ITER)
	{
		Old_EccAnom = EccAnom;
		EccAnom     = MeanAnom + e1*sin(EccAnom);

		if(fabs(EccAnom - Old_EccAnom) < KEPLER_EPSILON)
		{
			converged = 1;
		}
	}

	// Obtains true anomaly:
	a = sqrt(1-e2)* sin(EccAnom);
	b = cos(EccAnom)-e1;
	if(fabs(a) < SATEPSILON && fabs(b) < SATEPSILON)
	{
		trueAnom = 0.;
	} else {
		trueAnom = atan2(a,b);
	}
	GoodAnom = EccAnom;
	EccAnom  = acos((e1+cos(trueAnom))/(1+e1*cos(trueAnom)));

	// Argument of Latitude:
	// ---------------------
	ArgLat = trueAnom + ephem->d_w;

	// Corrections:
	// ------------
	TwoArgLat = 2.*ArgLat;
	cosT      = cos(TwoArgLat);
	sinT      = sin(TwoArgLat);
	Delta_Uk  = ephem->dCuc*cosT + ephem->dCus*sinT;
	Delta_rk  = ephem->dCrc*cosT + ephem->dCrs*sinT;
	Delta_ik  = ephem->dCic*cosT + ephem->dCis*sinT;

	// Corrected radius, argument of latitude,
	// inclination and ascending node:
	// ---------------------------------------
	Radius = SemiAxis*(1 -e1*cos(EccAnom)) + Delta_rk;

	ArgLat = ArgLat + Delta_Uk;

	ik = ephem->d_io + Delta_ik + ephem->dIDOT*Delta_T;

	AscNode = ephem->dOMEGAo + Delta_T*(ephem->dOMEGADOT - OMEGA_EARTH_BEI) - OMEGA_EARTH_BEI*ephem->dtoe_bei;

	// Position in Orbital Plane:
	// --------------------------
	cosAr = cos(ArgLat);
	sinAr = sin(ArgLat);
	Xk 	  = Radius*cosAr;
	Yk    = Radius*sinAr;

	// WGS-84 Position coordinates:
	// ----------------------------
	cosAs  = cos(AscNode);
	sinAs  = sin(AscNode);
	cosi   = cos(ik);
	sini   = sin(ik);

	// Satellite position at Tx time expressed in WGS-84 axes as they are at Tx time
	Rsv0[0] = Xk*cosAs - Yk*sinAs*cosi;
	Rsv0[1] = Xk*sinAs + Yk*cosAs*cosi;
	Rsv0[2] = Yk*sini;

	// Satellite position at Tx time expressed in WGS-84 axes as they are at Rx time
	pephsat->pos[0] = Rsv0[0] + OMEGA_EARTH_BEI*travelTime*Rsv0[1];
	pephsat->pos[1] = Rsv0[1] - OMEGA_EARTH_BEI*travelTime*Rsv0[0];
	pephsat->pos[2] = Rsv0[2];

	// Now the satellite velocity is computed using numerical differentiation between
	// two close satellite positions.

	/////////////////////////////////////////////////////////////////////////////////
	// This method to calculate satellite velocity is faster in computational terms
	// but is more unprecise. If satellite propagation is activated, this method is
	// not recommended.
	/////////////////////////////////////////////////////////////////////////////////
/*
 	//double Vxk, Vyk, derpar, rdot, rLatdot;

	// Calculate satellite velocity using keplerian orbit parameters.
	if(fabs(1.0-e2) < SATEPSILON) {
		derpar = 0.0;
	} else {
		derpar = sqrt(MU_EARTH/SemiAxis/(1.-e2));
	}
	rdot 	= derpar*e1*sin(trueAnom);
	rLatdot = derpar*(1 + e1*cos(trueAnom));
	Vxk 	= rdot*cosAr - rLatdot*sinAr;
	Vyk 	= rdot*sinAr + rLatdot*cosAr;

	Vsv[0] = Vxk*cosAs - Vyk*sinAs*cosi;
	Vsv[1] = Vxk*sinAs + Vyk*cosAs*cosi;
	Vsv[2] = Vyk*sini;

	// CMVV rotate sat velocity due to user reception time
	aux[0] = Vsv[0] + OMEGA_EARTH*TransitTime*Vsv[1];
	aux[1] = Vsv[1] - OMEGA_EARTH*TransitTime*Vsv[0];

	Vsv[0] = aux[0];
	Vsv[1] = aux[1];
*/
	/////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////

	// The velocity calculated with keplerian orbital parameters is ignored. It is
	// considered to be much more precise the computation of the satellite velocity
	// using numerical derivative. To do so, we calculate the satellite position
	// At = 1msec after.

	TxTime += DELTA_T_DIFF;

	// Compute lapsed time since toe when required for end of week cross-over
	Delta_T = compute_tow_difference(TxTime, ephem->dtoe);

	// Mean Anomaly:
	// -------------
	MeanAnom = ephem->dMo + MeanMotion*Delta_T;

	// Solves Kepler for eccentric anomaly:
	// ------------------------------------
	e1 = ephem->dEcc;
	e2 = e1*e1;
	EccAnom = MeanAnom;
	converged = 0;
	count     = 0;

	while(converged == 0 && count < MAX_KEPLER_ITER)
	{
		Old_EccAnom = EccAnom;
		EccAnom = MeanAnom + e1*sin(EccAnom);
		count++;
		if(fabs(EccAnom - Old_EccAnom) < KEPLER_EPSILON) { converged = 1; }
	}

	// Obtains true anomaly:
	a = sqrt(1-e2)* sin(EccAnom);
	b = cos(EccAnom)-e1;
	if(fabs(a)<SATEPSILON && fabs(b)<SATEPSILON) {
		trueAnom = 0.;
	} else {
		trueAnom = atan2(a,b);
	}
	GoodAnom = EccAnom;
	EccAnom = acos((e1+cos(trueAnom))/(1+e1*cos(trueAnom)));

	// Argument of Latitude:
	// ---------------------
	ArgLat = trueAnom + ephem->d_w;

	// Corrections:
	// ------------
	TwoArgLat = 2.*ArgLat;
	cosT = cos(TwoArgLat);
	sinT = sin(TwoArgLat);
	Delta_Uk = ephem->dCuc*cosT + ephem->dCus*sinT;
	Delta_rk = ephem->dCrc*cosT + ephem->dCrs*sinT;
	Delta_ik = ephem->dCic*cosT + ephem->dCis*sinT;

	// Corrected radius, argument of latitude,
	// inclination and ascending node:
	// ---------------------------------------
	Radius = SemiAxis*(1 -e1*cos(EccAnom)) + Delta_rk;

	ArgLat = ArgLat + Delta_Uk;

	ik = ephem->d_io + Delta_ik + ephem->dIDOT*Delta_T;

	AscNode = ephem->dOMEGAo + Delta_T*(ephem->dOMEGADOT - OMEGA_EARTH_BEI) - OMEGA_EARTH_BEI*ephem->dtoe_bei;

	// Position in Orbital Plane:
	// --------------------------
	cosAr = cos(ArgLat);
	sinAr = sin(ArgLat);
	Xk 	  = Radius*cosAr;
	Yk    = Radius*sinAr;

	// WGS-84 Position coordinates:
	// ----------------------------
	cosAs  = cos(AscNode);
	sinAs  = sin(AscNode);
	cosi   = cos(ik);
	sini   = sin(ik);

	Rsv1[0] = Xk*cosAs - Yk*sinAs*cosi;
	Rsv1[1] = Xk*sinAs + Yk*cosAs*cosi;
	Rsv1[2] = Yk*sini;

	// Satellite inertial velocity at Tx time expressed
	// in WGS-84 axes as they are at Tx time.
	pephsat->vel[0] = (Rsv1[0] - Rsv0[0])/DELTA_T_DIFF - OMEGA_EARTH_BEI*Rsv0[1];
	pephsat->vel[1] = (Rsv1[1] - Rsv0[1])/DELTA_T_DIFF + OMEGA_EARTH_BEI*Rsv0[0];
	pephsat->vel[2] = (Rsv1[2] - Rsv0[2])/DELTA_T_DIFF;

	// Satellite inertial velocity at Tx time expressed
    // in WGS-84 axes as they are at Rx time.
	aux[0] = pephsat->vel[0] + OMEGA_EARTH_BEI*travelTime*pephsat->vel[1];
	aux[1] = pephsat->vel[1] - OMEGA_EARTH_BEI*travelTime*pephsat->vel[0];

	pephsat->vel[0] = aux[0];
	pephsat->vel[1] = aux[1];

	/////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////

	// Recompute delta_t with respect to time of clock.
	Delta_T = compute_tow_difference(TxTime, ephem->dTOC);

	// Calculate bias and drift values for the GPS satellite's clock.
	//
	// Satellite Clock Bias.
	pephsat->clk_bias = SPEED_OF_LIGHT*(ephem->dAfo + ephem->dAf1*Delta_T +
						ephem->dAf2*Delta_T*Delta_T - coef_TGD*ephem->dTGD_B1B3);

	pephsat->clk_bias += RELATIV_CONST_GAL*ephem->dEcc*ephem->dAroot*sin(GoodAnom)*SPEED_OF_LIGHT;

	// Satellite Clock Drift.
	pephsat->clk_drift = SPEED_OF_LIGHT*(ephem->dAf1 + 2*ephem->dAf2*Delta_T);

	// Store Clock Drift Rate.
	pephsat->clk_driftrate = ephem->dAf2;

	return 1;

} // END of function GICSRxBEIPos

//////////////////////////////////////////////////////////////////////////////////
// Name         : r_part
// Purpose      : Builds vec_out = d/dt (vec_in)
// Argument I/O : NAME     TYPE          I/O   DESCRIPTION
//
//                vec_in   (double[6])    I    Vector that contains the
//                                             position and velocity of the
//                                             satellite in km and km/sec
//                acc      (double[3])    I    Perturbation acceleration in
//                                             km/sec2
//                vec_out  (double[6])    O    Vector that represents the
//                                             derivative of vec_in
//
// I/O files    : None
// Returns      : char (TRUE if everything OK, FALSE if position could
//                not be obtained)
// History      :
//////////////////////////////////////////////////////////////////////////////////
static char r_part(double vec_in[6], double vec_out[6], double acc[3])
{
	// Definition of variables.
	double r2, r3, r5, a, b, c;
	char   result = 1;

	// Square of modulus of position: X2+Y2+Z2
	r2 = vec_in[0]*vec_in[0] + vec_in[1]*vec_in[1] + vec_in[2]*vec_in[2];

	// No valid position.
	if(fabs(r2) < 1.) {	result = 0;	}

	// Third and fifth order powers of the modulus.
	r3 = r2 * sqrt(r2);
	r5 = r3 * r2;

	// Several coefficients required are computed separately.
	a = 5.0*vec_in[2]*vec_in[2]/r2;

	b = -MU_EARTH_GLO/r3 + 1.5*C20_GLO*MU_EARTH_GLO*AE_GLO*AE_GLO*(1-a)/r5 +
			OMEGA_EARTH_GLO*OMEGA_EARTH_GLO;

	c = -MU_EARTH_GLO/r3 + 1.5*C20_GLO*MU_EARTH_GLO*AE_GLO*AE_GLO*(3-a)/r5;

	vec_out[0] = vec_in[3];
	vec_out[1] = vec_in[4];
	vec_out[2] = vec_in[5];

	// Acceleration in earth-fixed axes (PZ-90) in accordance with ICD.
	vec_out[3] = vec_in[0]*b + 2.0*OMEGA_EARTH_GLO*vec_in[4] + acc[0];
	vec_out[4] = vec_in[1]*b - 2.0*OMEGA_EARTH_GLO*vec_in[3] + acc[1];
	vec_out[5] = vec_in[2]*c + acc[2];

	return result;
}

/////////////////////////////////////////////////////////////////////////////////
// Name         : step_rk4
// Purpose      : Performs a single step of the Runge-Kutta integration
//
// Argument I/O : NAME     TYPE          I/O   DESCRIPTION
//
//                vec_io   (double[6])   I/O   Vector that contains the
//                                             position and velocity of the
//                                             satellite
//                acc      (double[3])    I    Perturbation acceleration in
//                                             km/sec2
//                step     (double)       I    Integration step
//
// I/O files    : None
// Returns      : bool (TRUE if everything OK, FALSE if position could
//                not be obtained)
// History      :
//////////////////////////////////////////////////////////////////////////////////
static char step_rk4(double vec_io[6], double acc[3], double step)
{
	// Definition of variables.
	double vec_in[6], vec_out[6], sum[6];
	char   result = 1;
	int    i;

	for(i = 0; i < 6; i++) { vec_in[i] = vec_io[i]; }

	// The derivative of vec_in is obtained, that is, employing the
	// equations of the acceleration in earth-fixed coordinates.
	result = (r_part(vec_in, vec_out, acc) != 0);

	// The first term of the R-K integration is obtained. If the discrete
	// differential equation is expressed as follows:
	// V(n+1) = V(n) + step/6 * (k1+2*k2+2*k3+k4), and calling the acceleration
	// F(tn,V), the terms of the integration are:
	// k1 = F (tn, V(n))
	// k2 = F (tn+step/2, V(n)+step/2*k1)
	// k3 = F (tn+step/2, V(n)+step/2*k2)
	// k4 = F (tn+step, V(n)+step*k3)
	// Note that to obtain k2, k3 and k4, the only dependencies of the
	// acceleration with time are present in the position and velocity

	for(i = 0; i < 6; i++)
	{
		sum[i]    = vec_out[i];  					  // k1 = vec_out
		vec_in[i] = vec_io[i] + 0.5*step*vec_out[i];  // vec_in = V(n)+step/2*k1
	}

	result = (r_part(vec_in, vec_out, acc) != 0);  	  // k2 = vec_out

	for(i = 0; i < 6; i++)
	{
		sum[i]    += 2.0*vec_out[i];  				  // sum = k1 + 2*k2
		vec_in[i]  = vec_io[i] + 0.5*step*vec_out[i]; // vec_in = V(n)+step/2*k2
	}

	result = (r_part(vec_in, vec_out, acc) != 0);     // k3 = vec_out

	for(i = 0; i < 6; i++)
	{
		sum[i]    += 2.0*vec_out[i];  				  // sum = k1 + 2*k2 + 2*k3
		vec_in[i]  = vec_io[i] + step*vec_out[i];  	  // vec_in = V(n)+step*k3
	}

	result = (r_part(vec_in, vec_out, acc) != 0);  	  // k4 = vec_out

	for(i = 0; i < 6; i++)
	{
		vec_io[i] += (sum[i] + vec_out[i])*step/6.0;  // V(n+1) = vec_io
	}

	return result;
}

//////////////////////////////////////////////////////////////////////////////////
// Name         : rk4
// Purpose      : Performs the 4th order Runge-Kutta integration algorithm for
//                the state vector propagation
// Argument I/O : NAME     TYPE          I/O   DESCRIPTION
//
//                vec_io   (double[6])   I/O   Vector that contains the
//                                             position and velocity
//                acc      (double[3])    I    Perturbation acceleration in
//                                             km/sec2
//                t0       (double)       I    Origin of integration interval
//                tend     (double)       I    Integration time (desired time
//                                             for output data)
//                step     (double)       I    Integration step
//
// I/O files    : None
// Returns      : char (TRUE if everything OK, FALSE if position could
//                not be obtained)
// History      :
//////////////////////////////////////////////////////////////////////////////////
static char rk4(double vec_io[6], double acc[3], double t0, double tend, double step)
{
	double t, h;
	char result = 1;

	// If the condition tend < t0 is fulfilled, step is negative, backwards
	// integration: h starts being the step
	h = (tend >= t0 ? +step : -step);

	// The integration is performed step by step, calling function step_rk4
	// successively until the integration point is closer to the end than the
	// nominal step.
	for (t = t0; fabs(tend-t) > fabs(h); t = t+h)
	{
		result = (step_rk4(vec_io, acc, h) != 0);
	}

	// Last step in the integration: h is less than the R-K nominal step.
	h = tend - t;

	// Last step in the integration: h is less than the R-K nominal step.
	result = (step_rk4(vec_io, acc, h) != 0);  	  // k4 = vec_out

	return result;
}

//////////////////////////////////////////////////////////////////////////////////
// Name         : GICSRxGLOPos
// Purpose      : Computes position and velocity in PZ-90 reference frame for
//                a GLONASS Satellite in accordance with GLONASS ICD, using
//                broadcast ephemeris, and converts them to WGS-84.
//                The method employed consists of calculating the acceleration
//                of the SV in an earth-fixed frame and integrating it using a
//                fourth-degree Runge-Kutta technique.
//
// Argument I/O : NAME     TYPE          I/O
//   DESCRIPTION
//
//                ephem	   (TEPHGLO)  	  I    GLONASS Navigation message
//                                             structure
//                TxTime (double)         I    Desired time for output data (time
//                                             of day)
//                Rsv      (double[3])    O    Satellite position in WGS-84 (m)
//                Vsv      (double[3])    O    Satellite velocity in WGS-84 (m/s)
// I/O files    : None
// Returns      : char (1 if everything OK, 0 if position could not be obtained)
// History      :
//////////////////////////////////////////////////////////////////////////////////
char GICSRxGLOPos(ephglo_t *ephem, double RxTime, double TxTime, ephemsat_t *pephsat, double coef_TGD)
{
	// Declaration of variables.
	double SV_state[6], acc[3], aux[2], Delta_T;
	double Rsv_PZ90[3], Vsv_PZ90[3];
	char result = 1;

	// Transmission time of day.
	double TxTOD = fmod(TxTime, DAYSECS);

	// Store satellite status at tb:
	//
	// Position.
	SV_state[0] = ephem->X_tb;
	SV_state[1] = ephem->Y_tb;
	SV_state[2] = ephem->Z_tb;

	// Velocity.
	SV_state[3] = ephem->Xdot_tb;
	SV_state[4] = ephem->Ydot_tb;
	SV_state[5] = ephem->Zdot_tb;

	// Acceleration (do not vary during 15min, constant).
	acc[0]      = ephem->Xdot2_tb;
	acc[1]      = ephem->Ydot2_tb;
	acc[2]      = ephem->Zdot2_tb;

	// Calculate time gap between transmission time and GLONASS 'toe'. Both variables are in Time of Day.
	Delta_T = TxTOD - ephem->tb;

	// Correct day rollovers.
	if (Delta_T < -HALFDAY)
	{
		TxTOD   += DAYSECS;
		Delta_T += DAYSECS;
	}
	else if (Delta_T > HALFDAY)
	{
		TxTOD   -= DAYSECS;
		Delta_T -= DAYSECS;
	}

	// If the time gap is not included in the ephemeris interval, return.
	if (fabs(Delta_T) > (QUARTHOUR_GLO + MARGIN_TB_GLO))
	{
		result = 0;
	}

	// Apply 4th order Runge-Kutta to calculate the current position.
	if(rk4(SV_state, acc, ephem->tb, TxTOD, RK_STEP_GLO) == 0)
	{
		// Cannot be computed.
		result = 0;
	}
	else
	{
		// Store calculated position and velocity.
		//
		// Position in PZ-90, convert from km to meters.
		Rsv_PZ90[0] = SV_state[0] * 1000;
		Rsv_PZ90[1] = SV_state[1] * 1000;
		Rsv_PZ90[2] = SV_state[2] * 1000;

		// Velocity in PZ-90.
		Vsv_PZ90[0] = SV_state[3] * 1000;
		Vsv_PZ90[1] = SV_state[4] * 1000;
		Vsv_PZ90[2] = SV_state[5] * 1000;

		// Convert relative velocity to inertial in PZ-90.
		Vsv_PZ90[0] -= OMEGA_EARTH_GLO * Rsv_PZ90[1];
		Vsv_PZ90[1] += OMEGA_EARTH_GLO * Rsv_PZ90[0];

		// Use PZ90 position as WGS-84.
		pephsat->pos[0] = Rsv_PZ90[0]; pephsat->pos[1] = Rsv_PZ90[1]; pephsat->pos[2] = Rsv_PZ90[2];
		pephsat->vel[0] = Vsv_PZ90[0]; pephsat->vel[1] = Vsv_PZ90[1]; pephsat->vel[2] = Vsv_PZ90[2];

		aux[0] = pephsat->pos[0] + OMEGA_EARTH_GLO*(RxTime - TxTime)*pephsat->pos[1];
		aux[1] = pephsat->pos[1] - OMEGA_EARTH_GLO*(RxTime - TxTime)*pephsat->pos[0];

		pephsat->pos[0] = aux[0];
		pephsat->pos[1] = aux[1];

		aux[0] = pephsat->vel[0] + OMEGA_EARTH_GLO*(RxTime - TxTime)*pephsat->vel[1];
		aux[1] = pephsat->vel[1] - OMEGA_EARTH_GLO*(RxTime - TxTime)*pephsat->vel[0];

		pephsat->vel[0] = aux[0];
		pephsat->vel[1] = aux[1];

		// Calculate clock bias and drift in meters.
		pephsat->clk_bias = SPEED_OF_LIGHT*(-ephem->SV_timedev + ephem->SV_freqdev*Delta_T - coef_TGD*ephem->dTGD);

		pephsat->clk_drift = SPEED_OF_LIGHT*ephem->SV_freqdev;

		// Store clock drift rate.
		pephsat->clk_driftrate = 0;

		/* Ephemeris refresh: the inclusion of this refinement in the algorithm
		 * permits saving time, as the origin of the interval is moved to the
		 * point computed last. In principle, this moves the tb reference time
		 * so it should not be included unless process time restrictions, which
		 * is precisely the case when using TESEO II.
		 */
		ephem->SV_timedev -= ephem->SV_freqdev*(TxTOD - ephem->tb);

		ephem->tb      = TxTOD;
		ephem->X_tb    = SV_state[0];
		ephem->Y_tb    = SV_state[1];
		ephem->Z_tb    = SV_state[2];

		ephem->Xdot_tb = SV_state[3];
		ephem->Ydot_tb = SV_state[4];
		ephem->Zdot_tb = SV_state[5];
	}

	return result;
}

///////////////////////////////////////////////////////////////////////////////
//
// Function : GICSRxSatRange
// Purpose  : Computes the user-satellite range, and the unit vector
//            from the user to the satellite.
// Args I   : dUserPos  (double[3])  User      position in WGS84 in meters.
//            dSVPos    (double[3])  Satellite position in WGS84 in meters.
// Args I/O :
// Args O   : unit      (float[3])   Unit vector from user to sat in WGS84.
// Returns  : double                 Range in meters.
// Depends  :
// Calls    :
// Comments : GPS ICD-200 models for satellite position computation
//            using the broadcast navigation message.
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//          |         |          |      |                                    |
//          |   1.0   | 97/04/29 | GNSS |            First version           |
//	        |   1.1   | 08/01/21 | SSKK | Sat Range calculated due to user   |
//	        |         |          |      | reception time 	   				 |
//////////////////////////////////////////////////////////////////////////////
double GICSRxSatRange(const double dUserPos[3], const double dSVPos[3], float unit[3])
{
	double range, Xga, Yga, Zga;

	// comment out earth rotation because calculation is made due to
	// user reception time
	Xga = dSVPos[0] - dUserPos[0];
	Yga = dSVPos[1] - dUserPos[1];
	Zga = dSVPos[2] - dUserPos[2];

	range = sqrt (Xga*Xga + Yga*Yga + Zga*Zga);

	if (range > 0.)
	{
		unit[0] = (float)(Xga/range);
		unit[1] = (float)(Yga/range);
		unit[2] = (float)(Zga/range);
	}
	else
	{
		unit[0] = 1.;
		unit[1] = 0.;
		unit[2] = 0.;
	}

	return range;

} // END of function GICSRxSatRange

#if CALC_ATMSPHCORR == 1

//////////////////////////////////////////////////////////////////////////////////////
//
// Function : GICSRxSatElevAzim
// Purpose  : Computes the elevation and azimuth of a GNSS satellite.
// Args I   : LatLonTrig  (float[5])  Lat and Lon cosines and sines:
//									  {cos(Lat), sin(Lat), cos(Lon), sin(Lon)}
//            Rsu         (float[3])  LoS vector in WGS84 in meters.
// Args I/O :
// Args O   : El          (float * )  Elevation angle in degrees.
//          : Az          (float * )  Azimuth   angle in degrees.
// Returns  :
// Depends  :
// Calls    :
// Comments :
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE                       |
//          |         |          |      |                                            |
//          |   1.0   | 97/04/29 | GNSS |            First  version                  |
//          |   2.0   | 12/01/17 | cmvv |            Second version                  |
//
//////////////////////////////////////////////////////////////////////////////////////
void GICSRxSatElevAzim(float LatLonTrig[4], float Rsu[3], float *El, float *Az)
{
	float aux, Nor[3], StoN[3], WtoE[3], Hor[3];

	// Remember:
	// LatLonTrig[0] = cos(UserLat)
	// LatLonTrig[1] = sin(UserLat)
	// LatLonTrig[2] = cos(UserLon)
	// LatLonTrig[3] = sin(UserLon)

	Nor[0] = LatLonTrig[0]*LatLonTrig[2];
	Nor[1] = LatLonTrig[0]*LatLonTrig[3];
	Nor[2] = LatLonTrig[1];

	StoN[0] = -LatLonTrig[2]*LatLonTrig[1];
	StoN[1] = -LatLonTrig[3]*LatLonTrig[1];
	StoN[2] = +LatLonTrig[0];

	WtoE[0] = -LatLonTrig[3];
	WtoE[1] = +LatLonTrig[2];
	WtoE[2] = 0.0;

	aux = Nor[0]*Rsu[0] + Nor[1]*Rsu[1] + Nor[2]*Rsu[2];

	Hor[0] = Rsu[0] - aux*Nor[0];
	Hor[1] = Rsu[1] - aux*Nor[1];
	Hor[2] = Rsu[2] - aux*Nor[2];

	// Elevation angle is returned in degrees:
	// ---------------------------------------
	*El = RADTODEG*(PI_G/2.0-acos(aux));

	// Azimuth angle is returned in degrees:
	// ------------------------------------
	*Az = RADTODEG*(atan2(WtoE[0]*Hor[0] + WtoE[1]*Hor[1] + WtoE[2]*Hor[2],
						  StoN[0]*Hor[0] + StoN[1]*Hor[1] + StoN[2]*Hor[2]));
}

#else

//////////////////////////////////////////////////////////////////////////////////////
//
// Function : GICSRxSatElev
// Purpose  : Computes the elevation of a GNSS satellite.
// Args I   : LatLonTrig  (float[5])  Lat and Lon cosines and sines:
//									  {cos(Lat), sin(Lat), cos(Lon), sin(Lon)}
//            Rsu         (float[3])  LoS vector in WGS84 in meters.
// Args I/O :
// Args O   : El          (float * )  Elevation angle in degrees.
// Returns  :
// Depends  :
// Calls    :
// Comments :
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE                       |
//          |         |          |      |                                            |
//          |   1.0   | 12/01/17 | cmvv |            Second version                  |
//
//////////////////////////////////////////////////////////////////////////////////////
void GICSRxSatElev(float LatLonTrig[4], float Rsu[3], float *El)
{
	float aux, Nor[3];

	Nor[0] = LatLonTrig[0]*LatLonTrig[2];
	Nor[1] = LatLonTrig[0]*LatLonTrig[3];
	Nor[2] = LatLonTrig[1];

	aux = Nor[0]*Rsu[0] + Nor[1]*Rsu[1] + Nor[2]*Rsu[2];

	// Elevation angle is returned in degrees:
	// ---------------------------------------
	*El = RADTODEG*(PI_G/2.0-acos(aux));
}
#endif

///////////////////////////////////////////////////////////////////////////////
//
// Function : GICSRxSatDoppler
// Purpose  : Computes the user-satellite velocity, and the unit vector
//            from the user to the satellite.
// Args I   : dUserVel  (double[3])  User      velocity in m/s.
//            dSVVel    (double[3])  Satellite velocity in m/s.
//            unit      (float [3])  Unit vector from user to sat in WGS84.
//            dUserPos  (double[3])  User      position  in meters.
// Args I/O :
// Args O   :
// Returns  : double        velocity in m/s
// Depends  :
// Calls    :
// Comments : GPS ICD-200 models for satellite position computation
//            using the broadcast navigation message.
// History  :
//          | VERSION |   DATE     | NAME |               CHANGE               |
//          |         |            |      |                                    |
//          |   1.0   | 08/06/2001 | rrrl |            First version           |
//			|   2.0   | 19/10/2012 | cmvv | 		  Second version           |
//
///////////////////////////////////////////////////////////////////////////////
float GICSRxSatDoppler(const double dUserVel[3], const double dSVVel[3], const float unit[3], const double dUserPos[3])
{
	float dopvel, vel[3];

	// CMVV: implementation of doppler velocity:
	// dopvel = unit*(R(At)*dSVVel - (dUserVel + OMEGA_MAT*dUserPos))

	vel[0] = dSVVel[0] - dUserVel[0] + OMEGA_EARTH*dUserPos[1];
	vel[1] = dSVVel[1] - dUserVel[1] - OMEGA_EARTH*dUserPos[0];
	vel[2] = dSVVel[2] - dUserVel[2];

	dopvel  = (vel[0]*unit[0] + vel[1]*unit[1] + vel[2]*unit[2]);

	return dopvel;
}

/******************************************************************************
 * This function is the Lagrange interpolation, it interpolates, for all
 * satellites, in the tables X[][], Y[][], Z[][] for time tIntrp. The size of
 * these tables is needed (nepoch) to avoid using zeroes at the end. Output is
 * pOut[MAXSAT][3] which are position vectors for all the satellites at time
 * tIntrp.
 *****************************************************************************/
static int GICSRxInterpolate(double tIntrp,
		int    nepoch,
		double T[MAXSP3INTERVALS],
		double X[MAXSP3INTERVALS],
		double Y[MAXSP3INTERVALS],
		double Z[MAXSP3INTERVALS],
		double Clock[MAXSP3INTERVALS],
		double pOut[4],
		double pDeltaOut[3])
{
	int i0 = 0,i,j;
	double facTop, facBot, fac[MININTERP], deltaFac[MININTERP];

	// Find interpolation point in time table.
	while(i0 < nepoch && T[i0+ (int)(MININTERP/2)]<tIntrp)
		i0++;

	while (i0 + MININTERP -1 >= nepoch)
		i0--;

	if (i0 < 0 || nepoch < MININTERP) {
		return 1;
	}

	// Compute interpolation coefficients.
	for (i = 0; i < MININTERP; i++)
	{
		facBot = facTop = 1.0;
		deltaFac[i] = 0;

		for(j = 0; j < MININTERP; j++)
		{
			if (i != j)
			{
				if(((double)(i-j)) * (T[i0+i] -T[i0+j]) <= 0)
				{
					// logMessage("Warning: Inconsistency detected while interpolating. "
					// "Probably input SP3-C files overlap",
					//  false,false,false);
					return(-1);
				}

				facTop *= tIntrp  - T[i0+j];
				facBot *= T[i0+i] - T[i0+j];
			}
		}

		fac[i] = facTop/facBot;

		for(j = 0; j < MININTERP; j++)
		{
			if (i != j)
			{
				deltaFac[i] += (fac[i]/(tIntrp - T[i0+j]));
			}
		}
	}

	//  Interpolate vectors for all satellites.

	pOut[0] = pOut[1] = pOut[2] = pOut[3] = 0.0;
	pDeltaOut[0] = pDeltaOut[1] = pDeltaOut[2] = 0.0;

	for (j = 0; j < MININTERP; j++)
	{
		pOut[0] += fac[j] * X[i0+j];
		pOut[1] += fac[j] * Y[i0+j];
		pOut[2] += fac[j] * Z[i0+j];
		pOut[3] += fac[j] * Clock[i0+j];

		pDeltaOut[0] += deltaFac[j]* X[i0+j];
		pDeltaOut[1] += deltaFac[j]* Y[i0+j];
		pDeltaOut[2] += deltaFac[j]* Z[i0+j];
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
//
// Function : GICSRxEastNorthUpToXYZ
// Purpose  : Change velocity from north, east, up into WGS84 coord system
// Args I   : DOUBLE     lat       Latitude  (rad)
//            DOUBLE     lon       Longitude (rad)
//            DOUBLE     h         Ellipsoidal height (m)
//                                             DOUBLE     east    velocity_east  (m/s)
//            DOUBLE     north   velocity_north (m/s)
//            DOUBLE     up      velocity_up    (m/s)
// Args I/O :
// Args O   : DOUBLE velocity_XYZ[3]    velocity in WGS84 (m)
// Returns  :
// Depends  :
// Calls    :
// Comments :
// History  :
//          |  Version   |     Date     |   Name   |        Change        |
//          |            |              |          |                      |
//
///////////////////////////////////////////////////////////////////////////////
static void GICSRxEastNorthUpToXYZ(double Lat, double Lon,
		double east, double north, double up, double XYZ[3])
{
	double ENU2XYZ[3][3] = {{-sin(Lon), -sin(Lat)*cos(Lon), cos(Lat)*cos(Lon)},
							{ cos(Lon), -sin(Lat)*sin(Lon), cos(Lat)*sin(Lon)},
			{0, cos(Lat), sin(Lat)}};

	double ENU[3] = {east, north, up};
	int i, k;

	for(i = 0; i < 3; i++)
	{
		XYZ[i] = 0;
		for(k = 0; k < 3; k++)
		{
			XYZ[i] += ENU2XYZ[i][k] * ENU[k];
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
//
// Function : GICSRxPos_SP3
// Purpose  : Computes Position, Velocity in WGS84 and
//            eccentric anomaly for a GPS satellite using SP3 ephemeris.
//
// Args I   :  TxTime   (double)          Desired time for output data.
//	      	   RxTime   (double)          User reception time.
//			   pSp3		(SP3STRUCT_D*)    Pointer to SP3 ephemeris.
//			   PRN      (int)             PRN of the satellite (multi-constellation)
// Args I/O :
// Args O   :
// Returns  : char                   	  1 if everything is OK
// Depends  :
// Calls    :
// Comments : GPS SP3 model for satellite position computation
//            using the broadcast navigation message.
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//
///////////////////////////////////////////////////////////////////////////////
char GICSRxSatPos_SP3(SP3STRUCT_D *pSp3, int PRN, double RxTime, double TxTime, ephemsat_t *pephsat)
{
	int i = 0, k = pSp3->psiPrnToIndex[PRN];

	double dt = 0;
	double RelCorr;
	double aux[2];

	double  geod_pos[3];

	double  NGS_INTERP[4]     = {0,0,0,0},
		    NGS_INTERP_m[3]   = {0,0,0},
			NGS_INTERP_Vel[3] = {0,0,0},
			XYZant[3] 		  = {0,0,0};

	if(GICSRxInterpolate(TxTime,pSp3->nrec[k], pSp3->timeTable[k],
			pSp3->Pos_X_Table[k], pSp3->Pos_Y_Table[k], pSp3->Pos_Z_Table[k],
			pSp3->Clock_Table[k], NGS_INTERP, NGS_INTERP_Vel) != 0)
	{
		return 0;
	}

	NGS_INTERP_m[0] = NGS_INTERP[0]*1000;
	NGS_INTERP_m[1] = NGS_INTERP[1]*1000;
	NGS_INTERP_m[2] = NGS_INTERP[2]*1000;

	// Compute the unit vector: Centre of Earth -> Satellite
	ECEFtoNAV_pos(NGS_INTERP_m, geod_pos);

	GICSRxEastNorthUpToXYZ(geod_pos[0], geod_pos[1], 0, 0, 0, XYZant);

	for(i = 0; i < 3; i++)
	{
		NGS_INTERP[i] = NGS_INTERP[i] - XYZant[i]/1000.;
	}

	// WGS-84 Position coordinates
	pephsat->pos[0] = NGS_INTERP[0]*1000;
	pephsat->pos[1] = NGS_INTERP[1]*1000;
	pephsat->pos[2] = NGS_INTERP[2]*1000;

	// approxTransit = PR/c + ck_bias_receptor - ck_bias_satelite
	dt = RxTime - TxTime + NGS_INTERP[3]*1e-6;

	aux[0] = pephsat->pos[0] + pephsat->pos[1] * OMEGA_EARTH * dt;
	aux[1] = pephsat->pos[1] - pephsat->pos[0] * OMEGA_EARTH * dt;

	pephsat->pos[0] = aux[0];
	pephsat->pos[1] = aux[1];

	// Inertial Velocity in WGS-84 coordinates
	pephsat->vel[0] = NGS_INTERP_Vel[0] * 1000;
	pephsat->vel[1] = NGS_INTERP_Vel[1] * 1000;
	pephsat->vel[2] = NGS_INTERP_Vel[2] * 1000;

	aux[0] = pephsat->vel[0] - OMEGA_EARTH*pephsat->pos[1] - OMEGA_EARTH*OMEGA_EARTH*dt*pephsat->pos[0];
	aux[1] = pephsat->vel[1] + OMEGA_EARTH*pephsat->pos[0] - OMEGA_EARTH*OMEGA_EARTH*dt*pephsat->pos[1];

	pephsat->vel[0] = aux[0];
	pephsat->vel[1] = aux[1];

	// Calculate satellite clock bias and drift.
	//
	// Relativistic correction.
	RelCorr = -2*(pephsat->pos[0]*pephsat->vel[0] + pephsat->pos[1]*pephsat->vel[1] + pephsat->pos[2]*pephsat->vel[2])/SPEED_OF_LIGHT;

	// Satellite Clock Bias.
	pephsat->clk_bias = NGS_INTERP[3]*1e-6*SPEED_OF_LIGHT + RelCorr;

	// Satellite Clock Drift.
	pephsat->clk_drift = 0;

	// Satellite Clock Drift Rate.
	pephsat->clk_driftrate = 0;

	return 1;

} // GICSRxSatPos_SP3

char GICSRxReadSP3File(SP3STRUCT_D *sp3_P, char *sp3_filename)
{
	int i			  = 0,
		j			  = 0,
		k			  = 0,
		SVN			  = 0,
		svnIndex	  = 0,
		intervalIndex = 0,
		year		  = 0,
		month		  = 0,
		day			  = 0,
		hour		  = 0,
		min	          = 0,
		siNEpoch	  = 0,
		WN			  = 0,
		bEndOfFile    = 0,
		deltaIndexSp3 = 0,
		IDsat[MAX_SATS_INSP3_SI],
		*psiPRNRead   = NULL,
		psiPosStatus[MAX_SATS_INSP3_SI][MAXSP3INTERVALS];

	double  sec     = 0,
			secGPS  = 0,
			stepSec = 0,
			TOW		= 0,
			X[MAX_SATS_INSP3_SI][MAXSP3INTERVALS],
			Y[MAX_SATS_INSP3_SI][MAXSP3INTERVALS],
			Z[MAX_SATS_INSP3_SI][MAXSP3INTERVALS],
			C[MAX_SATS_INSP3_SI][MAXSP3INTERVALS],
			T[MAXSP3INTERVALS]; // MJD (GPS time)

	char    sLine[200],
			SatType = 0;

    long    GPSweek = 0;

	FILE   *sp3File = fopen(sp3_filename, "rt"); // SP3 file descriptor.

	/************************/
	/* INITIALIZE VARIABLES */
	/************************/

	for (i = 0;i < 200; i++)
	{
		sLine[i] = 0;
	}

	for (i=0; i < MAX_SATS_INSP3_SI; i++)
	{
		IDsat[i]   = 0;
	}

	for (i = 0; i < MAX_SATS_INSP3_SI; i++)
	{
		for (j=0; j < MAXSP3INTERVALS; j++)
		{
			X[i][j]            = 0;
			Y[i][j]            = 0;
			Z[i][j]            = 0;
			psiPosStatus[i][j] = 0;
			C[i][j]            = 0;
		}
	}

	for (j = 0; j < MAXSP3INTERVALS; j++)
	{
		T[j] = 0;
	}

	for (i=0; i < MAX_SATS_INSP3_SI; i++)
	{
		sp3_P->psiPrnToIndex[i] = 0;
	}

	/***************/
	/* READ HEADER */
	/***************/

	// Read Header line 1
	fgets(sLine,SP3_LINE_LENGTH_SI,sp3File);
	sscanf(&sLine[3],"%d %d %d %d %d %lf %d", &year, &month, &day, &hour, &min, &sec, &siNEpoch);

	CalToGPS_G(year, month, day, hour, min, sec, &GPSweek, &TOW);

	sp3_P->dStartAGPS = TOW ;

	// Read Header line 2
	fgets(sLine, SP3_LINE_LENGTH_SI, sp3File);
	sscanf(&sLine[3],"%d %lf %lf",&WN,&secGPS,&stepSec);

	sp3_P->dStepTime = stepSec;

	// Read Header lines 3 to 7

	// Satellite type and number
	// There are five lines with the satellites included in the file.
	fgets (sLine,SP3_LINE_LENGTH_SI,sp3File);
	sscanf(&sLine[3],"%d",&(sp3_P->siNumberSats));

	psiPRNRead = IDsat;

	for (i = 0; i < MAX_SATS_INSP3_SI; i++)
	{
		IDsat[i] = 0;
	}

	for (i=0; i < 5; i++)
	{
		for (j = 9; j < 60; j+=3)
		{
			sscanf(&sLine[j],"%c%2d", &SatType, psiPRNRead);

			if (SatType=='R')
			{
				*psiPRNRead += 37;
			}

			if (SatType=='E')
			{
				*psiPRNRead += 61;
			}


			if (SatType=='1')
			{
				*psiPRNRead += 100;
			}

			if (SatType=='2')
			{
				*psiPRNRead += 200;
			}

			psiPRNRead++;
		}

		fgets(sLine,SP3_LINE_LENGTH_SI,sp3File);
	}

	// Check if all required satellites are present in the file
	// and store in the PRN to index mask array (PRN is stored in IDsat[prnIndex])
	for (svnIndex=0; svnIndex < sp3_P->siNumberSats; svnIndex++)
	{
		sp3_P->psiPrnToIndex[IDsat[svnIndex]] = svnIndex;
	}

	// Read the Position and Clock of the different satellites.
	bEndOfFile    = 0;

	intervalIndex = -1;
	while (!bEndOfFile)
	{
		if (fgets(sLine,SP3_LINE_LENGTH_SI,sp3File) == NULL)
		{
			bEndOfFile = 1;
		}
		else if (strncmp(sLine,"EOF",3) == 0)
		{
			bEndOfFile = 1;
		}
		else if ((0 == strncmp(sLine,"*  19",5)) || (0 == strncmp(sLine,"*  20",5)))
		{
			intervalIndex++;

			if (intervalIndex >= MAXSP3INTERVALS)
			{
				fprintf(stderr, "Number of NGS registers is larger than MAXSP3INTERVALS macro \n");
				return(-1);
			}

			sscanf(&sLine[3],"%d %d %d %d %d %lf", &year, &month, &day, &hour, &min, &sec);

			CalToGPS_G(year, month, day, hour, min, sec, &GPSweek, &TOW);

			T[intervalIndex] = TOW;

			svnIndex = 0;

			for (i = 0;i<MAX_SATS_INSP3_SI;i++)
			{
				sp3_P->validSP3Flag[i] = 63;
			}

			for (i=0; i < sp3_P->siNumberSats; i++)
			{
				if(sLine[0] != 'P')
				{
					fgets(sLine,SP3_LINE_LENGTH_SI,sp3File);
				}

				if ('P' == sLine[0])
				{
					if(sLine[1]=='G')  // GPS
					{
						deltaIndexSp3 = 0;
					}
					else if(sLine[1] == 'R')
					{
						deltaIndexSp3 = 37;
					}
					else if(sLine[1] == 'E')
					{
						deltaIndexSp3 = 61;
					}

					if (0 != IDsat[i])
					{
						sscanf(&sLine[2],"%d %lf %lf %lf %lf",
							   &SVN, &X[svnIndex][intervalIndex],&Y[svnIndex][intervalIndex],&Z[svnIndex][intervalIndex], &C[svnIndex][intervalIndex]);

						if (( X[svnIndex][intervalIndex] != BAD_POS_VALUE_SP3) &&
							( Y[svnIndex][intervalIndex] != BAD_POS_VALUE_SP3) &&
							( Z[svnIndex][intervalIndex] != BAD_POS_VALUE_SP3) &&
							( C[svnIndex][intervalIndex] < BAD_CLK_VALUE_SP3))
						{
							psiPosStatus[svnIndex][intervalIndex]    = 1;
							sp3_P->validSP3Flag[SVN-1+deltaIndexSp3] = 1;
						}
						svnIndex++;
					}
				} // end checking of Position Line.

				fgets(sLine, SP3_LINE_LENGTH_SI, sp3File);
			}
		}
	}
	fclose(sp3File);

	sp3_P->siNRec = intervalIndex+1;

	for( k = 0; k < svnIndex; k++)
	{
		sp3_P->nrec[k]=0;

		for(j=0;j<sp3_P->siNRec;j++)
		{
			if(psiPosStatus[k][j])
			{
				sp3_P->timeTable[k][sp3_P->nrec[k]]   = T[j];

				sp3_P->Pos_X_Table[k][sp3_P->nrec[k]] = X[k][j];
				sp3_P->Pos_Y_Table[k][sp3_P->nrec[k]] = Y[k][j];
				sp3_P->Pos_Z_Table[k][sp3_P->nrec[k]] = Z[k][j];

				sp3_P->Clock_Table[k][sp3_P->nrec[k]] = C[k][j];

				sp3_P->nrec[k]++;
			}
		}
	}

	return 0;
}

