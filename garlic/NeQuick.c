//*******************************************************************************************//
// Copyright: GMV, S.A. 2015. Property of GMV, S.A.; all rights reserved.
//
// Project:                  PRESENCE2
// File:                     NeQuick.h
// Version:                  1.0
// Date (YY/MM/DD):          2015/04/27
// Component:
//
// Purpose:   NeQuick G ionosphere prediction model implementation according to the document
//            "Ionospheric Correction Algorithm for Galileo Single Frequency Users"
//
// Language:  C
//
// History:
//   1.0   -   2015/04/27   -   CMVV   -   Initial Release
//*******************************************************************************************//

#include <stdio.h>
#include <math.h>

#include "NeQuick.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Declaration of internal functions used in the NeQuick TEC computation algorithm.

static char   NeqCheckInputs(NeQuickInputData_st *pstNeQuickInputData);
static void   DoTECIntegration(IntegrateData_st *pstIntegrateData, char bVert, SPoint_st *stP0, double *pdNmax, LayerProperties_st *pstLayers, double *pdTEC);
static double NeqCalcModip(double dLat, double dLng, MODIP_st *pstModip);
static double NeqInterpolate(double pdZ[4], double dDeltaX);
static double NeqModipToAz(double dModip, int siNumCoeff, double *pdCoeff);
static char   NeqGetRayProperties(SPoint_st *pstP1, SPoint_st *pstP2, SPoint_st *pstRay, double *pdZeta, double *pdSinSig, double *pdCosSig);
static void   NeqCalcRayProperties1(SPoint_st *pstP1, SPoint_st *pstP2, SPoint_st *pstRay, double *pdZeta);
static void   NeqCalcRayProperties2(SPoint_st *pstP1, SPoint_st *pstP2, SPoint_st *pstRay, double *pdSinSig, double *pdCosSig);
static double NeqIntegrate(IntegrateData_st *pstIntegrateData, double dH1, double dH2, int siCurrentLevel, SPoint_st *pstPactual, double *pdNmax, LayerProperties_st *pstLayers);
static double NeqGetNeOnVertRay(double dH, LayerProperties_st *pstLayers, double *pdNmax);
static double NeqGetNeOnSlantRay(double dS, NeQuickInputData_st *pstNeQuickInputData, GeometryData_st *pstGeom, double *pdNmax, SPoint_st *pstPactual, LayerProperties_st *pstLayers, CurrentCCIR_st *pstCurrCCIR);
static void   NeqCalcLLHOnRay(double dS, SPoint_st *pstRay, SPoint_st *pstP1, double dSinSig, double dCosSig, SPoint_st *pstPactual);
static void   NeqCalcEpstParams(NeQuickInputData_st *pstNeQuickInputData, SPoint_st *pstPactual, double dSinDelta, double dCosDelta, double *pdNmax, LayerProperties_st *pstLayers, CurrentCCIR_st *pstCurrCCIR);
static void   NeqCalcSphLegCoeffs(double dUT, CurrentCCIR_st *pstCurrCCIR);
static double NeqGetF2FreqFromCCIR(double dCosLat, double dLng, double *pdLegCoeffs, double pdSinModipToN[12], int siMode);
static double NeqCriticalFreqToNe(double dF0);
static double NeqCalcF2PeakHeight(double dM3000, double dF0E, double dF0F2);
static double NeqEpstein(double dNmax, double dHmax, double dB, double dH);
static double NeqCalcTopsideNe(double dH, LayerProperties_st *pstLayers, double *pdNmax);
static double NeqCalcBottomsideNe(double dHH, LayerProperties_st *pstLayers);
static double NeqJoin(double dF1, double dF2, double dAlpha, double dX);
static double NeqClipExp(double dPower);
static double NeqSquared(double dValue);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// NeQuick G algorithm. Calculates the Slant Total Electron Content (STEC) taking as
// input data the user and satellite positions, plus broadcast ionosphere parameters.
// Return an element of the NeQuick_status enumeration (OK or ERROR).
NeQuick_status NeQuick(IntegrateData_st *pstIntegrateData, double *pdSTEC)
{
	char bError;
	NeQuick_status retval = E_ERROR;
	
	// Internal data structure.
	LayerProperties_st pstLayers;

	// Call NeqCheckInputs to check that inputs are within ranges.
	bError = NeqCheckInputs(&pstIntegrateData->pstNeQuick);
	
	// No error is detected in the inputs.
	if(bError == 0)
	{
		// Set P1 (receiver position) and P2 (satellite position).
		pstIntegrateData->pstGeom.stP1.dLat = pstIntegrateData->pstNeQuick.pdGssPosLLH[0];
		pstIntegrateData->pstGeom.stP1.dLng = pstIntegrateData->pstNeQuick.pdGssPosLLH[1];
		pstIntegrateData->pstGeom.stP1.dH   = pstIntegrateData->pstNeQuick.pdGssPosLLH[2];
		pstIntegrateData->pstGeom.stP1.dR   = pstIntegrateData->pstGeom.stP1.dH + RE;
		
		pstIntegrateData->pstGeom.stP2.dLat = pstIntegrateData->pstNeQuick.pdSatPosLLH[0];
		pstIntegrateData->pstGeom.stP2.dLng = pstIntegrateData->pstNeQuick.pdSatPosLLH[1];
		pstIntegrateData->pstGeom.stP2.dH   = pstIntegrateData->pstNeQuick.pdSatPosLLH[2];
		pstIntegrateData->pstGeom.stP2.dR   = pstIntegrateData->pstGeom.stP2.dH + RE;
		
		// Call function NeqCalcModip to calculate modified dip latitude at the receiver position P1.
		double dModipRx = NeqCalcModip(pstIntegrateData->pstGeom.stP1.dLat, pstIntegrateData->pstGeom.stP1.dLng, &pstIntegrateData->pstNeQuick.pstModip);

		// Store modified dip latitude. Used later to calculate the geo-magnetic latitude and the variance associated to the
		// ionospheric delay correction, based on the algorithm described in MOPS.
		pstIntegrateData->dModipRx = dModipRx;

		// Calculate F10.7 solar flux.
		pstIntegrateData->pstNeQuick.dAzBase = NeqModipToAz(dModipRx, pstIntegrateData->pstNeQuick.siNumCoeff, pstIntegrateData->pstNeQuick.pdCoeff);
		
		// Check that solar flux is within ranges. If not, define maximum and minimum values.
		if(pstIntegrateData->pstNeQuick.pdCoeff[0] == 0 && pstIntegrateData->pstNeQuick.pdCoeff[1] == 0 && pstIntegrateData->pstNeQuick.pdCoeff[2] == 0)
		{
			pstIntegrateData->pstNeQuick.dAzBase = 63.7;
		}
		if(pstIntegrateData->pstNeQuick.dAzBase < 0)
		{
			pstIntegrateData->pstNeQuick.dAzBase = 0;
		}
		if(pstIntegrateData->pstNeQuick.dAzBase > 400)
		{
			pstIntegrateData->pstNeQuick.dAzBase = 400;
		}
		
		// Calculate ray properties and check if it is valid.
		bError = NeqGetRayProperties(&pstIntegrateData->pstGeom.stP1, &pstIntegrateData->pstGeom.stP2, &pstIntegrateData->pstGeom.stRay, &pstIntegrateData->pstGeom.dZeta, &pstIntegrateData->pstGeom.dSinSig, &pstIntegrateData->pstGeom.dCosSig);
		
		// If ray is valid.
		if(bError == 0)
		{
			// Calculate slant distance of each point.
			pstIntegrateData->pstGeom.stP1.dS = sqrt(NeqSquared(pstIntegrateData->pstGeom.stP1.dR) - NeqSquared(pstIntegrateData->pstGeom.stRay.dR));
			pstIntegrateData->pstGeom.stP2.dS = sqrt(NeqSquared(pstIntegrateData->pstGeom.stP2.dR) - NeqSquared(pstIntegrateData->pstGeom.stRay.dR));
			
			// Start point for integration.
			SPoint_st stP0;
			
			stP0.dLat = pstIntegrateData->pstGeom.stP1.dLat;
			stP0.dLng = pstIntegrateData->pstGeom.stP1.dLng;
			if(pstIntegrateData->pstGeom.stP1.dH < 0)
			{
				stP0.dH = 0;
				stP0.dR = RE;
			}
			else
			{
				stP0.dH = pstIntegrateData->pstGeom.stP1.dH;
				stP0.dR = pstIntegrateData->pstGeom.stP1.dR;
			}
			stP0.dS = sqrt(NeqSquared(stP0.dR) - NeqSquared(pstIntegrateData->pstGeom.stRay.dR));
			
			// Initialize current position.
			pstIntegrateData->pstGeom.stPactual.dLat = pstIntegrateData->pstGeom.stP1.dLat;
			pstIntegrateData->pstGeom.stPactual.dLng = pstIntegrateData->pstGeom.stP1.dLng;

			pstIntegrateData->pstGeom.stPactual.dH   = pstIntegrateData->pstGeom.stP1.dH;
			pstIntegrateData->pstGeom.stPactual.dR   = pstIntegrateData->pstGeom.stP1.dH + RE;
			pstIntegrateData->pstGeom.stPactual.dS   = sqrt(NeqSquared(pstIntegrateData->pstGeom.stPactual.dR) - NeqSquared(pstIntegrateData->pstGeom.stRay.dR));

			// Calculate sine and cosine of delta (solar declination)
			double amrad = (0.9856*(pstIntegrateData->pstNeQuick.siMonth*30.5 - 15 + (18 - pstIntegrateData->pstNeQuick.dUT)/24) - 3.289)*DR;
			
			pstIntegrateData->pstGeom.dSinDelta = (0.39782*sin(amrad + (1.916*sin(amrad) + 0.02*sin(2*amrad) + 282.634)*DR));
			pstIntegrateData->pstGeom.dCosDelta = sqrt(1 - NeqSquared(pstIntegrateData->pstGeom.dSinDelta));

			double dNmax = 0, dTEC = 0;
			if(pstIntegrateData->pstGeom.stRay.dR < 0.1)
			{
				// Calculate ionosphere parameters for vertical ray.
				NeqCalcEpstParams(&pstIntegrateData->pstNeQuick, &pstIntegrateData->pstGeom.stPactual, pstIntegrateData->pstGeom.dSinDelta, pstIntegrateData->pstGeom.dCosDelta, &dNmax, &pstLayers, &pstIntegrateData->pstCurrCCIR);
				pstIntegrateData->bVert = 1;
			}
			else
			{
				pstIntegrateData->bVert = 0;
			}

			// Call function DoTECIntegration to perform TEC integration along ray (dTEC).
			DoTECIntegration(pstIntegrateData, pstIntegrateData->bVert, &stP0, &dNmax, &pstLayers, &dTEC);
			
			retval = E_OK;
			
			// Convert internal TEC value to correct units for output from NeQuick function, factor of 1000 since integration is done based in km.
			*pdSTEC = 1000*dTEC;
		}
		else
		{
			retval = E_ERROR;
		}
	}
	else
	{
		retval = E_ERROR;
	}
	return retval;
}

static char NeqCheckInputs(NeQuickInputData_st *pstNeQuickInputData)
{
	// This function checks some of the input data to NeQuick to determine if values are
	// within range and will allow a valid TEC value to be computed or not. Note that
	// MODIP values, CCIR maps and Kronrod tolerances are not checked because in
	// NeQuick G these values should already be checked before being passed to NeQuick.

	char bError;
	
	// User and satellite latitudes.
	double lat_u = pstNeQuickInputData->pdGssPosLLH[0];
	double lat_s = pstNeQuickInputData->pdSatPosLLH[0];
	
	// Latitude flags.
	char user_lat_check = (lat_u >= -90 && lat_u <= +90);
	char ssat_lat_check = (lat_s >= -90 && lat_s <= +90);
	
	// Time flags.
	char month_check = (pstNeQuickInputData->siMonth >= 1 && pstNeQuickInputData->siMonth <= 12);
	char utime_check = (pstNeQuickInputData->dUT     >= 0 && pstNeQuickInputData->dUT     <= 24);
	
	// Coefficients vector flag.
	char coeff_check = (pstNeQuickInputData->siNumCoeff >= 1 && pstNeQuickInputData->pdCoeff != NULL);
	
	// Check whether all the flags are valid. Otherwise set error flag to 1.
	if(user_lat_check && ssat_lat_check && month_check && utime_check && coeff_check)
	{
		bError = 0;
	}
	else
	{
		bError = 1;
	}
	return bError;
}

static void DoTECIntegration(IntegrateData_st *pstIntegrateData, char bVert, SPoint_st *stP0, double *pdNmax, LayerProperties_st *pstLayers, double *pdTEC)
{
	// This function checks whether the ray is vertical or slant, and where the start and
	// end points are located in regards to the different integration points, before passing
	// the appropriate information to the integration function.

	const double H1a = 1000;
	const double H1b = 2000;
	
	double height_P0, height_P1, height_P2, S1a, S1b;

	// Get height values that define the integration ranges.
	if(bVert == 1)
	{
		height_P0 = stP0->dH;
		height_P1 = pstIntegrateData->pstGeom.stP1.dH;
		height_P2 = pstIntegrateData->pstGeom.stP2.dH;

		// Slant distances are equal to the H1a and H1b heights for vertical ray.
		S1a = H1a;
		S1b = H1b;
	}
	else
	{
		height_P0 = stP0->dS;
		height_P1 = pstIntegrateData->pstGeom.stP1.dS;
		height_P2 = pstIntegrateData->pstGeom.stP2.dS;

		// Get slant distances using the radius of the perigee ray.
		S1a = sqrt(NeqSquared(H1a + RE) - NeqSquared(pstIntegrateData->pstGeom.stRay.dR));
		S1b = sqrt(NeqSquared(H1b + RE) - NeqSquared(pstIntegrateData->pstGeom.stRay.dR));
	}
	
	// Initialize TEC value to 0 before integration to add different ionosphere contributions
	*pdTEC = 0;
	
	// Get pointer to stPactual point.
	SPoint_st *pstPactual = &pstIntegrateData->pstGeom.stPactual;

	// Check if ray path crosses either of the integration break points and split up integration accordingly
	// (it is assumed that P1 is always lower than P2).
	//
	// Start and end point both below 1st break point.
	if(pstIntegrateData->pstGeom.stP2.dH <= H1a)
	{
		// Set tolerance for heights below 1st break point.
		pstIntegrateData->dTolerance = pstIntegrateData->pstNeQuick.pdKronrodTol[0];
		// Call NeqIntegrate with start point of height_P0, end point of height_P2.
		*pdTEC += NeqIntegrate(pstIntegrateData, height_P0, height_P2, 1, pstPactual, pdNmax, pstLayers);
	}
	else
	{
		// Set slant distance of 1st integration breakpoint S1a at 1000km.
		//
		// End point below 2nd break point.
		if(pstIntegrateData->pstGeom.stP2.dH <= H1b)
		{
			// Start point above 1st break point.
			if(pstIntegrateData->pstGeom.stP1.dH >= H1a)
			{
				// Start and end points are both between 1st and 2nd break points.

				// Set tolerance for heights above 1st break point.
				pstIntegrateData->dTolerance = pstIntegrateData->pstNeQuick.pdKronrodTol[1];
				// Call NeqIntegrate with start point of height_P1, end point of height_P2.
				*pdTEC += NeqIntegrate(pstIntegrateData, height_P1, height_P2, 1, pstPactual, pdNmax, pstLayers);
			}
			else
			{
				// Ray path crosses 1st integration break point.
				// Sum the TEC values from the two NeqIntegrate calls.

				// Set tolerance for heights below 1st break point.
				pstIntegrateData->dTolerance = pstIntegrateData->pstNeQuick.pdKronrodTol[0];
				// Call NeqIntegrate with start point of height_P0, end point S1a.
				*pdTEC += NeqIntegrate(pstIntegrateData, height_P0, S1a, 1, pstPactual, pdNmax, pstLayers);

				// Set tolerance for heights above 1st break point.
				pstIntegrateData->dTolerance = pstIntegrateData->pstNeQuick.pdKronrodTol[1];
				// Call NeqIntegrate with start point of S1a, end point of height_P2.
				*pdTEC += NeqIntegrate(pstIntegrateData, S1a, height_P2, 1, pstPactual, pdNmax, pstLayers);
			}
		}
		else
		{
			// Start point above 2nd break point.
			if(pstIntegrateData->pstGeom.stP1.dH >= H1b)
			{
				// Start and end points are both above 2nd break point.

				// Set tolerance for heights above 1st break point.
				pstIntegrateData->dTolerance = pstIntegrateData->pstNeQuick.pdKronrodTol[1];
				// Call NeqIntegrate with start point of height_P1, end point of height_P2.
				*pdTEC += NeqIntegrate(pstIntegrateData, height_P1, height_P2, 1, pstPactual, pdNmax, pstLayers);
			}
			else
			{
				// Set slant distance of 2nd integration breakpoint S1b at 2000km.
				if(pstIntegrateData->pstGeom.stP1.dH >= H1a)
				{
					// Ray path crosses 2nd integration break point.
					// Sum the TEC values from the two NeqIntegrate calls.

					// Set tolerance for heights above 1st break point.
					pstIntegrateData->dTolerance = pstIntegrateData->pstNeQuick.pdKronrodTol[1];

					// Call NeqIntegrate with start point of height_P1, end point of S1b.
					*pdTEC += NeqIntegrate(pstIntegrateData, height_P1, S1b, 1, pstPactual, pdNmax, pstLayers);
					// Call NeqIntegrate with start point of S1b, end point of height_P2.
					*pdTEC += NeqIntegrate(pstIntegrateData, S1b, height_P2, 1, pstPactual, pdNmax, pstLayers);
				}
				else
				{
					// Ray path crosses 1st and 2nd integration break points.
					// Sum the TEC values from the three NeqIntegrate calls.

					// Set tolerance for heights below 1st break point.
					pstIntegrateData->dTolerance = pstIntegrateData->pstNeQuick.pdKronrodTol[0];
					// Call NeqIntegrate with start point of height_P0, end point of S1a.
					*pdTEC += NeqIntegrate(pstIntegrateData, height_P0, S1a, 1, pstPactual, pdNmax, pstLayers);

					// Set tolerance for heights above 1st break point.
					pstIntegrateData->dTolerance = pstIntegrateData->pstNeQuick.pdKronrodTol[1];

					// Call NeqIntegrate with start point of S1a, end point of S1b.
					*pdTEC += NeqIntegrate(pstIntegrateData, S1a, S1b, 1, pstPactual, pdNmax, pstLayers);
					// Call NeqIntegrate with start point of S1b, end point of height_P2.
					*pdTEC += NeqIntegrate(pstIntegrateData, S1b, height_P2, 1, pstPactual, pdNmax, pstLayers);
				}
			}
		}
	}
}

static double NeqCalcModip(double dLat, double dLng, MODIP_st *pstModip)
{
	// This function uses the current latitude and longitude to calculate the corresponding
	// modified dip latitude (MODIP) value. The MODIP grip should be ‘pre-wrapped’ at
	// edges and poles.
	double dResult;
	
	if(dLat <= -90)
	{
		dResult = -90;
	}
	else if(dLat >= +90)
	{
		dResult = 90;
	}
	else
	{
		// Define properties of MODIP grid.
		const double lngp = 36, dlatp = 5, dlngp = 10;
		
		// Obtain lon grid square (sj) and position in that square (dj).
		double lng1 = (dLng + 180)/dlngp;
		
		int    sj = floor(lng1) - 2;
		double dj = lng1 - floor(lng1);
		
		// Adjust for sign and wrap to grid if required.
		if(sj < 0)
		{
			sj = sj + lngp;
		}
		if(sj > (lngp - 3))
		{
			sj = sj - lngp;
		}
		
		// obtain lat grid square (si) and position in that square (di).
		double lat1 = (dLat + 90)/dlatp + 1;
		
		int    si = floor(lat1 - 1e-6) - 2;
		double di = lat1 - si - 2;
		
		// Interpolate across lat grid to obtain values at 4 lon points on lat line.
		double z1[4], z[4];

		int j, k;
		for(k = 1; k <= 4; k++)
		{
			for(j = 1; j <= 4; j++)
			{
				z1[j-1] = pstModip->pdModip[si+j][sj+k+1];
			}
			z[k-1] = NeqInterpolate(z1, di);
		}
		// Interpolate for lon value using these 4 points.
		dResult = NeqInterpolate(z, dj);
	}
	return dResult;
}

static double NeqInterpolate(double pdZ[4], double dDeltaX)
{
	// This function performs third order interpolation. It is used when calculating the modified
	// dip latitude value. Input z[4] is -1,0,1,2 values, x is position to interpolate to.
	double dIntZ = 0;
	
	if(dDeltaX < 1e-10)
	{
		dIntZ = pdZ[1];
	}
	else
	{
		double g0 = (pdZ[2] + pdZ[1]);
		double g1 = (pdZ[2] - pdZ[1]);
		double g2 = (pdZ[3] + pdZ[0]);
		double g3 = (pdZ[3] - pdZ[0])/3;
		
		double a[4];
		
		a[0] = 9*g0 - g2;
		a[1] = 9*g1 - g3;
		a[2] =   g2 - g0;
		a[3] =   g3 - g1;
		
		double Ax = 2*dDeltaX - 1;
		
		int k;
		for(k = 3; k >= 0; k--)
		{
			dIntZ = dIntZ*Ax + a[k];
		}
		dIntZ = dIntZ/16;
	}

	return dIntZ;
}

static double NeqModipToAz(double dModip, int siNumCoeff, double *pdCoeff)
{
	// This function calculates Az from the provided coefficients and modified dip latitude
	// value. If only one non-zero coefficient (a0) is provided then Az = a0.
	double dFlx = 0;
	
	int k;
	for(k = siNumCoeff-1; k >= 0; k--)
	{
		dFlx = dFlx*dModip + pdCoeff[k];
	}
	
	if(dFlx < 0)
	{
		dFlx = 0;
	}
	
	return dFlx;
}

static char NeqGetRayProperties(SPoint_st *pstP1, SPoint_st *pstP2, SPoint_st *pstRay, double *pdZeta, double *pdSinSig, double *pdCosSig)
{
	// This function obtains the properties of the ray and checks if it is a valid ray. 'Ray' is a
	// straight line passing trough p1 and p2. pstRay->dLat, pstRay->dLng and pstRay->dR are co-ordinates
	// of ray perigee, i.e. point on ray closest to center of Earth.
	char bError = 0;
	
	// Check if ray is vertical (point 2 is directly above point 1).
	if(fabs(pstP2->dLat - pstP1->dLat) < 1e-5 && fabs(pstP2->dLng - pstP1->dLng) < 1e-5)
	{
		// Set point 2 longitude to be exactly the same as point 1 longitude.
		pstP2->dLng = pstP1->dLng;
	}
	
	// Calculate ranges of points 1 and 2 from the center of the earth.
	pstP1->dR = pstP1->dH + RE;
	pstP2->dR = pstP2->dH + RE;
	
	// Call function NeqCalcRayProperties1 to calculate properties of the ray itself.
	NeqCalcRayProperties1(pstP1, pstP2, pstRay, pdZeta);
	
	// Check if ray is invalid.
	if(fabs(*pdZeta) > 90 && pstRay->dR < RE)
	{
		bError = 1;
	}
	
	// Check whether ray is not vertical.
	if(pstRay->dR >= 0.1)
	{
		// Call function NeqCalcRayProperties2 to calculate additional ray properties.
		NeqCalcRayProperties2(pstP1, pstP2, pstRay, pdSinSig, pdCosSig);
	}
	
	return bError;
}

static void NeqCalcRayProperties1(SPoint_st *pstP1, SPoint_st *pstP2, SPoint_st *pstRay, double *pdZeta)
{
	// This function calculates the properties of the ray. It does not calculate as many properties if the ray
	// is vertical, as they are not needed. ‘Ray’ is straight line passing through p1 and p2. pstRay->dLat,
	// pstRay->dLng and pstRay->dR are coordinates of ray perigee, i.e. point on ray closest to centre of Earth.

	// Check if ray is vertical.
	if(fabs(pstP2->dLat - pstP1->dLat) < 1e-5 && fabs(pstP2->dLng - pstP1->dLng) < 1e-5)
	{
		// Set the ray latitude and longitude to be the same as point 1.
		pstRay->dLat = pstP1->dLat;
		pstRay->dLng = pstP1->dLng;
		
		// Set the ray to have no slant.
		pstRay->dR = 0;
		*pdZeta = 0;
	}
	else
	{
		// Calculate and store sine and cosine of point 1 and point 2 latitudes.
		pstP1->dSinLat = sin(pstP1->dLat*DR);
		pstP1->dCosLat = cos(pstP1->dLat*DR);
		pstP2->dSinLat = sin(pstP2->dLat*DR);
		pstP2->dCosLat = cos(pstP2->dLat*DR);
		
		// Calculate temporary variables.
		double cosDl12 = cos((pstP2->dLng - pstP1->dLng)*DR);
		double sinDl12 = sin((pstP2->dLng - pstP1->dLng)*DR);
		double cosDel  = pstP1->dSinLat*pstP2->dSinLat + pstP1->dCosLat*pstP2->dCosLat*cosDl12;
		double sinDel  = sqrt(1 - NeqSquared(cosDel));
		
		// Calculate and store pdZeta (zenith angle of p2, seen from p1).
		*pdZeta = atan2(sinDel, cosDel - pstP1->dR/pstP2->dR);
		
		// Calculate temporary variables.
		double sinSigp = sinDl12*pstP2->dCosLat/sinDel;
		double cosSigp = (pstP2->dSinLat - cosDel*pstP1->dSinLat)/(sinDel*pstP1->dCosLat);
		double delp    = PI/2 - *pdZeta;
		double sinDelp = sin(delp);
		double cosDelp = cos(delp);
		
		// Calculate ray perigee latitude.
		double sinPhp  = pstP1->dSinLat*cosDelp - pstP1->dCosLat*sinDelp*cosSigp;
		double cosPhp  = sqrt(1 - NeqSquared(sinPhp));
		pstRay->dLat   = atan2(sinPhp, cosPhp)*RD;
		
		// Calculate ray perigee longitude.
		double sinLamp = -sinSigp*sinDelp/cosPhp;
		double cosLamp = (cosDelp - pstP1->dSinLat*sinPhp)/(pstP1->dCosLat*cosPhp);
		pstRay->dLng   = atan2(sinLamp, cosLamp)*RD + pstP1->dLng;
		
		// Calculate radius of ray perigee.
		pstRay->dR = pstP1->dR*sin(*pdZeta);
		pstRay->dS = 0;
		
		*pdZeta *= RD;
	}
}

static void NeqCalcRayProperties2(SPoint_st *pstP1, SPoint_st *pstP2, SPoint_st *pstRay, double *pdSinSig, double *pdCosSig)
{
	// This function calculates the sine and cosine of end point latitudes and azimuth. It is only called for slanted rays.

	// Calculate sine and cosine of end point latitudes (using ray perigee latitude for point 1).
	pstP1->dSinLat = sin(pstRay->dLat*DR);
	pstP1->dCosLat = cos(pstRay->dLat*DR);
	pstP2->dSinLat = sin(pstP2->dLat*DR);
	pstP2->dCosLat = cos(pstP2->dLat*DR);
	
	// Calculate difference in longitude of ray end points.
	double deltaLong = (pstP2->dLng - pstRay->dLng)*DR;
	
	// Check if latitude of lower end point is +-90 degrees - would cause divide by zero error later on.
	if(fabs(pstRay->dLat - 90) < 1e-10)
	{
		// Set sine and cosine of azimuth.
		*pdSinSig = 0;
		*pdCosSig = (pstRay->dLat > 0 ? -1 : +1);
	}
	else
	{
		// Calculate sine and cosine of angular distance between ends of ray (psi).
		double cosPsi = pstP1->dSinLat*pstP2->dSinLat + pstP1->dCosLat*pstP2->dCosLat*cos(deltaLong);
		double sinPsi = sqrt(1 - NeqSquared(cosPsi));
		
		// Calculate sine and cosine of azimuth.
		*pdSinSig = pstP2->dCosLat*sin(deltaLong)/sinPsi;
		*pdCosSig = (pstP2->dSinLat - pstP1->dSinLat*cosPsi)/(sinPsi*pstP1->dCosLat);
	}
}

static double NeqIntegrate(IntegrateData_st *pstIntegrateData, double dH1, double dH2, int siCurrentLevel, SPoint_st *pstPactual, double *pdNmax, LayerProperties_st *pstLayers)
{
	// Integration function for calculating TEC along rays using Kronrod G7-K15 adaptive quadrature method. This method involves
	// sampling values at 15 points and calculating the integration from them. At the same time it misses out half of the
	// points to see what difference it makes and therefore the likely error contained in the result, before deciding whether
	// to accept the result, or to split the portion into two and try again in order to improve accuracy.
	// Note that this method is recursive but has appropriate safeguards in the form of the recursion limit passed in from configuration.
	double dResult;
	
	// Set Kronrod integration coefficients (G7-K15).

	const double wi[15] = {0.022935322010529224963732008058970, 0.063092092629978553290700663189204, 0.104790010322250183839876322541518,
                           0.140653259715525918745189590510238, 0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
                           0.204432940075298892414161999234649, 0.209482141084727828012999174891714, 0.204432940075298892414161999234649,
                           0.190350578064785409913256402421014, 0.169004726639267902826583426598550, 0.140653259715525918745189590510238,
                           0.104790010322250183839876322541518, 0.063092092629978553290700663189204, 0.022935322010529224963732008058970};
						   
	const double wig[7] = {0.129484966168869693270611432679082, 0.279705391489276667901467771423780, 0.381830050505118944950369775488975,
                           0.417959183673469387755102040816327, 0.381830050505118944950369775488975, 0.279705391489276667901467771423780,
                           0.129484966168869693270611432679082};					   
						   
	const double xi[15] = {-0.991455371120812639206854697526329,-0.949107912342758524526189684047851,-0.864864423359769072789712788640926,
                           -0.741531185599394439863864773280788,-0.586087235467691130294144838258730,-0.405845151377397166906606412076961,
                           -0.207784955007898467600689403773245, 0,                                   0.207784955007898467600689403773245,
                            0.405845151377397166906606412076961, 0.586087235467691130294144838258730, 0.741531185599394439863864773280788,
                            0.864864423359769072789712788640926, 0.949107912342758524526189684047851, 0.991455371120812639206854697526329};

	// Calculate the midpoint, hh and the half difference h2.
	double h2 = (dH2 - dH1)/2;
	double hh = (dH2 + dH1)/2;
	
	// Initialize integration results intk and intg, and G7 counter Gind.
	double intk = 0;
	double intg = 0;
	int    Gind = 0;
	
	// Loop through the G15 and K7 integration points.
	int i;
	for(i = 0; i < 15; i++)
	{
		double x = h2*xi[i] + hh;
		double y;
		
		// If ray is vertical.
		if(pstIntegrateData->bVert == 1)
		{
			y = NeqGetNeOnVertRay(x, pstLayers, pdNmax);
		}
		else
		{
			y = NeqGetNeOnSlantRay(x, &pstIntegrateData->pstNeQuick, &pstIntegrateData->pstGeom, pdNmax, pstPactual, pstLayers, &pstIntegrateData->pstCurrCCIR);
		}
		
		// Accumulate on to the k15 total.
		intk = intk + y*wi[i];
		
		// If this is a G7 point.
		if(i % 2 == 1)
		{
			intg = intg + y*wig[Gind++];
		}
	}
	
	// Complete the calculation of the integration results.
	intk *= h2;
	intg *= h2;
	
	// Check if the result is within tolerance.
	if(fabs((intk - intg)/intk) <= pstIntegrateData->dTolerance || fabs(intk - intg) <= pstIntegrateData->dTolerance)
	{
		// Result is within tolerance so set return value equal to intk.
		dResult = intk;
	}
	else if(siCurrentLevel == pstIntegrateData->pstNeQuick.siMaxRecurse)
	{
		// Can do not further integration. Set return value intk as best guess.
		dResult = intk;
	}
	else
	{
		// Result is not within tolerance.
		// Split portion int two equal halves (from dH1 to dH1+h2 and from dH1+h2 to dH2 with
		// h2 = (dH2-dH1)/2 and call NeqIntegrate on each new portion.
		double dResult1 = NeqIntegrate(pstIntegrateData,    dH1, dH1+h2, siCurrentLevel+1, pstPactual, pdNmax, pstLayers);
		double dResult2 = NeqIntegrate(pstIntegrateData, dH1+h2,    dH2, siCurrentLevel+1, pstPactual, pdNmax, pstLayers);
		
		// Sum the return values from the two NeqIntegrate calls and set as return value.
		dResult = dResult1 + dResult2;
	}
	
	return dResult;
}

static double NeqGetNeOnVertRay(double dH, LayerProperties_st *pstLayers, double *pdNmax)
{
	// This function returns electron density at a specified point along a vertical ray.
	double dResult;
	
	// Check if specified height is above F2 peak.
	if(dH > pstLayers->pdPeakHeight[F2])
	{
		dResult = NeqCalcTopsideNe(dH, pstLayers, pdNmax);
	}
	else
	{
		dResult = NeqCalcBottomsideNe(dH, pstLayers);
	}
	return dResult;
}

static double NeqGetNeOnSlantRay(double dS, NeQuickInputData_st *pstNeQuickInputData, GeometryData_st *pstGeom, double *pdNmax, SPoint_st *pstPactual, LayerProperties_st *pstLayers, CurrentCCIR_st *pstCurrCCIR)
{
	// This function returns electron density at the specified point along a slanted ray.
	double dResult;
	
	// Call function NeqCalcLLHOnRay to adjust position information for current position along ray, pstPactual.
	NeqCalcLLHOnRay(dS, &pstGeom->stRay, &pstGeom->stP1, pstGeom->dSinSig, pstGeom->dCosSig, pstPactual);
	
	// Call function NeqCalcEpstParams to recalculate ionosphere information now that the latitude and longitude have changed.
	NeqCalcEpstParams(pstNeQuickInputData, pstPactual, pstGeom->dSinDelta, pstGeom->dCosDelta, pdNmax, pstLayers, pstCurrCCIR);
	
	if(pstPactual->dH > pstLayers->pdPeakHeight[F2])
	{
		dResult = NeqCalcTopsideNe(pstPactual->dH, pstLayers, pdNmax);
	}
	else
	{
		dResult = NeqCalcBottomsideNe(pstPactual->dH, pstLayers);
	}
	return dResult;
}

static void NeqCalcLLHOnRay(double dS, SPoint_st *pstRay, SPoint_st *pstP1, double dSinSig, double dCosSig, SPoint_st *pstPactual)
{
	// This function sets the latitude, longitude and height of the current position along the ray stPactual according
	// to the specified slant position s (the distance along the slanted ray).

	// Calculate trig of angle at center of earth between lines to ray perigee and point on ray (Del).
	double tanDel = dS/pstRay->dR;
	double cosDel = 1/sqrt(1 + NeqSquared(tanDel));
	double sinDel = tanDel*cosDel;
	
	// Calculate latitude.
	double arg = pstP1->dSinLat*cosDel + pstP1->dCosLat*sinDel*dCosSig;
	pstPactual->dLat = atan2(arg, sqrt(1 - NeqSquared(arg)))*RD;
	
	// Calculate longitude.
	double cLong = atan2(sinDel*dSinSig*pstP1->dCosLat, cosDel - pstP1->dSinLat*arg)*RD;
	pstPactual->dLng = cLong + pstRay->dLng;
	
	// Calculate height.
	pstPactual->dH = sqrt(NeqSquared(dS) + NeqSquared(pstRay->dR)) - RE;
	pstPactual->dR = sqrt(NeqSquared(dS) + NeqSquared(pstRay->dR));
	pstPactual->dS = dS;
}

static void NeqCalcEpstParams(NeQuickInputData_st *pstNeQuickInputData, SPoint_st *pstPactual, double dSinDelta, double dCosDelta, double *pdNmax, LayerProperties_st *pstLayers, CurrentCCIR_st *pstCurrCCIR)
{
	// This function calculates the values of ionospheric properties for the current latitude, longitude, time, etc.

	// Calculate MODIP at current longitude and latitude.
	double modip = NeqCalcModip(pstPactual->dLat, pstPactual->dLng, &pstNeQuickInputData->pstModip);
	
	// Retrieve F10.7 flux at user receiver position.
	double Az = pstNeQuickInputData->dAzBase;

	// Calculate R12 sunspot number from solar flux (Az).
	double R12 = sqrt(1123.6*(Az - 63.7) + 167273) - 408.99;
	
	// Initialize pdNmax to -1.0 to force NeqCalcTopsideNe function to call NeqCalcBottomsideNe.
	*pdNmax = -1.0;
	
	// Check if month or solar flux has changed since spherical Legendre coefficients were last computed.
	if(pstCurrCCIR->siMonth != pstNeQuickInputData->siMonth || pstCurrCCIR->dR12 != R12)
	{
		// Load working matrices with CCIR coefficients.
		double RR2 = R12/100;
		double RR1 = 1 - RR2;
		
		// Blend high and low activity cases in ratio RR2:RR1
		int i, j;
		for(i = 0; i < 13; i++)
		{
			for(j = 0; j < 76; j++)
			{
				pstCurrCCIR->pdF0F2[j*13 + i] = pstNeQuickInputData->pstCCIR.pdF2[pstNeQuickInputData->siMonth-1][i][j][0]*RR1 +
												pstNeQuickInputData->pstCCIR.pdF2[pstNeQuickInputData->siMonth-1][i][j][1]*RR2;
			}
		}
		for(i = 0; i < 9; i++)
		{
			for(j = 0; j < 49; j++)
			{
				pstCurrCCIR->pdM3000F2[j*9 + i] = pstNeQuickInputData->pstCCIR.pdM3000[pstNeQuickInputData->siMonth-1][i][j][0]*RR1 +
												  pstNeQuickInputData->pstCCIR.pdM3000[pstNeQuickInputData->siMonth-1][i][j][1]*RR2;
			}
		}
		
		// Set current R12, UT and month.
		pstCurrCCIR->dR12    = R12;
		pstCurrCCIR->dUT     = pstNeQuickInputData->dUT;
		pstCurrCCIR->siMonth = pstNeQuickInputData->siMonth;
		
		// Call NeqCalcSphLegCoeffs with current time and blended CCIR information
		NeqCalcSphLegCoeffs(pstCurrCCIR->dUT, pstCurrCCIR);
	}
	else if(pstCurrCCIR->dUT != pstNeQuickInputData->dUT)
	{
		// Set current UT.
		pstCurrCCIR->dUT = pstNeQuickInputData->dUT;

		// Call NeqCalcSphLegCoeffs with current time and blended CCIR information
		NeqCalcSphLegCoeffs(pstCurrCCIR->dUT, pstCurrCCIR);
	}
	
	// Calculate sin^n(modip) array for n = 0 to 11, and cosine of latitude.
	pstPactual->dSinLat = sin(pstPactual->dLat*DR);
	pstPactual->dCosLat = cos(pstPactual->dLat*DR);
	
	double pdSinModipToN[12];
	pdSinModipToN[0] = 1;
	
	int n;
	for(n = 1; n < 12; n++)
	{
		pdSinModipToN[n] = pdSinModipToN[n-1]*sin(modip*DR);
	}
	
	// Get F2 layer data from CCIR.
	//
	// For f0F2, call NeqGetF2FreqsFromCCIR function with foF2 working arrays.
	double f0F2  = NeqGetF2FreqFromCCIR(pstPactual->dCosLat, pstPactual->dLng, pstCurrCCIR->pdLegCoeffs_F0, pdSinModipToN, 0);
	//
	// For M3000, call NeqGetF2FreqsFromCCIR function with M3000 working arrays.
	double M3000 = NeqGetF2FreqFromCCIR(pstPactual->dCosLat, pstPactual->dLng, pstCurrCCIR->pdLegCoeffs_M3000, pdSinModipToN, 1);
	
	// Calculate local time and map back into 24 hour period.
	double xlt = pstCurrCCIR->dUT + pstPactual->dLng/15.0;
	if(xlt < 0)
	{
		xlt = xlt + 24.0;
	}
	else if(xlt >= 24)
	{
		xlt = xlt - 24.0;
	}
	
	// Calculate Solar Zenith angle.
	double cosChi = sin(pstPactual->dLat*DR)*dSinDelta + cos(pstPactual->dLat*DR)*dCosDelta*cos(PI*(12 - xlt)/12);
	double chi    = atan2(sqrt(1 - NeqSquared(cosChi)), cosChi)*RD;
	double chi0   = 86.23;
	
	// Set season flag.
	double seas = 0;
	switch(pstCurrCCIR->siMonth)
	{
	case 1:
	case 2:
	case 11:
	case 12:
		seas = -1;
	break;
	
	case 5:
	case 6:
	case 7:
	case 8:
		seas = +1;
	break;
	
	default:
		seas = 0;
	break;
	}
	
	// Estimate foE and f0F1.
	// The model for foE adopted is based on the solar zenith angle law. An exponential transition
	// between "day" and "night" is used which ensures differentiability.
	// For daytime the model takes foF1 = 1.4foE, for night-time foF1 = 0, using the same exponential
	// day-night transition as for foE. In addition foF1 is reduced by 15% if too close to foF2.
	double ee = NeqClipExp(0.3*pstPactual->dLat);
	
	seas = seas*(ee - 1)/(ee + 1);
	
	double chin = NeqJoin(90.0 - 0.24*NeqClipExp(20.0 - 0.2*chi), chi, 12, chi - chi0);
	double sfac = (1.112 - 0.019*seas)*sqrt(sqrt(Az));
	double fa   = sfac*NeqClipExp(log(cos(chin*DR))*0.3);
	
	// Calculate E peak plasma frequency.
	double f0E = sqrt(fa*fa + 0.49);
	
	// Calculate F1 peak plasma frequency and set to zero if negligible.
	double f0F1;
	f0F1 = NeqJoin(1.4*f0E, 0, 1000.0, f0E - 2); // Titheridge's formula f0F1 = 1.4f0F2.
	f0F1 = NeqJoin(0, f0F1, 1000.0, f0E - f0F1);
	f0F1 = NeqJoin(f0F1, 0.85*f0F1, 60.0, 0.85*f0F2 - f0F1);
	
	if(f0F1 < 1e-6)
	{
		f0F1 = 0;
	}
	
	// Calculate peak electron densities from critical frequencies for the F2, F1 and E layers.
	double Nm[3], Amp[3];
	
	// F2, F1 and E correspond to indices 0, 1 and 2 respectively.
	Nm[F2] = NeqCriticalFreqToNe(f0F2);
	Nm[F1] = NeqCriticalFreqToNe(f0F1);
	Nm[E]  = NeqCriticalFreqToNe(f0E);
	
	// Calculate height of electron density peaks for the layers. F2 peak is calculated each time, 
	// E layer peak is fixed at 120km and the F1 peak is set halfway between them.
	double PeakHeight[3];
	
	PeakHeight[F2] = NeqCalcF2PeakHeight(M3000, f0E, f0F2);
	PeakHeight[E ] = 120;
	PeakHeight[F1] = (PeakHeight[F2] + PeakHeight[E])/2;
	
	// Calculate density gradient at base of F2 layer (10^9 m^-3 km^-1).
	double NdHmx = 0.01*exp(-3.467 + 0.857*log(f0F2*f0F2) + 2.02*log(M3000));
	
	// Calculate bottom-side thickness parameters.
	double BotThick[3], TopThick[3];
	
	BotThick[F2] = 0.385*Nm[F2]/NdHmx;
	TopThick[F1] = 0.3*(PeakHeight[F2] - PeakHeight[F1]);
	BotThick[F1] = 0.5*(PeakHeight[F1] - PeakHeight[E ]);
	TopThick[E ] = BotThick[F1];
	if(TopThick[E] < 7)
	{
		TopThick[E] = 7;
	}
	BotThick[E] = 5;
	
	Amp[F2] = 4*Nm[F2];
	Amp[F1] = 4*Nm[F1];
	Amp[E ] = 4*Nm[E ];

	// Calculate Epstein function amplitudes.
	// The construction of the vertical profile is based on "anchor" points related to the ionospheric
	// characteristics of the main layers routinely scaled from the ionograms: f0F2, M(3000)F2, f0F1 and f0E.
	if(f0F1 < 0.5)
	{
		Amp[F1] = 0;
		Amp[E ] = 4*(Nm[E] - NeqEpstein(Amp[F2], PeakHeight[F2], BotThick[F2], PeakHeight[E]));
	}
	else
	{
		int i;
		for(i = 0; i < 5; i++)
		{
			Amp[F1] = 4*(Nm[F1] - NeqEpstein(Amp[F2], PeakHeight[F2], BotThick[F2], PeakHeight[F1])
							    - NeqEpstein(Amp[E ], PeakHeight[E ], TopThick[E ], PeakHeight[F1]));

			Amp[F1] = NeqJoin(Amp[F1], 0.8*Nm[F1], 1, Amp[F1] - 0.8*Nm[F1]);

			Amp[E ] = 4*(Nm[E ] - NeqEpstein(Amp[F1], PeakHeight[F1], BotThick[F1], PeakHeight[E])
							    - NeqEpstein(Amp[F2], PeakHeight[F2], BotThick[F2], PeakHeight[E]));
		}
	}
	Amp[E] = NeqJoin(Amp[E], 0.05, 60.0, Amp[E] - 0.005);
	
	// Calculate shape factor for topside F2 region.
	double b2k;
	if(pstCurrCCIR->siMonth > 3 && pstCurrCCIR->siMonth < 10) // April to September.
	{
		b2k = 6.705 - 0.014*R12 - 0.008*PeakHeight[F2];
	}
	else	// October to May.
	{
		b2k = -7.77 + 0.097*NeqSquared(PeakHeight[F2]/BotThick[F2]) + 0.153*Nm[F2];
	}
	b2k = NeqJoin(b2k, 2, 1, b2k - 2.0);
	b2k = NeqJoin(8, b2k, 1, b2k - 8.0);
	
	// Adjust the vertical TEC value to take into account exosphere electron density.
	TopThick[F2] = b2k*BotThick[F2];
	
	double x = (TopThick[F2] - 150.0)/100.0;
	double v = (0.041163*x - 0.183981)*x + 1.424472;
	
	TopThick[F2] = TopThick[F2]/v;
		
	// Store frequencies, M3000 and the other parameters in LayerProperties_st data structure.
	pstLayers->pdF0[F2] = f0F2;
	pstLayers->pdF0[F1] = f0F1;
	pstLayers->pdF0[E ] = f0E;
	pstLayers->dM3000   = M3000;
	
	pstLayers->pdAmp[F2] = Amp[F2];
	pstLayers->pdAmp[F1] = Amp[F1];
	pstLayers->pdAmp[E]  = Amp[E ];
	
	pstLayers->pdBotThick[F2] = BotThick[F2];
	pstLayers->pdBotThick[F1] = BotThick[F1];
	pstLayers->pdBotThick[E ] = BotThick[E ];
	
	pstLayers->pdTopThick[F2] = TopThick[F2];
	pstLayers->pdTopThick[F1] = TopThick[F1];
	pstLayers->pdTopThick[E ] = TopThick[E ];
	
	pstLayers->pdPeakHeight[F2] = PeakHeight[F2];
	pstLayers->pdPeakHeight[F1] = PeakHeight[F1];
	pstLayers->pdPeakHeight[E ] = PeakHeight[E ];
}

static void NeqCalcSphLegCoeffs(double dUT, CurrentCCIR_st *pstCurrCCIR)
{
	// This function calculates the spherical Legendre coefficients which are used in the
	// calculations for foF2 or M(3000)F2 frequencies calculated from CCIR map file data.

	const int numMaxHarm = 6;
	
	double t = (dUT*15 - 180)*DR;
	
	double sinHarm[6], cosHarm[6];
	
	sinHarm[0] = sin(t);
	cosHarm[0] = cos(t);
	
	int i, k, N;
	for(k = 1; k < numMaxHarm; k++)
	{
		sinHarm[k] = sinHarm[k-1]*cosHarm[0] + cosHarm[k-1]*sinHarm[0];
		cosHarm[k] = cosHarm[k-1]*cosHarm[0] - sinHarm[k-1]*sinHarm[0];
	}
	
	for(i = 0, N = 13; i < 76; i++)
	{
		pstCurrCCIR->pdLegCoeffs_F0[i] = pstCurrCCIR->pdF0F2[i*N];
		for(k = 0; k < 6; k++)
		{
			pstCurrCCIR->pdLegCoeffs_F0[i] += (pstCurrCCIR->pdF0F2[i*N+2*k+1]*sinHarm[k] + pstCurrCCIR->pdF0F2[i*N+2*k+2]*cosHarm[k]);
		}
	}
	
	for(i = 0, N = 9; i < 49; i++)
	{
		pstCurrCCIR->pdLegCoeffs_M3000[i] = pstCurrCCIR->pdM3000F2[i*N];
		for(k = 0; k < 4; k++)
		{
			pstCurrCCIR->pdLegCoeffs_M3000[i] += (pstCurrCCIR->pdM3000F2[i*N+2*k+1]*sinHarm[k] + pstCurrCCIR->pdM3000F2[i*N+2*k+2]*cosHarm[k]);
		}
	}
}

static double NeqGetF2FreqFromCCIR(double dCosLat, double dLng, double *pdLegCoeffs, double pdSinModipToN[12], int siMode)
{
	// This function returns foF2 or M(3000)F2 calculated from CCIR map file data.
	double dResult = 0;
	
	int QF[9] = {11,11,8,4,1,0,0,0,0};
	int QM[7] = {6,7,5,2,1,0,0};
	
	int k1, *nq;
	
	// Check mode: siMode=0 compute f0F2 value at current point. siMode compute M3000 value at current point.
	if(siMode == 0)
	{
		// Set constants used in spherical Legendre expansion.
		nq = QF;
		k1 = 9;
	}
	else
	{
		// Set constants used in spherical Legendre expansion.
		nq = QM;
		k1 = 7;
	}
	
	// Compute output value as sum of spherical Legendre functions.
	double pdSinModipToN_corr[12];
	
	int i, j, k, l, R;
	for(j = 0; j < 12; j++)
	{
		pdSinModipToN_corr[j] = (fabs(pdSinModipToN[j]) > 1e-30 ? pdSinModipToN[j] : 0);
	}
	
	for(j = 0; j <= nq[0]; j++)
	{
		dResult += pdLegCoeffs[j]*pdSinModipToN_corr[j];
	}
	
	for(k = 1; k <= (k1 - 1); k++)
	{
		for(i = 0; i <= nq[k]; i++)
		{
			R = 0;
			for(l = 0; l <= (k - 1); l++)
			{
				R += (nq[l] + 1);
			}
			R = 2*R - nq[0] - 1;
			
			dResult += pdSinModipToN_corr[i]*pow(dCosLat,k)*(pdLegCoeffs[R+2*i]*cos(k*dLng*DR) + pdLegCoeffs[R+2*i+1]*sin(k*dLng*DR));
		}
	}
	
	return dResult;
}

static double NeqCriticalFreqToNe(double dF0)
{
	// From critical frequency, calculates the associated electron density.
	return 0.124*dF0*dF0;
}

static double NeqCalcF2PeakHeight(double dM3000, double dF0E, double dF0F2)
{
	// This function calculates F2 layer peak height hm[F2] from foE, foF2 and M3000. It
	// is based on the method of Dudeney(1983), but modified such that the ratio foF2/foE
	// is clipped at 1.75 using NeqJoin. Note that the clipping is ‘soft’, the 1st
	// derivative is continuous but note the clipped value can be slightly below 1.75 at
	// the join (but note analysis indicate >1.73). Also, Dudeney uses a figure of 1470
	// rather than 1490 in the numerator of hmF2 and a figure of 1.296 rather than
	// 1.2967 in the denominator of MF.

	// Height at which electron density peaks in F2 layer.
	double PeakHeightF2, dAM;
	
	double MF = dM3000*sqrt((0.0196*dM3000*dM3000 + 1)/(1.2967*dM3000*dM3000 - 1));
	
	if(dF0E >= 1e-30)
	{
		double r  = dF0F2/dF0E;
		double r2 = (r*exp(20*(r - 1.75)) + 1.75)/(exp(20*(r - 1.75)) + 1);
		
		dAM = 0.253/(r2 - 1.215) - 0.012;
	}
	else
	{
		// Limit as r becomes very large.
		dAM = -0.012;
	}
	// Compute peak height in F2 (in km) later as:
	PeakHeightF2 = 1490*MF/(dM3000 + dAM) - 176;
	
	return PeakHeightF2;
}

static double NeqEpstein(double dNmax, double dHmax, double dB, double dH)
{
	// Evaluates Epstein layer function (without the normalization factor of 4).
	// This function only calculates a 1/4 of its true value because the factor of 4 is 
	// completed once other mathematical formulas have been done to decrease computation time.
	double dResult;
	
	double expTerm = NeqClipExp((dH - dHmax)/dB);
	
	dResult = dNmax*expTerm/NeqSquared(1 + expTerm);
	
	return dResult;
}

static double NeqCalcTopsideNe(double dH, LayerProperties_st *pstLayers, double *pdNmax)
{
	// This function calculates electron content at the specified height, in the top part of
	// the ionosphere above the F2 electron density peak point.
	double dResult;
	
	const double g    = 0.125;
	const double rfac = 100;
	
	double dH_PeakHF2 = dH - pstLayers->pdPeakHeight[F2];
	double TopThickF2 = pstLayers->pdTopThick[F2];
	
	double ee = NeqClipExp(dH_PeakHF2/(TopThickF2*(1 + rfac*g*dH_PeakHF2/(rfac*TopThickF2 + g*dH_PeakHF2))));
	double ep;
	
	if(ee > 1e11)
	{
		// Deal with limit when ee is very large.
		ep = 4/ee;
	}
	else
	{
		ep = 4*ee/NeqSquared(1 + ee);
	}
	
	// If pdNmax has not been calculated for this location call function NeqCalcBottomsideNe
	// to calculate pdNmax at PeakHeight[F2] (crossover point).
	if(*pdNmax < 0)
	{
		*pdNmax = NeqCalcBottomsideNe(pstLayers->pdPeakHeight[F2], pstLayers);
	}
	
	dResult = (*pdNmax)*ep;
	
	return dResult;
}

static double NeqCalcBottomsideNe(double dHH, LayerProperties_st *pstLayers)
{
	// This function calculates electron content at the specified height dHH, in the bottom
	// part of the ionosphere below the F2 peak height. The function sums semi-Epstein Layers
	// with a modification to reduce excessive Ne around F2 peak and 100km. 
	double dResult;
		
	const double h0 = 100;
	const double f1 = 10;
	const double f2 = 1;
		
	double B[3];
	
	B[F2]= pstLayers->pdBotThick[F2];
	
	if(dHH > pstLayers->pdPeakHeight[F1])
	{
		B[F1] = pstLayers->pdTopThick[F1];
	}
	else
	{
		B[F1] = pstLayers->pdBotThick[F1];
	}
	
	if(dHH > pstLayers->pdPeakHeight[E])
	{
		B[E] = pstLayers->pdTopThick[E];
	}
	else
	{
		B[E] = pstLayers->pdBotThick[E];
	}

	int j;
	double arg, s[3], ds[3];
	if(dHH < 100) // If height is below 100 km.
	{
		for(j = 0; j < 3; j++)	// Iterate for each ionospheric section, F2, F1, E
		{
			if(j == F2)
			{
				arg = (h0 - pstLayers->pdPeakHeight[F2])/B[F2];
			}
			else
			{
				arg = ((h0 - pstLayers->pdPeakHeight[j])/B[j])*exp(f1/(1 + f2*fabs(h0 - pstLayers->pdPeakHeight[F2])));
			}
			
			if(fabs(arg) > 25)
			{
				s[j] = ds[j] = 0;
			}
			else
			{
				s [j] = pstLayers->pdAmp[j]*exp(arg)/NeqSquared(1 + exp(arg));
				ds[j] = (1 - exp(arg))/(B[j]*(1 + exp(arg)));
			}
		}
		double Hd  = 10;
		double z   = (dHH - h0)/Hd;
		double aN0 = 1e11*(s[0] + s[1] + s[2]);
		
		double bf  = 1;
		if(s[0] + s[1] + s[2] != 0)
		{
			bf = 1 - Hd*(s[0]*ds[0] + s[1]*ds[1] + s[2]*ds[2])/(s[0] + s[1] + s[2]);
		}
		dResult = aN0*NeqClipExp(1 - bf*z - NeqClipExp(-z));
	}
	else
	{
		for(j = 0; j < 3; j++)
		{
			if(j == F2)
			{
				arg = (dHH - pstLayers->pdPeakHeight[F2])/B[F2];
			}
			else
			{
				arg = ((dHH - pstLayers->pdPeakHeight[j])/B[j])*exp(f1/(1 + f2*fabs(dHH - pstLayers->pdPeakHeight[F2])));
			}
			
			if(fabs(arg) > 25)
			{
				s[j] = 0;
			}
			else
			{
				s[j] = pstLayers->pdAmp[j]*exp(arg)/NeqSquared(1 + exp(arg));
			}
		}
		dResult = 1e11*(s[0] + s[1] + s[2]);
	}
	
	// Return the computed electron density at the given height in the bottom part of the
	// ionosphere below the F2 electron density peak.
	return dResult;
}

static double NeqJoin(double dF1, double dF2, double dAlpha, double dX)
{
	// Allows smooth joining of functions f1 and f2 (i.e. continuous first derivatives) at
	// origin. Alpha determines width of transition region. Calculates value of joined functions at x.
	double dResult;

	double ee;
	ee = NeqClipExp(dAlpha*dX);

	dResult = (dF1*ee + dF2)/(ee + 1);

	return dResult;
}

static double NeqClipExp(double dPower)
{
	// A clipped exponential function - always return valid output.
	double dResult;
	
	if(dPower > 80)
	{
		dResult = 5.5406e34;
	}
	else if(dPower < -80)
	{
		dResult = 1.8049e-35;
	}
	else
	{
		dResult = exp(dPower);
	}
	
	// Return clipped exponential value.
	return dResult;
}

static double NeqSquared(double dValue)
{
	// This function calculates the square of a number.
	double dResult;
	
	dResult = dValue*dValue;
	
	return dResult;
}

// Read CCIR files and store it in the table pstCCIR, according to the format specified in the NeQuick G
// ICD. It shall be invoked one time to initialize the input structure for the NeQuick G algorithm.
// Returns a flag error equal to 1 if any of the CCIR files has not been found, 0 if everything is OK.
// Although a char string including the path where to find the CCIR files, their name is assumed to be
// of the form ccirXX.asc, being XX the month number plus 10 (two digits).
NeQuick_status NeqReadCCIRFiles(char *filepath, CCIR_st *pstCCIR)
{
	char filename[256];

	FILE *fp;

	int j, k, l, month;
	for(month = 1; month <= 12; month++)
	{
		// Get path where CCIR file for current month is located.
		sprintf(filename, "%s/ccir%02d.asc", filepath, month+10);

		if((fp = fopen(filename, "rt")) != NULL)
		{
			// Read F2 and M3000 data structures from CCIR file according to ICD.
			for(l = 0; l < 2; l++)
			{
				for(k = 0; k < 76; k++)
				{
					for(j = 0; j < 13; j++)
					{
						fscanf(fp, "%lf", &pstCCIR->pdF2[month-1][j][k][l]);
					}
				}
			}

			for(l = 0; l < 2; l++)
			{
				for(k = 0; k < 49; k++)
				{
					for(j = 0; j < 9; j++)
					{
						fscanf(fp, "%lf", &pstCCIR->pdM3000[month-1][j][k][l]);
					}
				}
			}
			fclose(fp);
		}
		else
		{
			// Return error. CCIR files not properly read.
			return E_ERROR;
		}
	}
	return E_OK;
}

// Read MODIP files and store it in the table pstModip, according to the format specified in the NeQuick G
// ICD. It shall be invoked one time to initialize the input structure for the NeQuick G algorithm.
// Returns a flag error equal to 1 if the MODIP file has not been found, 0 if everything is OK.
// Although a char string including the path where to find the MODIP file, the name of this file is assumed
// to be modipNeQG_wrapped.asc.
NeQuick_status NeqReadMODIPFiles(char *filepath, MODIP_st *pstModip)
{
	NeQuick_status result = E_OK;

	char filename[256];

	FILE *fp;

	// Get path where MODIP file is located.
	sprintf(filename, "%s/modipNeQG_wrapped.asc", filepath);

	if((fp = fopen(filename, "rt")) != NULL)
	{
		// Read modified dip latitudes according to ICD. Notice that the first and last row and
		// columns are added as padding to ease the handling of this data structure.
		int j, k;
		for(j = 0; j < 39; j++)
		{
			for(k = 0; k < 39; k++)
			{
				fscanf(fp, "%lf", &pstModip->pdModip[j][k]);
			}
		}
		fclose(fp);
	}
	else
	{
		// Set error value. MODIP file not properly read.
		result = E_ERROR;
	}
	return result;
}
