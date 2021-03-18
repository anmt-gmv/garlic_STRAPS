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

#ifndef PI
#define PI 		(3.1415926535898)
#endif
#define DR 		(PI/180)
#define RD		(180/PI)
#define RE		(6371.2)
#define F2		0
#define F1		1
#define E		2

typedef enum
{
	E_OK    = 0,
	E_ERROR = 1
	
} NeQuick_status;

typedef struct
{
	double pdModip[39][39]; 		// Grid of modified dip latitude values. The grid is wrapped around the poles
} MODIP_st;

typedef struct
{
	double pdF2   [12][13][76][2];	// CCIR coefficients for computing FoF2, the critical frequency of the F2 layer
	double pdM3000[12][ 9][49][2];	// CCIR coefficients for computing M(3000)F2. The ratio of the maximum usable
									// frequency at a distance of 3000km to the F2 layer critical frequency foF2
} CCIR_st;

typedef struct
{
	int      siMaxRecurse;			// Maximum level of recursion allowed in Kronrod integration
	int      siMonth;				// Month during which STEC values is required
	int      siNumCoeff;			// Number of Az coefficients
	double   pdKronrodTol[2];		// Tolerances for Kronrod integration
	double   pdGssPosLLH[3];		// Receiver position (geodetics)
	double   pdSatPosLLH[3];		// Satellite position (geodetics)
	double   dUT;					// Time (UTC) at which STEC value is required
	double   pdCoeff[3];			// Az coefficients
	double   dAzBase;				// Az value at receiver locations
	MODIP_st pstModip;				// Structure containing grid of modified dip latitude values
	CCIR_st  pstCCIR;				// Structure containing CCIR coefficients for computing FoF2 and M(3000)F2
	
} NeQuickInputData_st;

typedef struct
{
	double dLat;					// Latitude of Point
	double dLng;					// Longitude of Point
	double dH;						// Height of Point
	double dR;						// Radius of Point
	double dS;						// Distance of Point to Ray Perigee
	double dSinLat;					// Sine of latitude of Point
	double dCosLat;					// Cosine of latitude of Point
	
} SPoint_st;

typedef struct
{
	int    siMonth;					// Month during which current STEC value has been computed
	double dLat;					// Latitude of Point
	double dR12;					// Current R12 index - twelve-month smoothed relative sunspot number
	double pdF0F2[988];				// Interpolated coefficients for computing FoF2 for current month and R12 conditions
	double pdM3000F2[441];			// Interpolated coefficients for computing M3000F2 for current month and R12 conditions
	double dUT;						// Time (UTC) at which current STEC value has been computed
	double pdLegCoeffs_F0[76];		// Spherical Legendre coefficients for calculating F0F2 for current month and R12 conditions
	double pdLegCoeffs_M3000[49];	// Spherical Legendre coefficients for calculating M(3000)F2 for current month and R12 conditions

} CurrentCCIR_st;

typedef struct
{
	double pdAmp[3];				// Epstein amplitude parameter
	double pdPeakHeight[3];			// Epstein peak height parameter
	double pdBotThick[3];			// Epstein bottom half-layer thickness parameter
	double pdTopThick[3];			// Epstein top half-layer thickness parameter
	double dM3000;					// Current M(3000)F2 value
	double pdF0[3];					// Current F0 (peak plasma frequency) for the F2, F1 and E layers
	
} LayerProperties_st;

typedef struct
{
	SPoint_st stP1;					// Information for point 1 (receiver)
	SPoint_st stP2;					// Information for point 2 (satellite)
	SPoint_st stRay;				// Information for ray
	SPoint_st stPactual;			// Information for current integration point
	double dZeta;					// Zenith angle of point 2 seen from point 1
	double dSinDelta;				// Sine of angle of declination of sun
	double dCosDelta;				// Cosine of angle of declination of sun
	double dSinSig;					// Sine of ray azimuth
	double dCosSig;					// Cosine of ray azimuth
	
} GeometryData_st;

typedef struct
{
	char 				bVert;		// Flag indicating whether ray is vertical or not
	NeQuickInputData_st pstNeQuick; // Input data to NeQuick Function
	GeometryData_st		pstGeom;	// Geometry data for ray
	CurrentCCIR_st		pstCurrCCIR;// foF2 and M(3000)F2 information for current month and R12	
	double				dTolerance;	// Tolerance for Kronrod integration
	
	double				dModipRx;   // Modified dip latitude at receiver point. Used as output
                                    // to calculate the geo-magnetic latitude and UIRE variance.
} IntegrateData_st;

// NeQuick G algorithm. Calculates the Slant Total Electron Content (STEC) taking as
// input data the user and satellite positions, plus broadcast ionosphere parameters.
// Return an element of the NeQuick_status enumeration (OK or ERROR).
NeQuick_status NeQuick(IntegrateData_st *pstIntegrateData, double *pdSTEC);

// Read CCIR files and store it in the table pstCCIR, according to the format specified in the NeQuick G
// ICD. It shall be invoked one time to initialize the input structure for the NeQuick G algorithm.
// Returns a flag error equal to 1 if any of the CCIR files has not been found, 0 if everything is OK.
// Although a char string including the path where to find the CCIR files, their name is assumed to be
// of the form ccirXX.asc, being XX the month number plus 10 (two digits).
NeQuick_status NeqReadCCIRFiles(char *filepath, CCIR_st *pstCCIR);

// Read MODIP files and store it in the table pstModip, according to the format specified in the NeQuick G
// ICD. It shall be invoked one time to initialize the input structure for the NeQuick G algorithm.
// Returns a flag error equal to 1 if the MODIP file has not been found, 0 if everything is OK.
// Although a char string including the path where to find the MODIP file, the name of this file is assumed
// to be modipNeQG_wrapped.asc.
NeQuick_status NeqReadMODIPFiles(char *filepath, MODIP_st *pstModip);
