#ifndef _DBTYPEDEF_HPP
#define _DBTYPEDEF_HPP

#include "dbdefine.h"
#include "GICSRx_defines.h"

typedef enum
{
	SYS_GPS    = 0,
	SYS_GAL    = 1,
	SYS_GLO    = 2,
	SYS_BEI    = 3
} sysFlag_t;

typedef enum
{
	FREQ_L1E1  = 0,
	FREQ_L5    = 1,
	FREQ_E5a   = 2,
	FREQ_E5b   = 3,
	FREQ_E6    = 4
} freqFlag_t;

typedef enum
{
	SIG_GPSL1   = 0,
	SIG_GALE1AD = 1,
	SIG_GALE1AP = 2,
	SIG_GALE1B  = 3,
	SIG_GALE1C  = 4,
	SIG_GPSL5I  = 5,
	SIG_GPSL5Q  = 6,
	SIG_GALE5aI = 7,
	SIG_GALE5aQ = 8,
	SIG_GALE5bI = 9,
	SIG_GALE5bQ = 10,
	SIG_GALE6AD = 11,
	SIG_GALE6AP = 12,
	SIG_GLOL1   = 13,
	SIG_BEIB1   = 14
} sigFlag_t;

//////////////////////////////////////////////////////////////////////////////
// External data

typedef struct
{
	double pos[3];
	double vel[3];
	double clk_bias;
	double clk_drift;
	double clk_driftrate;

} ephemsat_t;

//////////////////////////////////////////////////////////////////////////////
// GPS raw data

typedef enum
{
	OBS_VALID                  = 0,
	OBS_OK_BUT_UNUSED          = 1,
	REJ_INVALID_PRN            = 2,
	REJ_SV_NOT_TRACKED         = 3,
	REJ_SV_NOT_LOCKED          = 4,
	REJ_OUT_OF_RANGE           = 5,
	REJ_SV_STATE_NOT_AVAILABLE = 6,
	REJ_EPHEM_NOT_LOCKED       = 7,
	REJ_LOW_CN0                = 8,
	REJ_LOW_ELEV               = 9,
	REJ_CHISQUARE_TEST         = 10,
	REJ_FILTER_RESET           = 11,
	REJ_FILTER_INIT            = 12

} Obs_Status;

typedef struct
{
	sigFlag_t sigFlag;

	int    channel;
	char   lock;        // Available measurement at current epoch

	int    PRN;         // Satellite PRN.

	double C1;          // C/A code measurement (m)
	double L1;	    	// Carrier phase measurement (cycles)
	double D1;          // L1 Doppler measurement (m/s)
	float  S1;          // Carrier to noise density power ratio (dBHz)
	double lambda;		// Signal wavelength (m)
	float  elev;        // Elevation angle (deg)
	float  azim;        // Azimuth angle   (deg)

	char  iono_free;
	float atmsphCorr;  	// Fastrax atmospheric correction.
	char  atmsph_av;	// Atmospheric correction available.

	Obs_Status prange_status;
	Obs_Status cphase_status;
	Obs_Status doppler_status;

	float prange_residual;
	float doppler_residual;

	float prange_sigma;
	float doppler_sigma;

	ephemsat_t sat; // Satellite position (cartesian WGS84) and clock for ISOTROPY solution.

} obsdata_t;

typedef enum
{
	LNAV = 0,
	CNAV = 1

} gps_service;

typedef struct
{
	gps_service srv;

	char   vflg;        // Indicates if the ephemeris is available.
	int    prn;		    // PRN of the GPS SV.
	int    nWeekNumber; // Week Number.
	int    nIODE;       // Issue Of Data Ephemeris.
	int    nIODC;       // Issue Of Data Clock.
	int    nHealth;		// Health flag.
    char   fitInt;
	int    CodeOnL2;
	char   L2PFlag;
	double dtoe;        // time of ephemeris
	double accuracy;    // SV user range accuracy
	double dTOC;     	// time of clock (sec of week)
	int    nWeekCLK; 	// Week of clock (weeks)
	double dAfo;     	// 0 order clock correction (sec)
	double dAf1;     	// 1st order clock correction (sec/sec)
	double dAf2;     	// 2nd order clock correction (1/sec)
	double dCic;     	// inclination cosine term (rad)
	double dCis;     	// inclination sine term (rad)
	double dCrc;     	// radius cosine term (m)
	double dCrs;     	// radius sine term (m)
	double dCuc;     	// arg of lat cosine term (rad)
	double dCus;     	// arg of lat sine term (rad)
	double dAroot;   	// sqrt of major semiaxis (m)
	double dEcc;     	// eccentricity
	double d_w;      	// argument of perigee (rad)
	double dMo;      	// mean anomaly at ref time (rad)
	double dDn;      	// mean anomaly correction (rad/s)
	double dOMEGAo;  	// longitude of the ascending node (rad)
	double dOMEGADOT;	// rate of right ascension (rad/s)
	double d_io;     	// inclination angle at reference time (rad)
	double dIDOT;    	// inclination rate (rad/s)
	double dTGD;     	// broadcast time group delay (s)
	double updateTime;  // Time of update

} ephgps_t;

typedef struct
{
	char   vflg;        // Indicates if the ephemeris is available.
	int    prn;			// PRN of the GLONASS SV.
	int    freqslot;	// Frequency slot corresponding to the channel.

	int    nHealth;		// Health flag.
	double accuracy;	// SV user range accuracy.

	/* Only required in the RINEX version */
	int    nWeekNumber;
	double dtoe;

	double tb;			// Time of GLONASS ephemeris.
	unsigned int tk;	// Message frame time.

	double X_tb;		// Satellite PZ-90 X coordinate at tb (km).
	double Y_tb;		// Satellite PZ-90 Y coordinate at tb (km).
	double Z_tb;		// Satellite PZ-90 Z coordinate at tb (km).

	double Xdot_tb;		// Satellite PZ-90 X velocity at tb (km/sec).
	double Ydot_tb;		// Satellite PZ-90 Y velocity at tb (km/sec).
	double Zdot_tb;		// Satellite PZ-90 Z velocity at tb (km/sec).

	double Xdot2_tb;	// Satellite PZ-90 X acceleration at tb (km/sec2).
	double Ydot2_tb;	// Satellite PZ-90 Y acceleration at tb (km/sec2).
	double Zdot2_tb;	// Satellite PZ-90 Z acceleration at tb (km/sec2).

	double SV_freqdev;	// Relative frequency deviation (with respect to nominal)
	double SV_timedev;	// SV Clock bias correction (seconds).
	double dTGD;		// Time group delay between frequencies in SV (seconds).

	double tau_c;		// Time constant between GPS and GLONASS time systems (negligible).

	unsigned char P;    // Technological parameter from CS. Info for tau(GPS), tau(c).
	unsigned char P1;	// Immediate updating flag. Indicates a time interval between
						// two adjacent values of tb (minutes).
	unsigned char P2;	// Flag of oddness in tb parameter:  1-> odd, 0-> even.
	unsigned char P3;	// Flag indicating the number of satellite for which almanac
						// is transmitted: 1-> 5 sats, 0-> 4 sats.
	unsigned char P4;   // Flag to show that ephemeris are present (uploaded by CS).

	unsigned int ageDays;  // Age of the immediate information: time elapsed since the instant
						   // of its calculation (at CS) until the instant tb for n-satellite.
	unsigned int timeMark; // Time mark corresponding to the 0.3s message at the end of data.

} ephglo_t;

typedef enum
{
	GAL_INAV  = 0,
	GAL_PRS   = 1,
	GAL_FNAV  = 2,

} gal_service;

typedef struct
{
	char vflg;   	  // Indicates if the ephemeris are available.
	int  prn;         // PRN of the GPS SV.
	int  nWeekNumber; // Week Number.
	int  IODnav;      // Issue Of Data Clock.

	gal_service srv;  // Galileo service.

	// Galileo health information
	unsigned char E1B_hst; /* E1B health status */
	unsigned char E5a_hst; /* E5b health status */
	unsigned char E5b_hst; /* E5b health status */
	unsigned char E1A_hst; /* E1A health status */
	unsigned char E6A_hst; /* E6A health status */
    unsigned char E1B_dvs; /* E1B data validity */
    unsigned char E5a_dvs; /* E5b data validity */
    unsigned char E5b_dvs; /* E5b data validity */


	double sSISA;     // Signal-in-Space Accuracy

	int ToW;          // GST time of week
	int toeWeek;      // GST week number

	double dtoe;      // time of ephemeris
	double dTOC;      // time of clock (sec of week)
	int    nWeekCLK;  // Week of clock (weeks)
	double dAfo;      // 0 order clock correction (sec)
	double dAf1;      // 1st order clock correction (sec/sec)
	double dAf2;      // 2nd order clock correction (1/sec)
	double dCic;      // inclination cosine term (rad)
	double dCis;      // inclination sine term (rad)
	double dCrc;      // radius cosine term (m)
	double dCrs;      // radius sine term (m)
	double dCuc;      // arg of lat cosine term (rad)
	double dCus;      // arg of lat sine term (rad)
	double dAroot;    // sqrt of major semiaxis (m)
	double dEcc;      // eccentricity
	double d_w;       // argument of perigee (rad)
	double dMo;       // mean anomaly at ref time (rad)
	double dDn;       // mean anomaly correction (rad/s)
	double dOMEGAo;   // longitude of the ascending node (rad)
	double dOMEGADOT; // rate of right ascension (rad/s)
	double d_io;      // inclination angle at reference time (rad)
	double dIDOT;     // inclination rate (rad/s)

	// Group delay parameter (s)
	double bgd_E1E5a; /* Group delay parameter E1-E5a (s) */
	double bgd_E1E5b; /* Group delay parameter E1-E5b (s) */
	double bgd_E1E6A; /* Group delay parameter E1-E6A (s) */

	unsigned char rx_pages;

} ephgal_t;

typedef struct
{
	// TODO: implement definition of BeiDou ephemeris struct
	char   vflg;        // Indicates if the ephemeris is available.
	int    prn;		    // PRN of the BEI SV.
	int    ToW;         // GST time of week
	int    nWeekNumber; // Week Number.
	int    nIODE;       // Issue Of Data Ephemeris.
	int    nIODC;       // Issue Of Data Clock.
	int    nHealth;		// Health flag.
	double dtoe;        // time of ephemeris GPS time
	double dtoe_bei;	// time of ephemeris BEI time
	double accuracy;    // SV user range accuracy
	double dTOC;     	// time of clock (sec of week)
	int    nWeekCLK; 	// Week of clock (weeks)
	double dAfo;     	// 0 order clock correction (sec)
	double dAf1;     	// 1st order clock correction (sec/sec)
	double dAf2;     	// 2nd order clock correction (1/sec)
	double dCic;     	// inclination cosine term (rad)
	double dCis;     	// inclination sine term (rad)
	double dCrc;     	// radius cosine term (m)
	double dCrs;     	// radius sine term (m)
	double dCuc;     	// arg of lat cosine term (rad)
	double dCus;     	// arg of lat sine term (rad)
	double dAroot;   	// sqrt of major semiaxis (m)
	double dEcc;     	// eccentricity
	double d_w;      	// argument of perigee (rad)
	double dMo;      	// mean anomaly at ref time (rad)
	double dDn;      	// mean anomaly correction (rad/s)
	double dOMEGAo;  	// longitude of the ascending node (rad)
	double dOMEGADOT;	// rate of right ascension (rad/s)
	double d_io;     	// inclination angle at reference time (rad)
	double dIDOT;    	// inclination rate (rad/s)
	double dTGD_B1B3;   // broadcast time group delay (s)
	double dTGD_B2B3;   // broadcast time group delay (s)
} ephbei_t;

typedef union
{
	ephgps_t GPS;
	ephglo_t GLO;
	ephgal_t GAL;
	ephbei_t BEI;

} ephdata_t;



#define MAX_SAT_PRN_INSP3_SI (211)		  /* Maximum allowed PRN for a satellite in a SP3 File */

typedef struct
{
	//sp3 file
	int    psiPrnToIndex[MAX_SAT_PRN_INSP3_SI];
	double dStartAGPS;
	double dStepTime;
	int    siNRec;
	int    siNumberSats;

	//interpolation
	double timeTable[MAX_SAT_PRN_INSP3_SI][1000];
	double Pos_X_Table[MAX_SAT_PRN_INSP3_SI][1000];
	double Pos_Y_Table[MAX_SAT_PRN_INSP3_SI][1000];
	double Pos_Z_Table[MAX_SAT_PRN_INSP3_SI][1000];
	double Clock_Table[MAX_SAT_PRN_INSP3_SI][1000];
	int    nrec[MAX_SAT_PRN_INSP3_SI];
	int    validSP3Flag[MAX_SAT_PRN_INSP3_SI];

} SP3STRUCT_D;

//////////////////////////////////////////////////////////////////////////////
// satellite state propagation.
//
#define PROP_WINDOW 30

typedef struct
{
	double     tow_node [MAX_N_SAT_G];
	ephemsat_t sat      [MAX_N_SAT_G];

	int currChannel;

} propephem_t;
//////////////////////////////////////////////////////////////////////////////

typedef struct
{
	// Ionospheric Klobuchar model.
	double alpha[4];
	double beta [4];

	char   vflg;
	double A0;
	double A1;
	double tot;
	int    WNt;
	double dtls;
	int    WNlsf;
	double DN;
	double dtlsf;
	int    GPSweek;

} ionoutc_gps_t;

// Ionospheric and UTC data structure for I/NAV.
typedef struct
{
	// Position 0 reserved for NeQuick iono parameters (page 5 in OS, page 4 in PRS).
	// Position 1 reserved for UTC conversion parameters (page 6 in OS, page 5 in PRS).
	char vflg[2];

	/* Ionospheric model parameters */
	double         ai[3];            /* NeQuick parameters.  							*/
	unsigned char  StormFlags;       /* Ionospheric storm flags for 5 different regions */

	/* UTC conversion parameters. */
	double         A0, A1;           /* Coeffs for determining UTC time.                */
	unsigned long  tot;              /* Reference time for A0 & A1, sec of GPS week.    */
	short          dtls;             /* Cumulative past leap seconds.                   */
	unsigned       wnt;              /* Current UTC reference week number.              */
	unsigned       wnlsf;            /* Week number when dtlsf becomes effective.       */
	short          dn;               /* Day of week (1-7) when dtlsf becomes effective. */
	short          dtlsf;            /* Scheduled future leap second event.             */

	/* GPS to Galileo time offset (GGTO) and conversion parameters */
	double         A0g, A1g;         /* Coefficients for determining GGTO				*/
	unsigned long  tog;              /* Reference time for GGTO, sec of GAL week		*/
	unsigned       wng;              /* Reference week for for GGTO, GAL week number	*/

} ionoutc_gal_t;

typedef struct
{
	ionoutc_gps_t gps;
	ionoutc_gal_t gal;

} ionoutc_t;

//////////////////////////////////////////////////////////////////////////////
// Main structure for GPS/GLONASS (GNSS) computed data output

typedef struct
{
	int    week;
	double tow;

	// Arranged PRNs of satellites in the mask.
	int SAT_ID[MAX_CHANNELS_G];

	int noOfChannelsAv;

	ephdata_t EPH[MAX_N_SAT_G];
	obsdata_t OBS[MAX_CHANNELS_G];

	ionoutc_t iono_utc;

	char SP3_ON;

	// First-time PVT is computed.
	char first_fix;

	// According to consFlag_t definition
	int usedInPvt[NUM_SYS];

} gnssdata_t;

#endif

