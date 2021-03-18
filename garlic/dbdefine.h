#ifndef _DBDEFINE_HPP
#define _DBDEFINE_HPP

// Maximum number of channels available
#define MAX_CHANNELS_G  200

// Defined systems.
#define NUM_SYS		4
#define NUM_FREQ	4
#define TWO_RX		0

// Limits
#define MAXNUMGPS_G     (32) // maximum number of GPS     satellites
#define MAXNUMGAL_G		(36) // maximum number of Galileo satellites
#define MAXNUMGLO_G		(24) // maximum number of GLONASS satellites
#define MAXNUMBEI_G		(37) // maximum number of BeiDou  satellites

#define MAX_N_SAT_G     (MAXNUMGPS_G + MAXNUMGAL_G + MAXNUMGLO_G + MAXNUMBEI_G)

typedef enum
{
	ION_NONE   = 0,
	ION_AUTO   = 1,
	ION_KLOB   = 2,
	ION_NEQU   = 3,
	ION_FREE   = 4
} ionoModel_t;

// Frequency constants.
#define BANDFREQL1  1575420000.0
#define BANDFREQL5  1176450000.0
#define BANDFREQE5a 1176450000.0
#define BANDFREQE5b 1207140000.0
#define BANDFREQL6  1278750000.0
#define BANDFREQB1  1561098000.0
#define ONECODEFREQ 1023000.0

#define LAMBDA_GPS_L1  (SPEED_OF_LIGHT/BANDFREQL1)
#define LAMBDA_GPS_L5  (SPEED_OF_LIGHT/BANDFREQL5)
#define LAMBDA_GAL_E5a (SPEED_OF_LIGHT/BANDFREQE5a)
#define LAMBDA_GAL_E5b (SPEED_OF_LIGHT/BANDFREQE5b)
#define LAMBDA_GAL_E6  (SPEED_OF_LIGHT/BANDFREQL6)
#define LAMBDA_BEI_B1  (SPEED_OF_LIGHT/BANDFREQB1)

#define SBFACTOR_E1A  (1-15*ONECODEFREQ/BANDFREQL1)
#define SBFACTOR_E1BC (1-ONECODEFREQ/BANDFREQL1)

// L1 wavelengths in GLONASS, in terms of the K value (f_GLO = fc + K*Af0).
static const double LAMBDA_GLO_L1[14] = {0.187597455043216, 0.187531446086481, 0.187465483565873,
        							     0.187399567432411, 0.187333697637180, 0.187267874131334,
        							     0.187202096866097, 0.187136365792759, 0.187070680862681,
        							     0.187005042027290, 0.186939449238084, 0.186873902446626,
        							     0.186808401604549, 0.186742946663552};

// GPS constants
#define PI_G			3.1415926535898
#define R_EARTH			6378137.0
#define MU_EARTH		3.986005e14
#define SPEED_OF_LIGHT	299792458.0
#define OMEGA_EARTH		7.2921151467e-5
#define FLATTE_G		1.0/298.257223563
#define E2_G			FLATTE_G*(2.0-FLATTE_G)
#define E12_G			E2_G/(1.-E2_G)

// GLONASS constants
#define MU_EARTH_GLO      398600.4418 		/* mu GLONASS : km3/s2 */
#define OMEGA_EARTH_GLO   7.2921150e-5  /* omega GLONASS  7292115.1467e-11    */
#define QUARTHOUR_GLO     900. 				/* seconds in 15 minutes */
#define UTCSU_OFFSET_GLO  10800. 			/* seconds in 3 hours */
#define RK_STEP_GLO       60. 				/* step in runge-kutta */
#define MARGIN_TB_GLO     1800.    		    /* validity of GLONASS ephem */
#define AE_GLO            6378.136 			/* Earth radius GLONASS */
#define C20_GLO           -1082.62575e-6 	/* C20 GLONASS */

// Galileo constants
#define MU_EARTH_GAL      3.986004418e14

// BeiDou constants
#define MU_EARTH_BEI      3.986004418e14    /* mu BeiDou : m3/s2 */
#define OMEGA_EARTH_BEI   7.2921150e-5      /* omega BeiDou rad/s */

// Ancillary parameters
#define EPSILON			1e-10
#define DAYSECS         86400
#define HALFWEEK        302400
#define HALFDAY         43200
#define SECONDSONWEEK	604800
#define GRAVITY_G		9.801
#define DEGTORAD        (PI_G/180.0)
#define RADTODEG        (180.0/PI_G)

//////////////////////////////////////////////////////////////////////////////

#endif

