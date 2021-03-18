/*
 * load_ARAIM_parameters.c
 *
 *  Created on: 13 ago-> 2019
 *      Author: anmt
 */

#include "ARAIM_RINEX.h"

/*-------------------------------------------------------------------------------------------------------------------------------------------------
		                                                            load_ARAIM_parameters.c
-------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Function to load the parameters for ARAIM simulation
 *
 * INPUT PARAMETERS:
 *
 * 		- none
 *
 * OUTPUT PARAMETERS:
 *
 * 		- ARAIM:			struct containing all the ARAIM configuration parameters
 */

void load_ARAIM_parameters( configuration* ARAIM ){


	/* Values suggested by WG-C Advanced RAIM Technical Subgroup Reference Airborne Algorithm Description Document:
	 * 	- b_nom = .75m;
	 * 	- Psat = 1e-5;
	 * 	- Pconst = 1e-4 (GPS), 1e-8 (GAL); */

	ARAIM->Psat_GAL = 1e-5;
	ARAIM->Psat_GPS = 1e-5;
	ARAIM->Pconst_GAL = 1e-4;;
	ARAIM->Pconst_GPS = 1e-8;
	ARAIM->bnom_GAL = .75;
	ARAIM->bnom_GPS = .75;

	/* error model for user contribution to error (elevation dependent) */

	ARAIM->UERE_usr_GAL_deg[0] = 0; ARAIM->UERE_usr_GAL_deg[1] = 5; ARAIM->UERE_usr_GAL_deg[2] = 10; ARAIM->UERE_usr_GAL_deg[3] = 15; ARAIM->UERE_usr_GAL_deg[4] = 20;
	ARAIM->UERE_usr_GAL_deg[5] = 25; ARAIM->UERE_usr_GAL_deg[6] = 30; ARAIM->UERE_usr_GAL_deg[7] = 35; ARAIM->UERE_usr_GAL_deg[8] = 40; ARAIM->UERE_usr_GAL_deg[9] = 45;
	ARAIM->UERE_usr_GAL_deg[10] = 50; ARAIM->UERE_usr_GAL_deg[11] = 55; ARAIM->UERE_usr_GAL_deg[12] = 60; ARAIM->UERE_usr_GAL_deg[13] = 65; ARAIM->UERE_usr_GAL_deg[14] = 70;
	ARAIM->UERE_usr_GAL_deg[15] = 75; ARAIM->UERE_usr_GAL_deg[16] = 80; ARAIM->UERE_usr_GAL_deg[17] = 85; ARAIM->UERE_usr_GAL_deg[18] = 90;

	ARAIM->UERE_usr_GAL[0] = 0; ARAIM->UERE_usr_GAL[1] = 0.4529; ARAIM->UERE_usr_GAL[2] = 0.3553; ARAIM->UERE_usr_GAL[3] = 0.3063; ARAIM->UERE_usr_GAL[4] = 0.2638;
	ARAIM->UERE_usr_GAL[5] = 0.2593; ARAIM->UERE_usr_GAL[6] = 0.2555; ARAIM->UERE_usr_GAL[7] = 0.2504; ARAIM->UERE_usr_GAL[8] = 0.2438; ARAIM->UERE_usr_GAL[9] = 0.2396;
	ARAIM->UERE_usr_GAL[10] = 0.2359; ARAIM->UERE_usr_GAL[11] = 0.2339; ARAIM->UERE_usr_GAL[12] = 0.2302; ARAIM->UERE_usr_GAL[13] = 0.2295; ARAIM->UERE_usr_GAL[14] = 0.2278;
	ARAIM->UERE_usr_GAL[15] = 0.2297; ARAIM->UERE_usr_GAL[16] = 0.2310; ARAIM->UERE_usr_GAL[17] = 0.2274; ARAIM->UERE_usr_GAL[18] = 0.2277;

	/* Values suggested by WG-C Advanced RAIM Technical Subgroup Reference Airborne Algorithm Description Document:
	 * 	- sigma_URA = .5m, .75m, 1m, 1.5m, 2m (LPV-200), 2.5m (LPV-250);
	 * 	- sigma_URE = 2 / 3 * sigma_URA; */

	double sigma_ref = 2.5;

	ARAIM->SISE_GAL = ( 2.0 / 3.0 ) * sigma_ref;
	ARAIM->ura_GAL = sigma_ref;
	ARAIM->URE_GPS = ( 2.0 / 3.0 ) * sigma_ref;
	ARAIM->ura_GPS = sigma_ref;

	/* Values suggested by WG-C Advanced RAIM Technical Subgroup Reference Airborne Algorithm Description Document:
	 * 	- PHMI (total, VERT + HOR) = 1e-7;
	 * 	- PHMI_VERT = 9.8e-8 (LPV-200, LPV-250);
	 * 	- Pemt = 1e-5;
	 * 	- Pfa (total, VERT + HOR) = 4e-6;
	 * 	- Pfa_VERT = 3.9e-6 (LPV-200, LPV-250);
	 * 	- Pfa_HOR = 9e-8 (LPV-200, LPV-250), 1e-7 (RNP);
	 * 	- Psat_THRES = 8e-8 (LPV-200, LPV-250), 4e-8 (RNP);
	 * 	- TOL_PL = 5e-2 */

	ARAIM->TOL_PL = 5e-2;
	ARAIM->Psat_THRES = 4e-8; // default 4e-8, 1e-8 for N_f_max = 2, 1e-12 for N_f_max = 3, 1e-18 for N_f_max = 4, 1e-24 for N_f_max = 5
	ARAIM->Phmi_VERT = 2e-9;
	ARAIM->Phmi_HOR = 9.8e-8;
	ARAIM->Pfa_VERT = 3.9e-6;
	ARAIM->Pfa_HOR = 9e-8;

	/* Values suggested by WG-C Advanced RAIM Technical Subgroup Reference Airborne Algorithm Description Document:
	 * 	- KACC = 1.96;
	 * 	- KFF = 5.33; */

	ARAIM->KACC = 1.96;
	ARAIM->KFF = 5.33;

	ARAIM->Pemt = 1e-5;

	/* Values suggested by WG-C Advanced RAIM Technical Subgroup Reference Airborne Algorithm Description Document:
	 * 	- VAL = 35m (LPV-200), 50m (APV 1, LPV-250);
	 * 	- HAL = 40m (LPV-200, APV 1, LPV-250 ), 185m (RNP 0.1), 556m (RNP 0.3);
	 * 	- EMT = 15m (LPV-200);
	 * 	- sigma_acc = 1.87 m */

	ARAIM->std_v_req =  1.87;
	ARAIM->std_h_req =  2;
	ARAIM->ACCthres = 4; // (approx. 1.96 * 1.87)
	ARAIM->ACCffthres =   10; // (approx. 5.33 * 1.87)
	ARAIM->EMTthres = 15;
	ARAIM->HAL = 40;
	ARAIM->VAL = 35;

	ARAIM->Pfa_CHI2 = 1e-8;

#if ARAIM_ALG == 3 // For validation wrt Stanford Paper
	ARAIM->bnom_GAL = .5;
	ARAIM->bnom_GPS = .5;

	ARAIM->Psat_GAL = 1e-4;
	ARAIM->Psat_GPS = 1e-4;

	ARAIM->Pconst_GPS = 1e-5;
	ARAIM->Pconst_GAL = 1e-5;

	ARAIM->SISE_GAL = .5;
	ARAIM->ura_GAL = .75;
	ARAIM->URE_GPS = .5;
	ARAIM->ura_GPS = .75;

	ARAIM->Psat_THRES = 9e-8;
	ARAIM->Phmi_VERT = 9e-8;
	ARAIM->Phmi_HOR = 1e-8;
#endif

#if ARAIM_ALG == 4 // For validation wrt WG-C ARAIM Document
	ARAIM->bnom_GAL = .5;
	ARAIM->bnom_GPS = .5;

	ARAIM->Psat_GAL = 1e-5;
	ARAIM->Psat_GPS = 1e-5;

	ARAIM->Pconst_GPS = 1e-4;
	ARAIM->Pconst_GAL = 1e-4;

	ARAIM->SISE_GAL = .5;
	ARAIM->ura_GAL = .75;
	ARAIM->URE_GPS = .5;
	ARAIM->ura_GPS = .75;

	ARAIM->Psat_THRES = 8e-8;
	ARAIM->Phmi_VERT = 9.8e-8;
	ARAIM->Phmi_HOR = 2e-9;
#endif

#if ARAIM_ALG == 5 || ARAIM_ALG == 7 // for validation against SVS ARAIM (without and with approximate solution)
	ARAIM->bnom_GAL = .5;
	ARAIM->bnom_GPS = .5;
	ARAIM->SISE_GAL = 6;
	ARAIM->ura_GAL = 9;
	ARAIM->URE_GPS = 6;
	ARAIM->ura_GPS = 9;
	ARAIM->TOL_PL = 1e-2;
	ARAIM->Pfa_VERT = 9e-8;
	ARAIM->Pfa_HOR = 3.9e-6;
	ARAIM->HAL = 80;
	ARAIM->VAL = 0;
	ARAIM->ACCthres = 4; // (1.96 * 1.87)
	ARAIM->ACCffthres =  0; // (5.33 * 1.87)
	ARAIM->EMTthres = 0;
#endif

}

