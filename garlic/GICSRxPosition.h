///////////////////////////////////////////////////////////////////////////////
//
// Copyright: GMV S.A.
// Project  : SISNET
//
// Purpose  : Useful subroutines for getting user position
//            Declaration module.
//
// Functions: UserPos_WGS
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//          |         |          |      |                                    |
//          |  1.0.1  | 17/11/99 | ejgg |            First version           |
//          |  1.0.1  | 07/01/01 | rrrl |       Minor changes for ARMHADE    |
//          |  1.1... | 01/02/05 | ejgg |         Update for ISAGNSS2.       |
//          |  1.1... | 14/04/09 | iimp |         Update for Kalman.         |
//			|  2.0    | 12/01/11 | cmvv | 	   Migrated from UserPos.h       |
//////////////////////////////////////////////////////////////////////////////

#ifndef GICSRxUserPosH
#define GICSRxUserPosH

#include "dbtypedef.h"
#include "GICSRx_defines.h"

char GICSRxGPSPos(ephgps_t *ephem, double RxTime, double TxTime, ephemsat_t *pephsat, double coef_TGD);

char GICSRxGalPos(ephgal_t *ephem, double RxTime, double TxTime, ephemsat_t *pephsat, double coef_TGD);

char GICSRxGLOPos(ephglo_t *ephem, double RxTime, double TxTime, ephemsat_t *pephsat, double coef_TGD);

char GICSRxBEIPos(ephbei_t *ephem, double RxTime, double TxTime, ephemsat_t *pephsat, double coef_TGD);

double GICSRxSatRange(const double dUserPos[3], const double dSVPos[3], float unit[3]);

float  GICSRxSatDoppler(const double dUserVel[3], const double dSVVel[3], const float unit[3], const double dUserPos[3]);

#if CALC_ATMSPHCORR == 1
void GICSRxSatElevAzim(float LatLonTrig[4], float Rsu[3], float *El, float *Az);
#else
void GICSRxSatElev(float LatLonTrig[4], float Rsu[3], float *El);
#endif

char GICSRxSatPos_SP3(SP3STRUCT_D *pSp3, int PRN, double RxTime, double TxTime, ephemsat_t *pephsat);

char GICSRxReadSP3File(SP3STRUCT_D *sp3_P, char *sp3_filename);

#endif
