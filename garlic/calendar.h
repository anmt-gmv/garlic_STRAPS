///////////////////////////////////////////////////////////////////////////////
//                           
// Copyright: GMV, S.A.
// Project  : SISNET
//
//
// Purpose  : Generic utilities for time format transformation.
//            Declaration module.
//
// Functions: GPSToCal_G, CalToGPS, _CalToJD, JDToCal,
//            GPSToJD(internal use only), CheckGPSWeek (internal use only),
//            ChJDFormat, DateToDay, DayToDate, CheckDate, LeapYear, GPStoUTC 
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//          |         |          |      |                                    |
//          |  1.0    | 30/04/97 | GNSS |            First version           |
//          |  1.0.1  | 20/01/00 | ejgg |         Changes for ECUREV.        |
//          |  1.1... | 01/02/05 | ejgg |         Update for ISAGNSS2.       |
//
///////////////////////////////////////////////////////////////////////////////
#ifndef _CALENDAR_H
#define _CALENDAR_H

typedef enum {JD, MJD1950, MJD2000, MJDIERS} JDFORMAT;

#define JD2MJD1950      (2433282.5)
#define JD2MJD2000      (2451544.5)
#define JD2MJDIERS      (2400000.5)


// Functions for GPS calendar system:
// ----------------------------------
void GPSToCal_G(long week, double time, int *year, int *month,
              int *day, int *hour, int *minute, double *second);
void CalToGPS_G(int year, int month, int day, int hour, int min,
              double sec, long *wn, double *gpstime);

// Julian dates management:
// -----------------------
long double ChJDFormat(long double julian, JDFORMAT InFormat,
                       JDFORMAT OutFormat);
long double _CalToJD(int year, int month, int day, int hour, int min,
                     double sec, JDFORMAT OutFormat);
void JDToCal(long double julian, JDFORMAT InFormat, int *year,
             int *month, int *day, int *hour, int *minute, double *second);

// Miscellaneous:
// --------------
int  DateToDay(int year, int month, int day);
void DayToDate(int year, int numday,int *dd, int *mm, int *yy);
void CheckDate(int *year,int *month,int *day,int *hour,int *min,double *sec);

char LeapYear(int year);

void GPStoUTC(int gpsWn, double gpsT, int wnt, double dtot, double a0,
              double a1, int wnlsf, double ddn, double ddtls, double ddtlsf, 
              int *utcWn, double *utcT);

#endif //_CALENDAR_H
