///////////////////////////////////////////////////////////////////////////////
//                           
// Copyright: GMV, S.A.
// Project  : SISNET
//
//
// Purpose  : Generic utilities for time format transformation.
//
// Functions: GPSToCal_G, CalToGPS_G, _CalToJD, JDToCal,
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
#include <math.h>
#include <stdio.h>

#include "dbdefine.h"
#include "calendar.h"

static long        CheckGPSWeek(long week);
static long double GPSToJD(long week, double time);

///////////////////////////////////////////////////////////////////////////////
//
// Function : GPSToCal_G
// Purpose  : Change from GPS system time to calendar.
//
// Args I   : GPSweek   (long    )  GPS week number.           (0-  1024)
//            GPStime   (double  )  GPS time within the week.  (0-604800)
// Args I/O :
// Args O   : year      (int *   )  Year.                      (e.g.1996)
//            month     (int *   )  Month.                     (1-12)
//            day       (int *   )  Day.                       (1-31)
//            hour      (int *   )  Hours.                     (0-23)
//            minute    (int *   )  Minutes.                   (0-59)
//            second    (double *)  Seconds.                   (0-60)
// Returns  :           void
// Depends  :
// Calls    : GPSToJD, JDToCal.
// Comments : 
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//          |         |          |      |                                    |
//          |   1.0   | 97/04/30 | GNSS |            First version           |
//          |   1.1   | 97/10/13 | mafv |       Add CheckDate() function     |
//
///////////////////////////////////////////////////////////////////////////////
void GPSToCal_G( long GPSweek, double GPStime, int *year, int *month,
		int *day, int *hour, int *minute, double *second)
{
	// Performs conversion with rounded seconds and then add the fractional part
	// -------------------------------------------------------------------------
	// The goal is to prevent loss of precision in floating point ops.

	double timeInt = 0.0;
	long double jul = 0.0;
	double timeFrac = modf(GPStime, &timeInt);
	const double SMALLDELTA = 1.e-3;

	if(timeFrac > SMALLDELTA) {
		timeInt += SMALLDELTA;
	}
	jul = GPSToJD(GPSweek, timeInt);
	JDToCal(jul, JD, year, month, day, hour, minute, second);
	*second = floor(*second + 0.5) + timeFrac;
	if (timeFrac > SMALLDELTA){
		*second -= SMALLDELTA;
	}
	CheckDate(year, month, day, hour, minute, second);
}

///////////////////////////////////////////////////////////////////////////////
//
// Function : CalToGPS_G
// Purpose  : Change from calendar to GPS system time
//
// Args I   : year      (int      )  Year.                      (e.g.1996)
//            month     (int      )  Month.                     (1-12)
//            day       (int      )  Day.                       (1-31)
//            hour      (int      )  Hours.                     (0-23)
//            minute    (int      )  Minutes.                   (0-59)
//            second    (double   )  Seconds.                   (0-60)
// Args I/O :
// Args O   : GPSweek   (long *   )  GPS week number.           (0-  1024)
//            GPStime   (double * )  GPS time within the week.  (0-604800)

// Returns  :           void
// Depends  :
// Calls    : GPSToJD, JDToCal.
// Comments : 
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//          |         |          |      |                                    |
//          |   1.0   | 97/04/30 | GNSS |            First version           |
//          |   1.1   | 97/10/13 | mafv |       Add CheckDate() function     |
//
///////////////////////////////////////////////////////////////////////////////
void CalToGPS_G(int year, int month, int day, int hour, int min,
		double sec, long *GPSweek, double *GPStime)
{
	long double jul;

	jul      = _CalToJD(year, month, day, hour, 0, 0.0, JD);
	*GPSweek = (long)((jul - 2444244.5) / 7.0); //CCBG for linux (cast to long)
	*GPStime = ((int)(jul-2444244.5)-*GPSweek*7.0)*24*3600.0 + hour*3600.0+min*60.0+sec;
}

///////////////////////////////////////////////////////////////////////////////
//
// Function : _CalToJD
// Purpose  : Change from calendar to Julian date
//
// Args I   : year      (int      )  Year.                      (e.g.1996)
//            month     (int      )  Month.                     (1-12)
//            day       (int      )  Day.                       (1-31)
//            hour      (int      )  Hours.                     (0-23)
//            minute    (int      )  Minutes.                   (0-59)
//            second    (double   )  Seconds.                   (0-60)
//            OutFormat (JDFORMAT )  Target Julian date format.
// Args I/O :
// Args O   : 
// Returns  :           long double  Julian date in days.
// Depends  :
// Calls    : ChJDFormat
// Comments : Julian date is returned as a long double to provide enough
//            accuracy.
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//          |         |          |      |                                    |
//          |   1.0   | 97/04/30 | GNSS |            First version           |
//
///////////////////////////////////////////////////////////////////////////////
long double _CalToJD(int year, int month, int day, int hour, int min, double sec,
		JDFORMAT OutFormat)
{
	int          ja, jy=year, jm;
	const long   IGREG = (15+31L*(10+12L*1582));
	long  double Julian, dd;

	dd = (double)(day) +(double)(hour)/24.0 + (double)(min)/1440.0 + sec/86400.0;

	if(jy < 0) {
		++jy;
	}

	if(month > 2) {
		jm=month+1;
	}
	else
	{    --jy;
	jm=month+13; }

	Julian = floor(365.25*jy) + floor(30.6001*jm) + dd + 1720994.5;

	if(day+31L*(month+12L*year) >= IGREG)
	{
		ja      = (int)(0.01*jy);
		Julian += 2-ja+(int) (0.25*ja);
	}

	return(ChJDFormat(Julian, JD, OutFormat));
}

///////////////////////////////////////////////////////////////////////////////
//
// Function : JDToCal
// Purpose  : Change from Julian date to calendar.
//
// Args I   : Julian    (long double) Julian date
//            InFormat  (JDFORMAT   ) Source Julian date format.
// Args I/O :
// Args O   : year      (int    *  )  Year.                      (e.g.1996)
//            month     (int    *  )  Month.                     (1-12)
//            day       (int    *  )  Day.                       (1-31)
//            hour      (int    *  )  Hours.                     (0-23)
//            minute    (int    *  )  Minutes.                   (0-59)
//            second    (double *  )  Seconds.                   (0-60)
// Returns  : void
// Depends  :
// Calls    : ChJDFormat
// Comments : 
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//          |         |          |      |                                    |
//          |   1.0   | 30/04/97 | GNSS |            First version           |
//          |   2.0   | 28/08/00 | GNSS |        Upgrade used in ECIES       |
//
///////////////////////////////////////////////////////////////////////////////
void JDToCal(long double Julian, JDFORMAT InFormat, int *year, int *month, int *day,
		int *hour, int *minute, double *second)
{
	long int iDay, i, j, k, l;
	long double sFrac;

	Julian = ChJDFormat(Julian, InFormat, MJDIERS);

	iDay = (long int)Julian - 33282;
	i    = (4000*(iDay + 18204))/1461001;
	j    = iDay - (1461*i)/4 + 18234;
	k    = (80*j)/2447;
	l    = k/11;

	*year  = (int)(1900 + i + l);
	*month = (int)(k + 2 - 12*l);
	*day   = (int)(j - (2447*k)/80);

	sFrac   = (Julian - (int)Julian)*24.0;
	*hour   = (int)sFrac;
	sFrac   = (sFrac - *hour)*60.0;
	*minute = (int)sFrac;
	sFrac   = (sFrac - *minute)*60.0;
	*second = (double)sFrac;
	sFrac   = sFrac - *second;

	// Correction for floating point instability
	*second = floor(*second + 0.5);

	// Check date before returning values
	CheckDate(year, month, day, hour, minute, second);
}

///////////////////////////////////////////////////////////////////////////////
//
// Function : ChJDFormat
// Purpose  : Change Julian date format.
//
// Args I   : Julian    (long double) Julian date
//            InFormat  (JDFORMAT   ) Source Julian date format.
//            OutFormat (JDFORMAT )  Target Julian date format.
// Args I/O :
// Args O   : year      (int    *  )  Year.                      (e.g.1996)
//            month     (int    *  )  Month.                     (1-12)
//            day       (int    *  )  Day.                       (1-31)
//            hour      (int    *  )  Hours.                     (0-23)
//            minute    (int    *  )  Minutes.                   (0-59)
//            second    (double *  )  Seconds.                   (0-60)
// Returns  :           long double  Julian date in days.
// Depends  :
// Calls    :
// Comments : Ref: OAD Standards: Time and Coordinate Systems for ESOC Flight
//            Dynamics Operations. Issue 1, May 1994.
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//          |         |          |      |                                    |
//          |   1.0   | 97/07/16 | GNSS |            First version           |
//
///////////////////////////////////////////////////////////////////////////////
long double ChJDFormat(long double Julian, JDFORMAT InFormat, JDFORMAT OutFormat)
{
	if(InFormat == OutFormat) return(Julian);

	switch(InFormat)
	{
	case JD:
		switch(OutFormat)
		{
			case MJD1950:   return(Julian-JD2MJD1950);
			case MJD2000:   return(Julian-JD2MJD2000);
			case MJDIERS:   return(Julian-JD2MJDIERS);
			default:        return(0); //CCBG for linux
		}
		break;

	case MJD1950:
		switch(OutFormat)
		{
			case JD:        return(Julian+JD2MJD1950);
			case MJD2000:   return(Julian+JD2MJD1950-JD2MJD2000);
			case MJDIERS:   return(Julian+JD2MJD1950-JD2MJDIERS);
			default:        return(0); //CCBG for linux
		}
		break;

	case MJD2000:
		switch(OutFormat)
		{
			case JD:        return(Julian+JD2MJD2000);
			case MJD1950:   return(Julian+JD2MJD2000-JD2MJD1950);
			case MJDIERS:   return(Julian+JD2MJD2000-JD2MJDIERS);
			default:        return(0); //CCBG for linux
		}
		break;

	case MJDIERS:
		switch(OutFormat)
		{
			case JD:        return(Julian+JD2MJDIERS);
			case MJD1950:   return(Julian+JD2MJDIERS-JD2MJD1950);
			case MJD2000:   return(Julian+JD2MJDIERS-JD2MJD2000);
			default:        return(0); //CCBG for linux
		}
		break;

		default: return(0); //CCBG for linux
	}
}


///////////////////////////////////////////////////////////////////////////////
//
// Function : GPSToJD
// Purpose  : Change from GPS system time to Julian date.
//
// Args I   : GPSweek  (long    )    GPS week number.          (0-  1024)
//            GPStime  (double  )    GPS time within the week. (0-604800)
// Args I/O :
// Args O   :
// Returns  : static long double      Julian date in days.
// Depends  :
// Calls    : GheckGPSWeek
// Comments : Julian date is returned as a long double to provide enough
//            accuracy. Internal use only.
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//          |         |          |      |                                    |
//          |   1.0   | 97/04/30 | GNSS |            First version           |
//
///////////////////////////////////////////////////////////////////////////////
static long double GPSToJD(long week, double time)
{
	week = CheckGPSWeek(week);
	return (long double)(2444244.5) + week*7 + time/86400.0;
}

///////////////////////////////////////////////////////////////////////////////
//
// Function : CheckGPSWeek
// Purpose  : Checks GPS week and adds 1024 weeks if necessary (GPS roll over)
//
// Args I   : GPSweek  (long    )    GPS week number.          (0-  1024)
// Args I/O :
// Args O   :
// Returns  : static long       Input GPS week + 1024 if necessary.
// Depends  :
// Calls    :
// Comments : 
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//          |         |          |      |                                    |
//          |   1.0   | 97/04/30 | GNSS |            First version           |
//
///////////////////////////////////////////////////////////////////////////////
static long CheckGPSWeek(long GPSweek)
{
	// The GPS week is broadcast as modulo 1024
	// ----------------------------------------
	if (GPSweek < 750) GPSweek += 1024;
	return GPSweek;
}

///////////////////////////////////////////////////////////////////////////////
//
// Function : DateToDay
// Purpose  : Changes from date to day of the year.
//
// Args I   :  year     (int   )  Year.    (e.g.1996)
//             month    (int   )  Month.     (1-12)
//             day      (int   )  Day.       (1-31)
// Args I/O :
// Args O   :
// Returns  : int       Day of the year (1-366), counted from January 1st.
// Depends  :
// Calls    : LeapYear
// Comments : Takes into account if the input year is a leap-year. 
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//          |         |          |      |                                    |
//          |   1.0   | 97/04/30 | GNSS |            First version           |
//
///////////////////////////////////////////////////////////////////////////////
int DateToDay(int year, int month, int day)
{
	int i=0, yday=day;
	int monthDays[] =  { 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

	for(i=0; i<month; i++) yday += monthDays[i];
	if((LeapYear(year) == 1) && (month > 2)) yday++;

	return yday;
}

///////////////////////////////////////////////////////////////////////////////
//
// Function : DayToDate
// Purpose  : Changes from day of the year to date.
//
// Args I   :  year     (int   )  Year.                (e.g.1996)
//             numday   (int  )   Day of the year.     (1-366)
// Args I/O :
// Args O   : dd        (int *)   Year modulus 100.    (e.g.96)
//            mm        (int *)   Month.               (1-12)
//            yy        (int *)   Day.                 (1-31)
// Returns  : void
// Depends  :
// Calls    : LeapYear
// Comments : Takes into account if the input year is a leap-year. 
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//          |         |          |      |                                    |
//          |   1.0   | 97/04/30 | GNSS |            First version           |
//
///////////////////////////////////////////////////////////////////////////////
void DayToDate(int yy, int numday, int *dd, int *mm, int *year)
{
	int monthDays[] = { 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
	int yday=monthDays[1], i=1;

	if(yy<90) *year = 2000 + yy;
	else      *year = 1900 + yy;

	if(LeapYear(*year) == 1) monthDays[2]++;

	while(numday>yday) {
		i++;
		yday += monthDays[i];
	}

	*mm = i;
	*dd = numday - (yday - monthDays[i]);
}

///////////////////////////////////////////////////////////////////////////////
//
// Function : CheckDate
// Purpose  : Checks that year, month (etc) values are within the correct
//            ranges. If necessary, repairs the date.
//
// Args I   : year      (int *   )  Year.                      (e.g.1996)
//            month     (int *   )  Month.                     (1-12)
//            day       (int *   )  Day.                       (1-31)
//            hour      (int *   )  Hours.                     (0-23)
//            minute    (int *   )  Minutes.                   (0-59)
//            second    (double *)  Seconds.                   (0-60)
// Args I/O :
// Args O   : 
// Returns  : void
// Depends  :
// Calls    : LeapYear
// Comments : Takes into account if the input year is a leap-year. 
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//          |         |          |      |                                    |
//          |   1.0   | 97/07/15 | GNSS |            First version           |
//
///////////////////////////////////////////////////////////////////////////////
void CheckDate(int *year, int *month, int *day, int *hour, int *min, double *sec)
{
	int    monthDays[] = { 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
	int    idelta;
	double ddelta;

	// Checks and repairs date:
	// ------------------------
	if(*sec >= 60.0)
	{
		ddelta = *sec-60.0;
		*sec   = fmod(ddelta,60.0);
		*min  += (int)(1+ddelta/60.0);
	}
	if(*min >= 60)
	{
		idelta = *min-60;
		*min   = idelta%60;
		*hour += 1+idelta/60;
	}
	if(*hour >= 24)
	{
		idelta = *hour-24;
		*hour  = idelta%24;
		*day  += 1+idelta/24;
	}

	while(*day > monthDays[*month])
	{
		*day -=  monthDays[*month];
		(*month)++;

		if(*month > 12)
		{
			*month = 1;
			(*year)++;

			// Checks if the current year is a leap year:
			// ------------------------------------------
			if(LeapYear(*year) == 1) monthDays[2] = 29;
			else                      monthDays[2] = 28;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
//
// Function : LeapYear
// Purpose  : Checks if a year is a leap-year consisting of 366 days.
//
// Args I   : year      (int   )  Year.                      (e.g.1996)
// Args I/O :
// Args O   : 
// Returns  : bool      true if year is a leap-year, else returns false.
// Depends  :
// Calls    : 
// Comments : Takes into account if the input year is a leap-year. 
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//          |         |          |      |                                    |
//          |   1.0   | 97/07/15 | GNSS |            First version           |
//
///////////////////////////////////////////////////////////////////////////////
char LeapYear(int year)
{
	char leap = 0;

	if(((year%4 == 0) && (year%100 != 0)) || (year%400 == 0))
		leap = 1;
	else
		leap = 0;

	return leap;
}

///////////////////////////////////////////////////////////////////////////////
//
// Function : GPStoUTC
// Purpose  : Changes to UTC time.
//
//  Args I  : gpsWn     (i) gps week in GPS frame.
//            gpsT      (i) gps second in GPS frame (seconds of week).
//            wnt       (i) UTC reference week number (weeks).
//            dtot      (i) reference time for UTC data (sec).
//            a0        (i) Constant term for UTC correction (sec).
//            a1        (i) 1st order coefficient for UTC correction (sec/sec)
//            wnlsf     (i) Week number at which the leap second becomes
//                          effective (weeks).
//            ddn       (i) day number (days).
//            ddtls     (i) delta time due to leap seconds (sec).
//            ddtlsf    (i) effective delta time due to leap seconds (sec).
// Args I/O :
// Args O   : utcWn     (o) GPS week in the UTC frame
//            utcT      (o) GPS second in the UTC frame
// Returns  : 
// Depends  :
// Calls    : 
// Comments : See inside.
// History  :
//          | VERSION |   DATE   | NAME |               CHANGE               |
//          |         |          |      |                                    |
//          |   1.0   | 00/08/10 | ejgg |            First version           |
//          |   2.0   | 04/02/02 | aaag |         Upgrade used in ASQF       |
//
///////////////////////////////////////////////////////////////////////////////
void GPStoUTC(int gpsWn, 
		double gpsT,
		int wnt,
		double dtot,
		double a0,
		double a1,
		int wnlsf,
		double ddn,
		double ddtls,
		double ddtlsf,
		int *utcWn,
		double *utcT)
{

	/*
	 *  NOTE: The notation and formulation is followed from
	 *           "Navstar GPS Space Segment / Navigation User Interfaces."
	 *
	 * Depending upon the relationship of the effectivity date to the user's
	 * current GPS time, the following three different UTC/GPS-time
	 * relationships exist. As explained in NOTE 2, a simplified version has
	 * been implemented an only two different relationship are considered.
	 *
	 *   (a) Whenever the effectivity time indicated by the wnlsf and the ddn
	 *       values is not in the past (relative to the user's present time),
	 *       and the user's present time does not fall in the timespan
	 *       which starts at ddn + 3/4 and ends at ddn + 5/4.
	 *
	 *   (b) Whenever the user's current time falls within the timespan
	 *       of ddn + 3/4 to ddn + 5/4.

	 *   (c) Whenever the effectivity time of the leap second event, as indicated
	 *       by the wnlsf and ddn values, is in the "past" (relative to the
	 *       user's current time).
	 *
	 */

	// VARIABLE DECLARATION
	// Delta time to be apllicated to reansform from GPS time to UTC time.
	double deltaUTC = 0.0;

	// Compute time difference between user's time (gpsWn, gpsT) and (wnlsf, ddn)
	double timeDiff = 0.0;

	// Intermediate factor neccesary for the computations.
	double w = 0.0;

	// BEGINNING OF CODE
	// Compute time difference between user's time (gpsWn, gpsT) and (wnlsf, ddn)
	timeDiff = (double)( (gpsWn - wnlsf) * SECONDSONWEEK )   + gpsT -
			(ddn-1)*DAYSECS, /* ddn is input in days.
			 * Day 1 is the first day
			 * before the end/start of
			 * the week.*/

			// Convert this difference to days.
			timeDiff = timeDiff / (double)DAYSECS;

	// Initialise GPS week in the UTC frame
	*utcWn = gpsWn;

	// Check which case (a), (b) or (c) is given

	if ( timeDiff < 3.0/4.0 ) { // Case (a) is given

		// Compute deltaUTC
		deltaUTC = ddtls +
				a0 + a1 * ( gpsT - dtot + (double)(SECONDSONWEEK * (gpsWn - wnt)) );

		// Compute GPS second in the UTC frame (in seconds of week).
		*utcT = gpsT - deltaUTC;

	}// End of is case (a) is given

	else if ( timeDiff < 5.0/4.0 ) { // Case (b) is given

		// Compute deltaUTC
		deltaUTC = ddtls +
				a0 + a1 * ( gpsT - dtot + (double)(SECONDSONWEEK * (gpsWn - wnt)) );

		// Compute factor w
		w = fmod ( gpsT - deltaUTC - (double)(HALFDAY), (double)(DAYSECS) ) +
				HALFDAY;

		// Compute GPS second in the UTC frame (in seconds of day).
		*utcT = fmod ( w, (double)(DAYSECS) + ddtlsf - ddtls );

		// Compute GPS second in the UTC frame (in seconds of week).
		*utcT = fmod (*utcT,(double)(DAYSECS))+((int)(gpsT - deltaUTC) / DAYSECS ) * DAYSECS;

	}// End of is case (b) is given

	else { // Case (c) is given

		// Compute deltaUTC
		deltaUTC = ddtlsf + a0 + a1 * ( gpsT - dtot + (double)(SECONDSONWEEK * (gpsWn - wnt)) );

		// Compute GPS second in the UTC frame (in seconds of week).
		*utcT = gpsT - deltaUTC;

	}// End of is case (c) is given

	// Check if a week roll-over eas produced.
	if (*utcT < 0.0) {
		*utcT = *utcT + SECONDSONWEEK;
		*utcWn = *utcWn-1;
	}

} // End of function GPStoUTC

