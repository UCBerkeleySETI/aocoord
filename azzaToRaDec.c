/*******************************************************************************
* azzaToRaDec - convert az,za to ra dec J2000
*
* DESCRIPTION
*
* Given az,za, mjd and utcFraction of day, compute the corresponding
* ra,dec J2000 coordinates.  Prior to calling this routine you must
* call azzaToRaDecInit once to initialize the AZZA_TO_RADEC_INFO structure.
* This struct holds information on the observatory position, going
* utc to ut1, pointing  model information, and how often to update the 
* precession, nutation matrix.
*
* RETURN
* The J2000 ra,dec corresponding to the requested az,za, and time. The
* data is returned as radians.
*
*NOTE:
* The azFd is the azimuth of the feed. This is 180 degrees from the azimuth
*of the source. For the gregorian dome this is the same as the encoder
*azimuth. For the carriage house, the azFd is encoderAzimuth-180.
*
*  The model correction is not yet implemented.
*/ 

/*	includes	*/
#include	<stdio.h>
#include	<string.h>
#include	<azzaToRaDec.h>
#include	<mathCon.h>
#include	<math.h>
#include 	<stdlib.h>  //added for sV.v


void  azzaToRaDec
      (
	double	        azFd,			/* feed azimuth*/
	double		    za,         /* encoder za*/
    int             mjd,        /* modified julian date*/
	double			utcFrac,	/* fraction of utc day*/
	int				ofDate, /* if not zero return coord of date*/
	AZZA_TO_RADEC_INFO *pinfo,
    double         *praJ2_rd,  /* ra  j2000 or ofDate, radians*/
    double         *pdecJ2_rd  /* dec j2000 or ofDate, radians*/
	) 
{
	double  azEl_V[3];
    double  ut1Frac;                 /* ut1 fraction of a day*/
	double	raDApp_V[3],haAppD_V[3];	/* 3 vectors*/
	double  raDTrue_V[3],raDJ2000_V[3];
	double	lst;			/* localmean sidereal time*/
	double	last;			/* local apparent sidereal time*/
	double	aber_V[3];					/* aberation vector*/
	/*
	 * change azFd,za to az/el 3 vectors
	*/
	anglesToVec3((azFd+180.)*C_DEG_TO_RAD,(90.-za)*C_DEG_TO_RAD,azEl_V);
	azElToHa_V(azEl_V,pinfo->obsPosI.latRd,haAppD_V);
	ut1Frac=utcFrac + utcToUt1(mjd,utcFrac,&pinfo->utcI);
	if (ut1Frac >= 1.) {
       	ut1Frac-=1.;
       	mjd++;
    }
    else if (ut1Frac < 0.) {
           ut1Frac+=1.;
           mjd--;
    }
    /*	
	 * go from mjdUt to local  mean sidereal time 
	*/
    lst=ut1ToLmst(mjd,ut1Frac,pinfo->obsPosI.wLongRd);
	last=lst + pinfo->precNutI.eqEquinox;	/* local apparent sidereal time*/
	/*
	 * go ha/dec to  ra/dec apparent
	*/
	haToRaDec_V(haAppD_V,last,raDApp_V);
	/*
 	 * compute the annual aberation subtract from curr ra/dec
     * and then renormalize vector
	*/
    aberAnnual_V(mjd,ut1Frac,aber_V);/* compute the aberration offset*/
	VV3D_Sub(raDApp_V,aber_V,raDTrue_V);
	V3D_Normalize(raDTrue_V,raDTrue_V);/* re normalize it */
	/*
 	 *if not coord OfDate request,  apply inverse nutation,prec to go to j2000
	*/
	if (!ofDate) {
		precNut(mjd,ut1Frac,FALSE,&pinfo->precNutI,raDTrue_V,raDJ2000_V);
	} else {
		memcpy(raDJ2000_V,raDTrue_V,sizeof(raDTrue_V));	/* return of date*/
	}
	/*
 	* convert vectors back to angles
	*/
	vec3ToAngles(raDJ2000_V,TRUE,praJ2_rd,pdecJ2_rd);	/* go back to angles*/

	return;
}

/* Begin additions for SERENDIPV.v AS 2012*/
void scramAzZatoRaDec(long int AGC_SysTime, long int AGC_Time, double AGC_Az, double AGC_Za, double AlfaMotorPosition, int Beam, int DoNotPrecess, double * ra, double * dec, struct tm * hittime) {

static const double AST_TO_UTC_HOURS=4;
static const double D2R=(3.1415926535897932384626/180.0);  /* decimal degrees to radians conversion */


//=======================================================
//double seti_ao_timeMS2unixtime(long timeMs, time_t now) {
//=======================================================

long timeMs = AGC_Time;
time_t now = AGC_SysTime;

double secs_past_midnight, encoder_read_secs_past_midnight,fnow;
struct tm * now_tm;

fnow=now;
fnow-=AST_TO_UTC_HOURS*3600;                 // translate UTC to AST
now=(time_t)floor(fnow);
now_tm = gmtime(&now);  
*hittime = *now_tm;

secs_past_midnight = now_tm->tm_hour*3600 + now_tm->tm_min*60 + now_tm->tm_sec;
encoder_read_secs_past_midnight = (timeMs * 0.001);

fnow -= secs_past_midnight;              // just back up to midnight

fnow += encoder_read_secs_past_midnight;     // now now as the unix time of the encoder reading
fnow += AST_TO_UTC_HOURS*3600;                // translate AST to UTC


double coord_unixtime = fnow;

double zen_corr_coeff[13] = {-57.55,-95.56,-4.13,141.69,677.51,-10.41,-7.71,-10.39,0.08,0.43,-0.62, 0.03,-0.36};
double az_corr_coeff[13] = {-37,-6.05,92.35,-731.21,-1013.97,-24.53,-11.19,9.18,106.04,3.02,-1.74,-3.46,1.29};

double array_angle[7] = {0.0, 120.0, 180.0, -120.0, -60.0, 0.0, 60.0};

    
//    seti_time temp_time((time_t)0,coord_unixtime); 
//    coord_time = temp_time; 

double Za        = AGC_Za;    // don't change the header since we want to keep it as read.
double Az        = AGC_Az;    // don't change the header since we want to keep it as read.
double Rotation  = AlfaMotorPosition;
//    telescope_id Tel = channel_to_receiverid[channel];

    //ZenAzCorrection(Beam, &Za, &Az, Rotation);       
//	fprintf(stderr, "TIME %ld %ld %lf\n", scram_agc.Time, scram_agc.SysTime, seti_ao_timeMS2unixtime(scram_agc.Time, scram_agc.SysTime));



double CosAz, SinAz, SinZen, Cos2Az, Sin2Az, Cos3Az, Sin3Az, Cos6Az, Sin6Az, SinZen2, \
       Cos3Imb, Sin3Imb, AzRadians, ZenRadians, ZenCorrection, AzCorrection;


AzRadians  = Az  * D2R;
ZenRadians = Za * D2R;

CosAz  = cos(AzRadians);
SinAz  = sin(AzRadians);
Cos2Az = cos(2.0 * AzRadians);
Sin2Az = sin(2.0 * AzRadians);
Cos3Az = cos(3.0 * AzRadians);
Sin3Az = sin(3.0 * AzRadians);
Cos6Az = cos(6.0 * AzRadians);
Sin6Az = sin(6.0 * AzRadians);

SinZen  = sin(ZenRadians);
SinZen2 = SinZen * SinZen;

Cos3Imb = sin(ZenRadians - 0.1596997627) * Cos3Az;	// 0.159... = 9.15 deg.  This is related
Sin3Imb = sin(ZenRadians - 0.1596997627) * Sin3Az;    // to the balance of dome and CH. (via Phil)

ZenCorrection = zen_corr_coeff[ 0]         	+ \
                zen_corr_coeff[ 1] * CosAz 	+ \
                zen_corr_coeff[ 2] * SinAz 	+ \
                zen_corr_coeff[ 3] * SinZen	+ \
				  zen_corr_coeff[ 4] * SinZen2	+ \
				  zen_corr_coeff[ 5] * Cos3Az	+ \
				  zen_corr_coeff[ 6] * Sin3Az	+ \
				  zen_corr_coeff[ 7] * Cos3Imb	+ \
				  zen_corr_coeff[ 8] * Sin3Imb	+ \
				  zen_corr_coeff[ 9] * Cos2Az	+ \
				  zen_corr_coeff[10] * Sin2Az	+ \
				  zen_corr_coeff[11] * Cos6Az	+ \
				  zen_corr_coeff[12] * Sin6Az;

AzCorrection  = az_corr_coeff[ 0]                + \
                az_corr_coeff[ 1] * CosAz        + \
                az_corr_coeff[ 2] * SinAz        + \
                az_corr_coeff[ 3] * SinZen       + \
                az_corr_coeff[ 4] * SinZen2      + \
                az_corr_coeff[ 5] * Cos3Az       + \
                az_corr_coeff[ 6] * Sin3Az       + \
                az_corr_coeff[ 7] * Cos3Imb      + \
                az_corr_coeff[ 8] * Sin3Imb      + \
                az_corr_coeff[ 9] * Cos2Az       + \
                az_corr_coeff[10] * Sin2Az       + \
                az_corr_coeff[11] * Cos6Az       + \
                az_corr_coeff[12] * Sin6Az;





// Now correct for the offsets of the beams in the arrays.
// Arrays are distributed on an ellipse.
//  compute the az,za offsets to add to the az, za coordinates
//  these are great circle.
//  1. With zero rotation angle,  the th generated is counter clockwise.
//     The offsets in the azalfaoff.. table were generated from using
//     this orientation of the angle.
//  2. The rotangl of the array is positive clockwise (since that is the
//     way the floor rotates when given a positive number). We use a 
//     minus sign since it is opposite of the angle used to generate the 
//     offset array above.
//  3. After computing the angle of a beam and projecting it onto the 
//     major, minor axis of the ellipse, the values are subtracted 
//     from pixel 0 az,za. In other words, the positive direction is the
//     same as we've already used above for our correction, so we can just
//     add it to our correction.

  double array_az_ellipse;
  double array_za_ellipse;

  if(Beam != 0) {  
	  array_az_ellipse = 329.06;  //taken from seti_boinc/db/tools/s4_receivers.xml
	  array_za_ellipse = 384.005;
   } else {
	  array_az_ellipse = 0;
	  array_za_ellipse = 0;   
   }
  	

  double posrot=(array_angle[Beam] - Rotation)*D2R;
  AzCorrection+=array_az_ellipse*cos(posrot);
  ZenCorrection+=array_za_ellipse*sin(posrot);

  Za -= ZenCorrection / 3600.0;	  	        // Correction is in arcsec.
  Az  -= (AzCorrection  / 3600.0) / sin(Za * D2R); 	
  
  // Correction is in arcsec.

  // Sometimes the azimuth is not where the telescope is
  // pointing, but rather where it is physically positioned.
  Az    += 180.0; //taken from seti_boinc/db/tools/s4_receivers.xml

  if(Za < 0.0)                                // Have we corrected zen over
    {                                           //  the zenith?
    Za = 0.0 - Za;                          // If so, correct zen and
    Az += 180.0;                               //   swing azimuth around to
    }

  Az = wrap(Az, 0, 360, 1);              // az to 0 - 360 degrees


//double ra, dec;
seti_AzZaToRaDec(Az, Za, coord_unixtime, ra, dec, DoNotPrecess);

}




//__________________________________________________________________________
double wrap(double Val, long LowBound, long HighBound, int LowInclusive) {
//__________________________________________________________________________

	if(LowInclusive) {
		while(Val >= HighBound) Val -= HighBound; 
		while(Val <   LowBound) Val += HighBound; 
	} else {
		while(Val >  HighBound) Val -= HighBound; 
                while(Val <=  LowBound) Val += HighBound;
	}
	return Val;
}



//=======================================================
void seti_AzZaToRaDec(double Az, double Za, double coord_time, double * Ra, double * Dec, int DoNotPrecess) {
//=======================================================
// This calls AO Phil's code.
// Any desired model correction must be done prior to calling this routine.

    AZZA_TO_RADEC_INFO  azzaI; /* info for aza to ra dec*/

    const double d2r =  0.017453292;
    int dayNum, i_mjd;
    //DoNotPrecess - tells azzaToRaDec() to return coords of the date, ie not precessed
    double utcFrac;
    struct tm * coord_tm;
    time_t lcoord_time=(time_t)floor(coord_time);
    double fcoord_time=coord_time-floor(coord_time);

    coord_tm = gmtime(&lcoord_time);                 

    // arithmetic needed by AO functions
    coord_tm->tm_mon += 1;
    coord_tm->tm_year += 1900;
    dayNum = dmToDayNo(coord_tm->tm_mday,coord_tm->tm_mon,coord_tm->tm_year);
    i_mjd=gregToMjd(coord_tm->tm_mday, coord_tm->tm_mon, coord_tm->tm_year);
    utcFrac=(coord_tm->tm_hour*3600 + coord_tm->tm_min*60 + coord_tm->tm_sec + fcoord_time)/86400.;
    if (utcFrac >= 1.) {
        i_mjd++;
        utcFrac-=1.;
    }

    // call the AO functions
    if (azzaToRaDecInit(dayNum, coord_tm->tm_year, &azzaI) == -1) exit(1);
    // subtract 180 from Az because Phil wants dome azimuth
    azzaToRaDec(Az-180.0, Za, i_mjd, utcFrac, DoNotPrecess, &azzaI, Ra, Dec);

    // Ra to hours and Dec to degrees
    *Ra = (*Ra/d2r) / 15;
    *Dec = *Dec/d2r;

    // Take care of wrap situations in RA
    while (*Ra < 0) *Ra += 24;
    *Ra = fmod(*Ra,24);
}
