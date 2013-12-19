#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "azzaToRaDec.h"

/* test SCRAM az/za to ra/dec conversion routine */
/* function definition is: */ 
/* 

void scramAzZatoRaDec(long int AGC_SysTime, long int AGC_Time, double AGC_Az, double AGC_Za, \
	double AlfaMotorPosition, int Beam, int DoNotPrecess, double * ra, double * dec, \
	struct tm * hittime);

*/

/* AGC_SysTime, AGC_Time, AGC_Az, AGC_Za, AlfaMotorPosition are all taken from the SCRAM header */
/* Beam is the ALFA Beam number (0->6), 0 is center beam */

/* DoNotPrecess controls whether or not AO Phil's routine is commanded to precess back to J2000 */
/* DoNotPrecess = 0 returns J2000, DoNotPrecess = 1 returns epoch of the day */
/* ra, dec and hittime return populated with decimal hours, decimal degrees, AST time, respectively */

/* To Compile: */
/* cd to aocoord, type 'make' */
/* compile this code with: gcc -o azza_scram_test azza_scram_test.c -lm -I./aocoord -L./aocoord -lazzatoradec */
/* changing the path to 'aocoord' if necessary */




int main() {


/* return values for the conversion function */
struct tm asttime;
double ra;
double dec;


/* test values */
long AGC_SysTime = 1310902583; 
double AGC_Az = 126.624700; 
double AGC_Za = 12.432500;
long AGC_Time = 27383004;
double AlfaMotorPosition = 57.144062;
int Beam = 0;


/* CIMA log corresponding to above */
//2011-Jul-17 07:29:26 LOG4    log_alfabeam: ALFABEAM 0   2.076111111111117  25.414722222222217
//2011-Jul-17 07:29:26 LOG4    log_alfabeam: ALFABEAM 1   2.080337707241151  25.327191428718528
//2011-Jul-17 07:29:26 LOG4    log_alfabeam: ALFABEAM 2   2.083547486187920  25.407118095621545
//2011-Jul-17 07:29:26 LOG4    log_alfabeam: ALFABEAM 3   2.079374514047682  25.495161850111330
//2011-Jul-17 07:29:26 LOG4    log_alfabeam: ALFABEAM 4   2.071911656280776  25.502636621445280
//2011-Jul-17 07:29:26 LOG4    log_alfabeam: ALFABEAM 5   2.068652232439662  25.421407342929115
//2011-Jul-17 07:29:26 LOG4    log_alfabeam: ALFABEAM 6   2.072905154381599  25.334008876433305


/* loop counter */
int i;


/* test Atlantic Standard Time/Date return */
scramAzZatoRaDec(AGC_SysTime, AGC_Time, AGC_Az, AGC_Za, AlfaMotorPosition, Beam, 0,&ra, &dec, &asttime);

fprintf(stderr, "\ndate of observation: %d/%d/%d  time (AST): %02d:%02d:%02d\n\n", asttime.tm_mon+1, asttime.tm_mday,asttime.tm_year+1900,asttime.tm_hour, asttime.tm_min, asttime.tm_sec);

/* repeat conversion call for all seven beams, specifying J2000 coordinates */
fprintf(stderr, "\nEpoch J2000:\n");
for (i=0;i<7;i++) {
scramAzZatoRaDec(AGC_SysTime, AGC_Time, AGC_Az, AGC_Za, AlfaMotorPosition, i, 0, &ra, &dec, &asttime);
fprintf(stderr, "BEAM: %d RA: %f DEC: %f\n", i, ra, dec);
}

/* repeat conversion call for all seven beams, specifying EoD coordinates */
fprintf(stderr, "\nEpoch of the day:\n");
for (i=0;i<7;i++) {
scramAzZatoRaDec(AGC_SysTime, AGC_Time, AGC_Az, AGC_Za, AlfaMotorPosition, i, 1, &ra, &dec, &asttime);
fprintf(stderr, "BEAM: %d RA: %f DEC: %f\n", i, ra, dec);
}



}




