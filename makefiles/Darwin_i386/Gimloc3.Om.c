/* ../../src/navigation/Gimloc3.Om.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    integer incmax[2];
    real elvmax[2], scnmax[2], elvinc[2], scninc[2], elvln[2], scnpx[2];
} instco_;

#define instco_1 instco_

struct {
    doublereal xs[3], bt[9]	/* was [3][3] */, q3, pitch, roll, yaw;
    real pma, rma;
} elcomm_;

#define elcomm_1 elcomm_

/* Table of constant values */

static integer c__62 = 62;
static integer c__117 = 117;
static integer c__172 = 172;
static integer c__227 = 227;
static integer c__282 = 282;

/* *********************************************************************** */
/* ** */
/* **   INTEGRAL SYSTEMS, INC. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   ROUTINE   : SETCONS */
/* **   SOURCE    : F.SETCONS */
/* **   LOAD NAME : ANY */
/* **   PROGRAMMER: IGOR LEVINE */
/* ** */
/* **   VER.    DATA    BY   COMMENT */
/* **   ----  --------  ---  --------------------------------------------- */
/* **   A     02/16/89  IL   INITIAL CREATION */
/* ** */
/* **   B     05/19/94  NP   ADDED CALCULATION OF INSTRUMENT ELEVATION AND */
/* **                        SCAN ANGLE BIASES BASED ON USER INPUT */
/* *********************************************************************** */
/* ** */
/* **   THIS SUBROUTINE GENERATES CONSTANTS IN COMMON  INSTCOMM */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   CALLED BY       : ANY */
/* **   COMMONS MODIFIED: INSTCO */
/* **   INPUTS          : NONE */
/* **   OUTPUTS         : NONE */
/* **   ROUTINES CALLED : NONE */
/* ** */
/* *********************************************************************** */
/* *********************************************************************** */
/* Subroutine */ int setcon_(integer *instr, integer *ns_nad_cy__, integer *
	ns_nad_inc__, integer *ew_nad_cy__, integer *ew_nad_inc__)
{

/*     CALLING PARAMETERS */


/*     LOCAL VARIABLES */


/*     INCLUDE FILES */

/* *********************************************************************** */
/* ** */
/* **   INTEGRAL SYSTEMS, INC. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   NAME      : ELCONS */
/* **   TYPE      : DATA AREA */
/* **   SOURCE    : I.ELCONS */
/* ** */
/* **   VER.    DATA    BY                COMMENT */
/* **   ----  --------  ----------------  -------------------------------- */
/* **   A     01/09/89  I. LEVINE         INITIAL CREATION */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   DESCRIPTION */
/* **   MATHEMATICAL AND EARTH-RELATED CONSTANTS */
/* ** */
/* *********************************************************************** */
/* *********************************************************************** */

/*                    DEGREES TO RADIANS CONVERSION PI/180 */
/*                    NOMINAL RADIAL DISTANCE OF SATELLITE (km) */
/*                    EARTH EQUATORIAL RADIUS (km) */
/*                    EARTH FLATTENING COEFFICIENT = 1-(BE/AE) */
/* *********************************************************************** */
/* *********************************************************************** */
/* *********************************************************************** */
/* ** */
/* **   INTEGRAL SYSTEMS, INC. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   NAME      : INSTCOMM */
/* **   TYPE      : DATA AREA */
/* **   SOURCE    : I.INSTCOMM */
/* ** */
/* **   VER.    DATA    BY                COMMENT */
/* **   ----  --------  ----------------  -------------------------------- */
/* **   A     02/16/89  I. LEVINE         INITIAL CREATION */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   DESCRIPTION */
/* **   COMMON AREA FOR INSTRUMENT-RELATED CONTROL PARAMETERS */
/* ** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     VARIABLES */
/*        CONSTANTS NEEDED TO PERFORM TRANSFORMATIONS BETWEEN THE */
/*        LATITUDE/LONGITUDE, LINE/PIXEL AND INSTRUMENT CYCLES/INCREMENTS */
/*        COORDINATES. */

/*                       NUMBER OF INCREMENTS PER CYCLE */
/*                       BOUNDS IN ELEVATION (RADIANS) */
/*                       BOUNDS IN SCAN ANGLE (RADIANS) */
/*                       CHANGE IN ELEVATION ANGLE PER INCREMENT (RAD) */
/*                       CHANGE IN SCAN ANGLE PER INCREMENT (RADIANS) */
/*                       ELEVATION ANGLE PER DETECTOR LINE (RADIANS) */
/*                       SCAN ANGLE PER PIXEL (RADIANS) */
    instco_1.incmax[0] = 6136;
    instco_1.incmax[1] = 2805;
    instco_1.elvinc[0] = 8e-6f;
    instco_1.elvinc[1] = 1.75e-5f;
    instco_1.scninc[0] = 1.6e-5f;
    instco_1.scninc[1] = 3.5e-5f;
    instco_1.elvln[0] = 2.8e-5f;
    instco_1.elvln[1] = 2.8e-4f;
    instco_1.scnpx[0] = 1.6e-5f;
    instco_1.scnpx[1] = 2.8e-4f;
/*     ************************************************************ */
/*     COMMENTED OUT ELEVATION AND SCAN BIAS CONSTANTS SINCE INSTRUMENT */
/*     EARTH NADIR POSITION IS AVAILABLE IN GVAR DATA AND PERIODICALLY */
/*     UPDATED */


    instco_1.elvmax[0] = .220896f;
    instco_1.elvmax[1] = .22089375f;
    instco_1.scnmax[0] = .24544f;
    instco_1.scnmax[1] = .2454375f;
/*     RECOMPUTE ELEVATION AND SCAN BIASES BASED ON USER INPUTS OF */
/*     CYCLES & INCREMENTS OBTAINED FROM GVAR */
/*     ELVMAX(INSTR) = (NS_NAD_CY*INCMAX(INSTR)+NS_NAD_INC)*ELVINC(INSTR) */
/*      IF(INSTR.EQ.1)THEN */
/*        ELVMAX(INSTR)=(NS_NAD_CY*INCMAX(INSTR)+NS_NAD_INC)*ELVINC(INSTR) */
/*      ELSE */
/*        ELVMAX(INSTR)=((9-NS_NAD_CY)*INCMAX(INSTR)-NS_NAD_INC) */
/*     +                *ELVINC(INSTR) */
/*      ENDIF */
/*      SCNMAX(INSTR) = (EW_NAD_CY*INCMAX(INSTR)+EW_NAD_INC)*SCNINC(INSTR) */
/*     ************************************************************ */
    return 0;
} /* setcon_ */

/* *********************************************************************** */
/* ** */
/* **   INTEGRAL SYSTEMS, INC. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   ROUTINE   : LMODEL */
/* **   SOURCE    : F.LMODEL */
/* **   LOAD NAME : ANY */
/* **   PROGRAMMER: IGOR LEVINE */
/* ** */
/* **   VER.    DATA    BY   COMMENT */
/* **   ----  --------  ---  --------------------------------------------- */
/* **   1     01/09/89  IL   INITIAL CREATION */
/* **   2     06/02/89  IL   COORDINATE AXES CHANGED ACCORDING TO */
/* **                        FORD'S DEFINITION IN SDAIP, DRL 504-01 */
/* **   3     08/21/89  IL   CORRECTED ORBIT ANGLE COMPUTATIONS */
/* **   4     03/08/94  SC   S/C COMPENSATION APPLIED UNCONDITIONALLY; */
/* **                        REFERENCE RADIAL DISTANCE, LATITZDE AND */
/* **                        ORBIT YAW SET TO ZERO IF IMC DISABLED. */
/* **   5     03/08/94  SC   ADDED TRAP FOR SLAT=SYAW=0; CORRECTED */
/* **                        EXPRESSION FOR LAM. */
/* *********************************************************************** */
/* ** */
/* **   THIS SUBROUTINE COMPUTES THE POSITION OF THE SATELLITE AND THE */
/* **   ATTITZDE OF THE IMAGER OR SOUNDER.  THE CALCULATIONS ARE BASED */
/* **   ON THE OATS ORBIT AND ATTITZDE MODEL REPRESENTED BY THE O&A */
/* **   PARAMETER SET IN GVAR BLOCK 0. */
/* **        INPUTS: */
/* **          TIME,  EPOCH TIME, O&A PARAMETER SET, IMC STATUS. */
/* ** */
/* **        OUTPUTS: */
/* **          THE SPACECRAFT POSITION VECTOR IN EARTH FIXED COORDINATES; */
/* **          THE GEOMETRIC ROLL, PITCH, YAW ANGLES AND THE ROLL, */
/* **          PITCH MISALIGNMENTS FOR EITHER THE IMAGER OR THE SOUNDER; */
/* **          THE EARTH FIXED TO INSTRUMENT FRAME TRANSFORMATION MATRIX; */
/* **          GEOGRAPHIC LATITZDE AND LONGITZDE AT SUBSATELLITE POINT. */
/* ** */
/* **   DESCRIPTION */
/* **   LMODEL ACCEPTS AN INPUT DOUBLE PRECISION TIME IN MINUTES FROM */
/* **   1950, JAN.1.0 AND AN INPUT SET OF O&A PARAMETERS AND COMPUTES */
/* **   POSITION OF THE SATELLITE, THE ATTITZDE ANGLES AND ATTITZDE */
/* **   MISALIGNMENTS AND THE INSTRUMENT TO EARTH FIXED COORDINATES */
/* **   TRANSFORMATION MATRIX. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   CALLED BY       : ANY */
/* **   COMMONS MODIFIED: /ELCOMM/ XS,Q3,PITCH,ROLL,YAW,PMA,RMA,BT */
/* **   INPUTS          : NONE */
/* **   OUTPUTS         : NONE */
/* **   ROUTINES CALLED : INST2ER,GATT */
/* ** */
/* *********************************************************************** */
/* *********************************************************************** */
/* Subroutine */ int lmodel_(doublereal *t, doublereal *tu, doublereal *rec, 
	integer *imc, doublereal *rlat, doublereal *rlon)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal), sqrt(doublereal), atan2(
	    doublereal, doublereal), tan(doublereal), atan(doublereal);

    /* Local variables */
    static doublereal b[9]	/* was [3][3] */, r__, u, w, ca, ci, sa, dr, 
	    cu, te, wa, cw, si, ts, su, sw, cw1, c2w, cw3, sw1, s2w, sw3, asc,
	     lam, phi, psi, dlat;
    extern doublereal gatt_(integer *, doublereal *, doublereal *, doublereal 
	    *);
    static doublereal slat, dyaw, syaw;
    extern /* Subroutine */ int inst2er_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);


/*     CALLING ARGUMENTS */

/*                                   INPUT TIME FROM 1950, JAN 1.0 (MIN) */
/*                                   EPOCH TIME FROM 1950, JAN 1.0 (MIN) */
/*                                   INPUT O&A PARAMETER SET */
/*                                   INPUT IMC STATUS: 0 - ON, 1 - OFF */
/*                                   SUBSATELLITE GEODETIC LATITZDE (RAD) */
/*                                   SUBSATELLITE LONGITZDE IN RADIANS */

/*     LOCAL VARIABLES */

/*                    NORMALIZED SATELLITE DISTANCE (IN UNITS OF KMER9) */
/*                    TIME FROM EPOCH IN MINUTES */
/*                    SPACCRAFT TO EARTH FIXED COORDINATES TRANSFORMATION */
/*                    MATRIX */
/*                    EXPONENENTIAL TIME DELAY FROM EPOCH (IN MINUTES) */
/*                    SUBSATELLITE GEOCENTRIC LATITZDE IN RADIANS */
/*                    RADIAL DISTANCE FROM THE NOMINAL (KM) */
/*                    ORBITAL YAW (IN RADIANS) */
/*                    IMC LONGITZDE (IN RADIANS) */
/*                    ARGUMENT OF LATITZDE (IN RADIANS) */
/*                    SIN(U), COS(U) */
/*                    SINE AND COSINE OF THE ORBIT INCLINATION */
/*                    SINE OF GEOCENTRIC LATITZDE */
/*                    LONGITZDE OF THE ASCENDING NODE IN RADIANS */
/*                    SINE AND COSINE OF ASC */
/*                    SINE OF THE ORBIT YAW */
/*                    SOLAR ORBIT ANGLE IN RADIANS */
/*                    ORBIT ANGLE IN RADIANS */
/*                    SIN(W),  COS(W) */
/*                    SIN(2*W),  COS(2*W) */
/*                    SIN(0.927*W),  COS(0.927*W) */
/*                    SINE AND COSINE OF 1.9268*W */
/*                    CHANGE IN SINE OF GEOCENTRIC LATITZDE */
/*                    CHANGE IN SINE OF ORBIT YAW */
/*                    SUBROUTINE FUNCTION */
/*      REAL*8 A1,A2 */
/*                    WORK AREAS */

/*     INCLUDE FILES */

/* *********************************************************************** */
/* ** */
/* **   INTEGRAL SYSTEMS, INC. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   NAME      : ELCONS */
/* **   TYPE      : DATA AREA */
/* **   SOURCE    : I.ELCONS */
/* ** */
/* **   VER.    DATA    BY                COMMENT */
/* **   ----  --------  ----------------  -------------------------------- */
/* **   A     01/09/89  I. LEVINE         INITIAL CREATION */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   DESCRIPTION */
/* **   MATHEMATICAL AND EARTH-RELATED CONSTANTS */
/* ** */
/* *********************************************************************** */
/* *********************************************************************** */

/*                    DEGREES TO RADIANS CONVERSION PI/180 */
/*                    NOMINAL RADIAL DISTANCE OF SATELLITE (km) */
/*                    EARTH EQUATORIAL RADIUS (km) */
/*                    EARTH FLATTENING COEFFICIENT = 1-(BE/AE) */

/* *********************************************************************** */

/*     ASSIGN REFERENCE VALUES TO THE SUBSATELLITE LONGITZDE AND */
/*     LATITZDE, THE RADIAL DISTANCE AND THE ORBIT YAW. */

/* *********************************************************************** */
/* *********************************************************************** */
/* **   INTEGRAL SYSTEMS, INC. */
/* *********************************************************************** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   NAME      : ELCOMM */
/* **   TYPE      : DATA AREA */
/* **   SOURCE    : I.ELCOMM */
/* ** */
/* **   VER.    DATA    BY                COMMENT */
/* **   ----  --------  ----------------  -------------------------------- */
/* **   A     01/09/89  I. LEVINE         INITIAL CREATION */
/* *********************************************************************** */
/* **   DESCRIPTION */
/* **   INSTRUMENT POSITION AND ATTITUDE VARIABLES AND TRANSFORMATION */
/* **   MATRIX */
/* *********************************************************************** */

/*     COMMON VARIABLES */

/*                      NORMALIZED S/C POSITION IN ECEF COORDINATES */
/*                      ECEF TO INSTRUMENT COORDINATES TRANSFORMATION */
/*                      USED IN SUBROUTINE LPOINT */
/*                          PITCH,ROLL,YAW ANGLES OF INSTRUMENT (RAD) */
/*                          PITCH,ROLL MISALIGNMENTS OF INSTRUMENT (RAD) */
    /* Parameter adjustments */
    --rec;

    /* Function Body */
    lam = rec[5];
    dr = rec[6];
    phi = rec[7];
    psi = rec[8];

/*     ASSIGN REFERENCE VALUES TO THE ATTITZDES AND MISALIGNMENTS */

    elcomm_1.roll = rec[9];
    elcomm_1.pitch = rec[10];
    elcomm_1.yaw = rec[11];
    elcomm_1.rma = 0.f;
    elcomm_1.pma = 0.f;

/*     IF IMC IS OFF, COMPUTE CHANGES IN THE SATELLITE ORBIT */

    if (*imc != 0) {

/*     SET REFERENCE RADIAL DISTANCE, LATITZDE AND ORBIT YAW TO ZERO */

	dr = 0.;
	phi = 0.;
	psi = 0.;

/*     COMPUTE TIME SINCE EPOCH (IN MINUTES) */

	ts = *t - *tu;

/*     COMPUTES ORBIT ANGLE AND THE RELATED TRIGONOMETRIC FUNCTIONS. */
/*     EARTH ROTATIONAL RATE=.729115E-4 (RAD/S) */

	w = ts * .0043746900000000005;
	sw = sin(w);
	cw = cos(w);
	sw1 = sin(w * .927);
	cw1 = cos(w * .927);
	s2w = sin(w * 2.);
	c2w = cos(w * 2.);
	sw3 = sin(w * 1.9268);
	cw3 = cos(w * 1.9268);

/*     COMPUTES CHANGE IN THE IMC LONGITZDE FROM THE REFERENCE */

	lam = lam + rec[18] + (rec[19] + rec[20] * w) * w + (rec[27] * sw1 + 
		rec[28] * cw1 + rec[21] * sw + rec[22] * cw + rec[23] * s2w + 
		rec[24] * c2w + rec[25] * sw3 + rec[26] * cw3 + w * (rec[29] *
		 sw + rec[30] * cw)) * 2.;

/*     COMPUTES CHANGE IN RADIAL DISTANCE FROM THE REFERENCE (KM) */

	dr = dr + rec[31] + rec[32] * cw + rec[33] * sw + rec[34] * c2w + rec[
		35] * s2w + rec[36] * cw3 + rec[37] * sw3 + rec[38] * cw1 + 
		rec[39] * sw1 + w * (rec[40] * cw + rec[41] * sw);

/*     COMPUTES THE SINE OF THE CHANGE IN THE GEOCENTRIC LATITZDE */

	dlat = rec[42] + rec[43] * cw + rec[44] * sw + rec[45] * c2w + rec[46]
		 * s2w + w * (rec[47] * cw + rec[48] * sw) + rec[49] * cw1 + 
		rec[50] * sw1;

/*     COMPUTES GEOCENTRIC LATITZDE BY USING AN EXPANSION FOR ARCSINE */

	phi += dlat * (dlat * dlat / 6. + 1.);

/*     COMPUTES SINE OF THE CHANGE IN THE ORBIT YAW */

	dyaw = rec[51] + rec[52] * sw + rec[53] * cw + rec[54] * s2w + rec[55]
		 * c2w + w * (rec[56] * sw + rec[57] * cw) + rec[58] * sw1 + 
		rec[59] * cw1;

/*     COMPUTES THE ORBIT YAW BY USING AN EXPANSION FOR ARCSINE. */

	psi += dyaw * (dyaw * dyaw / 6. + 1.);

/*     CALCULATION OF CHANGES IN THE SATELLITE ORBIT ENDS HERE */

    }

/*     CONVERSION OF THE IMC LONGITZDE AND ORBIT YAW TO THE SUBSATELLITE */
/*     LONGITZDE AND THE ORBIT INCLINATION (REF: GOES-PCC-TM-2473, INPUTS */
/*     REQUIRED FOR EARTH LOCATION AND GRIDDING BY SPS,  JUNE 6, 1988) */

    slat = sin(phi);
    syaw = sin(psi);
/* Computing 2nd power */
    d__1 = slat;
/* Computing 2nd power */
    d__2 = syaw;
    si = d__1 * d__1 + d__2 * d__2;
    ci = sqrt(1. - si);
    si = sqrt(si);
    if (slat == 0. && syaw == 0.) {
	u = 0.;
    } else {
	u = atan2(slat, syaw);
    }
    su = sin(u);
    cu = cos(u);

/*     COMPUTES LONGITZDE OF THE ASCENDING NODE */

    asc = lam - u;
    sa = sin(asc);
    ca = cos(asc);

/*     COMPUTES THE SUBSATELLITE GEOGRAPHIC LATITZDE */

    *rlat = atan(tan(phi) * 1.0067391845079681f);

/*     COMPUTES THE SUBSATELLITE LONGITZDE */

    *rlon = asc + atan2(ci * su, cu);

/*     COMPUTES THE SPACECRAFT TO EARTH FIXED COORDINATES TRANSFORMATION */
/*     MATRIX: */
/*         (VECTOR IN ECEF COORDINATES) = B * (VECTOR IN S/C COORDINATES) */

    b[3] = -sa * si;
    b[4] = ca * si;
    b[5] = -ci;
    b[6] = -ca * cu + sa * su * ci;
    b[7] = -sa * cu - ca * su * ci;
    b[8] = -slat;
    b[0] = -ca * su - sa * cu * ci;
    b[1] = -sa * su + ca * cu * ci;
    b[2] = cu * si;

/*     COMPUTES THE NORMALIZED SPACECRAFT POSITION VECTOR IN EARTH FIXED */
/*     COORDINATES - XS. */

    r__ = (dr + 42164.365) / 6378.137;
    elcomm_1.xs[0] = -b[6] * r__;
    elcomm_1.xs[1] = -b[7] * r__;
    elcomm_1.xs[2] = -b[8] * r__;

/*     PRECOMPUTES Q3 (USED IN LPOINT) */

/* Computing 2nd power */
    d__1 = elcomm_1.xs[0];
/* Computing 2nd power */
    d__2 = elcomm_1.xs[1];
/* Computing 2nd power */
    d__3 = elcomm_1.xs[2];
    elcomm_1.q3 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 * 
	    1.0067391845079681f - 1.;

/*     COMPUTES THE ATTITZDES AND MISALIGNMENTS IF IMC IS OFF */

    if (*imc != 0) {

/*     COMPUTES THE SOLAR ORBIT ANGLE */

	wa = rec[60] * ts;

/*     COMPUTES THE DIFFERENCE BETWEEN CURRENT TIME, TS, AND THE */
/*     EXPONENTIAL TIME, REC(61). NOTE THAT BOTH TIMES ARE SINCE EPOCH. */

	te = ts - rec[61];

/*     COMPUTES ROLL + ROLL MISALIGNMENT */

	elcomm_1.roll += gatt_(&c__62, &rec[1], &wa, &te);

/*     COMPUTES PITCH + PITCH MISALIGNMENT */

	elcomm_1.pitch += gatt_(&c__117, &rec[1], &wa, &te);

/*     COMPUTES YAW */

	elcomm_1.yaw += gatt_(&c__172, &rec[1], &wa, &te);

/*     COMPUTES ROLL MISALIGNMENT */

	elcomm_1.rma = gatt_(&c__227, &rec[1], &wa, &te);

/*     COMPUTES PITCH MISALIGNMENT */

	elcomm_1.pma = gatt_(&c__282, &rec[1], &wa, &te);

/*     APPLY THE SPCECRAFT COMPENSATION */

	elcomm_1.roll += rec[15];
	elcomm_1.pitch += rec[16];
	elcomm_1.yaw += rec[17];
    }

/*     COMPUTES THE INSTRUMENT TO EARTH FIXED COORDINATES TRANSFORMATION */
/*     MATRIX - BT */

    inst2er_(&elcomm_1.roll, &elcomm_1.pitch, &elcomm_1.yaw, b, elcomm_1.bt);
    return 0;
} /* lmodel_ */

/* *********************************************************************** */
/* *********************************************************************** */
/* ** */
/* **   INTEGRAL SYSTEMS, INC. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   ROUTINE   : GPOINT */
/* **   SOURCE    : F.GPOINT */
/* **   LOAD NAME : ANY */
/* **   PROGRAMMER: IGOR LEVINE */
/* ** */
/* **   VER.    DATA    BY   COMMENT */
/* **   ----  --------  ---  --------------------------------------------- */
/* **   A     12/10/87  IL   INITIAL CREATION */
/* **   A     06/10/88  IL   REPLACED ASIN WITH ATAN TO SAVE TIME */
/* **   A     06/02/89  IL   COORDINATE AXES CHANGED ACCORDING TO */
/* **                        FORD'S DEFINITION IN SDAIP, DRL 504-01 */
/* **   4     03/08/94  SC   IMPLEMENTED NEW FORMULAE FOR SCAN ANGLE */
/* **                        CORRECTION DUE TO MISALIGNMENTS */
/* *********************************************************************** */
/* ** */
/* **   THIS SUBROUTINE CONVERTS GEOGRAPHIC LATITZDE AND LONGITZDE */
/* **   TO THE RELATED ELEVATION AND SCAN ANGLES. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   CALLED BY       : ANY */
/* **   COMMONS MODIFIED: NONE */
/* **   INPUTS          : NONE */
/* **   OUTPUTS         : NONE */
/* **   ROUTINES CALLED : NONE */
/* ** */
/* *********************************************************************** */
/* *********************************************************************** */
/* Subroutine */ int gpoint_(doublereal *rlat, doublereal *rlon, doublereal *
	alf, doublereal *gam, integer *ierr)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sin(doublereal), sqrt(doublereal), cos(doublereal), atan(
	    doublereal), tan(doublereal);

    /* Local variables */
    static doublereal f[3], u[3], w1, w2, ft[3], sing, slat;


/*     CALLING PARAMETERS */

/*                             GEOGRAPHIC LATITZDE IN RADIANS (INPUT) */
/*                             GEOGRAPHIC LONGITZDE IN RADIANS (INPUT) */
/*                             ELEVATION ANGLE IN RADIANS (OUTPUT) */
/*                             SCAN ANGLE IN RADIANS (OUTPUT) */
/*                             OUTPUT STATUS; 0 - SUCCESSFUL COMPLETION, */
/*                             1 - POINT WITH GIVEN LAT/LON IS INVISIBLE */

/*     LOCAL VARIABLES */

/*                         POINTING VECTOR IN EARTH CENTERED COORDINATES */
/*                         POINTING VECTOR IN INSTRUMENT COORDINATES */
/*                         COORDINATES OF THE EARTH POINT (KM) */
/*                                    WORK SPACE */

/*     INCLUDE FILES */

/* *********************************************************************** */
/* ** */
/* **   INTEGRAL SYSTEMS, INC. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   NAME      : ELCONS */
/* **   TYPE      : DATA AREA */
/* **   SOURCE    : I.ELCONS */
/* ** */
/* **   VER.    DATA    BY                COMMENT */
/* **   ----  --------  ----------------  -------------------------------- */
/* **   A     01/09/89  I. LEVINE         INITIAL CREATION */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   DESCRIPTION */
/* **   MATHEMATICAL AND EARTH-RELATED CONSTANTS */
/* ** */
/* *********************************************************************** */
/* *********************************************************************** */

/*                    DEGREES TO RADIANS CONVERSION PI/180 */
/*                    NOMINAL RADIAL DISTANCE OF SATELLITE (km) */
/*                    EARTH EQUATORIAL RADIUS (km) */
/*                    EARTH FLATTENING COEFFICIENT = 1-(BE/AE) */
/* *********************************************************************** */

/*     COMPUTES SINUS OF GEOGRAPHIC (GEODETIC) LATITZDE */

/* *********************************************************************** */
/* *********************************************************************** */
/* **   INTEGRAL SYSTEMS, INC. */
/* *********************************************************************** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   NAME      : ELCOMM */
/* **   TYPE      : DATA AREA */
/* **   SOURCE    : I.ELCOMM */
/* ** */
/* **   VER.    DATA    BY                COMMENT */
/* **   ----  --------  ----------------  -------------------------------- */
/* **   A     01/09/89  I. LEVINE         INITIAL CREATION */
/* *********************************************************************** */
/* **   DESCRIPTION */
/* **   INSTRUMENT POSITION AND ATTITUDE VARIABLES AND TRANSFORMATION */
/* **   MATRIX */
/* *********************************************************************** */

/*     COMMON VARIABLES */

/*                      NORMALIZED S/C POSITION IN ECEF COORDINATES */
/*                      ECEF TO INSTRUMENT COORDINATES TRANSFORMATION */
/*                      USED IN SUBROUTINE LPOINT */
/*                          PITCH,ROLL,YAW ANGLES OF INSTRUMENT (RAD) */
/*                          PITCH,ROLL MISALIGNMENTS OF INSTRUMENT (RAD) */
    sing = sin(*rlat);
    w1 = sing * -.013343333245450451f * sing;

/*     SINUS OF THE GEOCENTRIC LATITZDE */

    slat = ((w1 * .375 - .5) * w1 + 1.) * sing / 1.0067391845079681f;

/*     COMPUTES LOCAL EARTH RADIUS AT SPECIFIED POINT */

    w2 = slat * slat;
    w1 = w2 * .0067391845079680657f;
    w1 = (w1 * .375 - .5) * w1 + 1.;

/*     COMPUTES CARTESIAN COORDINATES OF THE POINT */

    u[2] = slat * w1;
    w2 = w1 * sqrt(1. - w2);
    u[0] = w2 * cos(*rlon);
    u[1] = w2 * sin(*rlon);

/*     POINTING VECTOR FROM SATELLITE TO THE EARTH POINT */

    f[0] = u[0] - elcomm_1.xs[0];
    f[1] = u[1] - elcomm_1.xs[1];
    f[2] = u[2] - elcomm_1.xs[2];
    w2 = u[0] * (real) f[0] + u[1] * (real) f[1] + u[2] * (real) f[2] * 
	    1.0067391845079681f;

/*     VERIFIES VISIBILITY OF THE POINT */

    if (w2 > 0.) {
/*                               INVISIBLE POINT ON THE EARTH */
	*ierr = 1;
	*alf = 99999.;
	*gam = 99999.;
	return 0;
    }

/*     CONVERTS POINTING VECTOR TO INSTRUMENT COORDINATES */

    ft[0] = elcomm_1.bt[0] * f[0] + elcomm_1.bt[1] * f[1] + elcomm_1.bt[2] * 
	    f[2];
    ft[1] = elcomm_1.bt[3] * f[0] + elcomm_1.bt[4] * f[1] + elcomm_1.bt[5] * 
	    f[2];
    ft[2] = elcomm_1.bt[6] * f[0] + elcomm_1.bt[7] * f[1] + elcomm_1.bt[8] * 
	    f[2];

/*     CONVERTS POINTING VECTOR TO SCAN AND ELEVATION ANGLES AND */
/*     CORRECTS FOR THE ROLL AND PITCH MISALIGNMENTS */

/* Computing 2nd power */
    d__1 = ft[1];
/* Computing 2nd power */
    d__2 = ft[2];
    *gam = atan(ft[0] / sqrt(d__1 * d__1 + d__2 * d__2));
    *alf = -atan(ft[1] / ft[2]);
    w1 = sin(*alf);
    w2 = cos(*gam);
    *alf = *alf + elcomm_1.rma * (1. - cos(*alf) / w2) + elcomm_1.pma * w1 * (
	    1. / w2 + tan(*gam));
    *gam -= elcomm_1.rma * w1;
    *ierr = 0;
    return 0;
} /* gpoint_ */

/* *********************************************************************** */
/* *********************************************************************** */
/* ** */
/* **   INTEGRAL SYSTEMS, INC. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   ROUTINE   : INST2ER */
/* **   SOURCE    : F.INST2ER */
/* **   LOAD NAME : ANY */
/* **   PROGRAMMER: IGOR LEVINE */
/* ** */
/* **   VER.    DATA    BY   COMMENT */
/* **   ----  --------  ---  --------------------------------------------- */
/* **   1     08/16/88  IL   INITIAL CREATION */
/* **   2     11/11/88  IL   TRIGONOMETRIC FUNCTIONS REPLACED WITH */
/* **                        SMALL ANGLE APPROXIMATIONS */
/* **   3    06/02/89   IL   COORDINATE AXES CHANGED ACCORDING TO */
/* **                        FORD'S DEFINITION IN SDAIP, DRL 504-01 */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   INST2ER ACCEPTS THE SINGLE PRECISION ROLL, PITCH AND YAW ANGLES */
/* **   OF AN INSTRUMENT AND RETURNS THE DOUBLE PRECISION INSTRUMENT TO */
/* **   EARTH COORDINATES TRANSFORMATION MATRIX. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   CALLED BY       : ANY */
/* **   COMMONS MODIFIED: NONE */
/* **   INPUTS          : NONE */
/* **   OUTPUTS         : NONE */
/* **   ROUTINES CALLED : NONE */
/* ** */
/* *********************************************************************** */
/* *********************************************************************** */
/* Subroutine */ int inst2er_(doublereal *r__, doublereal *p, doublereal *y, 
	doublereal *a, doublereal *at)
{
    static integer i__, j;
    static doublereal rpy[9]	/* was [3][3] */;


/*     CALLING PARAMETERS */

/*                                   ROLL ANGLE IN RADIANS */
/*                                   PITCH ANGLE IN RADIANS */
/*                                   YAW ANGLE IN RADIANS */
/*                                   SPACECRAFT TO ECEF COORDINATES */
/*                                   TRANSFORMATION MATRIX */
/*                                   INSTRUMENT TO ECEF COORDINATES */
/*                                   TRANSFORMATION MATRIX */

/*     LOCAL VARIABLES */

/*                                   INSTRUMENT TO BODY COORDINATES */
/*                                   TRANSFORMATION MATRIX */
/*                                   INDICES */

/*     INCLUDE FILES */

/* *********************************************************************** */

/*     WE COMPUTE INSTRUMENT TO BODY COORDINATES TRANSFORMATION */
/*     MATRIX BY USING A SMALL ANGLE APPROXIMATION OF TRIGONOMETRIC */
/*     FUNCTIONS OF THE ROLL, PITCH AND YAW. */

    /* Parameter adjustments */
    at -= 4;
    a -= 4;

    /* Function Body */
    rpy[0] = 1. - (*p * *p + *y * *y) * .5;
    rpy[3] = -(*y);
    rpy[6] = *p;
    rpy[1] = *y + *p * *r__;
    rpy[4] = 1. - (*y * *y + *r__ * *r__) * .5;
    rpy[7] = -(*r__);
    rpy[2] = -(*p) + *r__ * *y;
    rpy[5] = *r__ + *p * *y;
    rpy[8] = 1. - (*p * *p + *r__ * *r__) * .5;

/*     MULTIPLICATION OF MATRICES A AND RPY */

    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    at[i__ + j * 3] = a[i__ + 3] * rpy[j * 3 - 3] + a[i__ + 6] * rpy[
		    j * 3 - 2] + a[i__ + 9] * rpy[j * 3 - 1];
/* L10: */
	}
/* L20: */
    }
    return 0;
} /* inst2er_ */

/* *********************************************************************** */
/* *********************************************************************** */
/* ** */
/* **   INTEGRAL SYSTEMS, INC. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   ROUTINE   : LPOINT */
/* **   SOURCE    : F.LPOINT */
/* **   LOAD NAME : ANY */
/* **   PROGRAMMER: IGOR LEVINE */
/* ** */
/* **   VER.    DATA    BY   COMMENT */
/* **   ----  --------  ---  --------------------------------------------- */
/* **   A     01/09/89  IL   INITIAL CREATION */
/* **   A     06/02/89  IL   COORDINATE AXES CHANGED ACCORDING TO */
/* **                        FORD'S DEFINITION IN SDAIP, DRL504-01 */
/* **   3     03/08/94  SC   IMPLEMENTED NEW FORMULAE FOR SCAN ANGLE */
/* **                        CORRECTIONS DUE TO MISALIGNMENTS */
/* *********************************************************************** */
/* ** */
/* **   THIS SUBROUTINE CONVERTS THE INSTRUMENT ELEVATION AND SCAN */
/* **   ANGLES TO THE RELATED GEOGRAPHIC LATITZDE AND LONGITZDE. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   CALLED BY       : ANY */
/* **   COMMONS MODIFIED: NONE */
/* **   INPUTS          : NONE */
/* **   OUTPUTS         : NONE */
/* **   ROUTINES CALLED : NONE */
/* ** */
/* *********************************************************************** */
/* *********************************************************************** */
/* Subroutine */ int lpoint_(doublereal *alpha, doublereal *zeta, doublereal *
	rlat, doublereal *rlon, integer *ierr)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal), tan(doublereal), sqrt(doublereal)
	    , atan(doublereal), atan2(doublereal, doublereal);

    /* Local variables */
    static doublereal d__, g[3], h__, u[3], d1, g1[3], q1, q2, ca, da, sa, cz,
	     dz;


/*     CALLING PARAMETERS */

/*                             ELEVATION ANGLE (RAD) */
/*                             SCAN ANGLE (RAD) */
/*                             LATITZDE IN RADIANS (OUTPUT) */
/*                             LONGITZDE IN RADIANS (OUTPUT) */
/*                             OUTPUT STATUS; 0 - POINT ON THE EARTH */
/*                             FOUND, 1 - INSTRUMENT POINTS OFF EARTH */

/*     LOCAL VARIABLES */

/*                          POINTING VECTOR IN EARTH CENTERED COORDINATES */
/*                          SLANT DISTANCE TO THE EARTH POINT (KM) */
/*                          WORK SPACE */
/*                          POINTING VECTOR IN INSTRUMENT COORDINATES */
/*                          COORDINATES OF THE EARTH POINT (KM) */
/*                                     WORK SPACE */

/*     INCLUDE FILES */

/* *********************************************************************** */
/* ** */
/* **   INTEGRAL SYSTEMS, INC. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   NAME      : ELCONS */
/* **   TYPE      : DATA AREA */
/* **   SOURCE    : I.ELCONS */
/* ** */
/* **   VER.    DATA    BY                COMMENT */
/* **   ----  --------  ----------------  -------------------------------- */
/* **   A     01/09/89  I. LEVINE         INITIAL CREATION */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   DESCRIPTION */
/* **   MATHEMATICAL AND EARTH-RELATED CONSTANTS */
/* ** */
/* *********************************************************************** */
/* *********************************************************************** */

/*                    DEGREES TO RADIANS CONVERSION PI/180 */
/*                    NOMINAL RADIAL DISTANCE OF SATELLITE (km) */
/*                    EARTH EQUATORIAL RADIUS (km) */
/*                    EARTH FLATTENING COEFFICIENT = 1-(BE/AE) */
/* *********************************************************************** */
/* *********************************************************************** */
/* *********************************************************************** */
/* **   INTEGRAL SYSTEMS, INC. */
/* *********************************************************************** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   NAME      : ELCOMM */
/* **   TYPE      : DATA AREA */
/* **   SOURCE    : I.ELCOMM */
/* ** */
/* **   VER.    DATA    BY                COMMENT */
/* **   ----  --------  ----------------  -------------------------------- */
/* **   A     01/09/89  I. LEVINE         INITIAL CREATION */
/* *********************************************************************** */
/* **   DESCRIPTION */
/* **   INSTRUMENT POSITION AND ATTITUDE VARIABLES AND TRANSFORMATION */
/* **   MATRIX */
/* *********************************************************************** */

/*     COMMON VARIABLES */

/*                      NORMALIZED S/C POSITION IN ECEF COORDINATES */
/*                      ECEF TO INSTRUMENT COORDINATES TRANSFORMATION */
/*                      USED IN SUBROUTINE LPOINT */
/*                          PITCH,ROLL,YAW ANGLES OF INSTRUMENT (RAD) */
/*                          PITCH,ROLL MISALIGNMENTS OF INSTRUMENT (RAD) */
    *ierr = 1;

/*     COMPUTES TRIGONOMETRIC FUNCTIONS OF THE SCAN AND ELEVATION */
/*     ANGLES CORRECTED FOR THE ROLL AND PITCH MISALIGNMENTS */

    ca = cos(*alpha);
    sa = sin(*alpha);
    cz = cos(*zeta);
    da = *alpha - elcomm_1.pma * sa * (1. / cz + tan(*zeta)) - elcomm_1.rma * 
	    (1. - ca / cz);
    dz = *zeta + elcomm_1.rma * sa;
/*                              CORRECTED SCAN ANGLE */
    cz = cos(dz);

/*     COMPUTES POINTING VECTOR IN INSTRUMENT COORDINATES */

    g[0] = sin(dz);
    g[1] = -cz * sin(da);
    g[2] = cz * cos(da);

/*     TRANSFORMS THE POINTING VECTOR TO EARTH FIXED COORDINATES */

    g1[0] = elcomm_1.bt[0] * g[0] + elcomm_1.bt[3] * g[1] + elcomm_1.bt[6] * 
	    g[2];
    g1[1] = elcomm_1.bt[1] * g[0] + elcomm_1.bt[4] * g[1] + elcomm_1.bt[7] * 
	    g[2];
    g1[2] = elcomm_1.bt[2] * g[0] + elcomm_1.bt[5] * g[1] + elcomm_1.bt[8] * 
	    g[2];

/*     COMPUTES COEFFICIENTS AND SOLVES A QUADRATIC EQUATION TO */
/*     FIND THE INTERSECT OF THE POINTING VECTOR WITH THE EARTH */
/*     SURFACE */

/* Computing 2nd power */
    d__1 = g1[0];
/* Computing 2nd power */
    d__2 = g1[1];
/* Computing 2nd power */
    d__3 = g1[2];
    q1 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 * 1.0067391845079681f;
    q2 = elcomm_1.xs[0] * g1[0] + elcomm_1.xs[1] * g1[1] + elcomm_1.xs[2] * 
	    1.0067391845079681f * g1[2];
    d__ = q2 * q2 - q1 * elcomm_1.q3;
    if (abs(d__) < 1e-9) {
	d__ = 0.;
    }

/*     IF THE DISCIMINANTE OF THE EQUATION, D, IS NEGATIVE, THE */
/*     INSTRUMENT POINTS OFF THE EARTH */

    if (d__ < 0.) {
	*rlat = 999999.;
	*rlon = 999999.;
	return 0;
    }
    d__ = sqrt(d__);

/*     SLANT DISTANCE FROM THE SATELLITE TO THE EARTH POINT */

    h__ = -(q2 + d__) / q1;

/*     CARTESIAN COORDINATES OF THE EARTH POINT */

    u[0] = elcomm_1.xs[0] + h__ * g1[0];
    u[1] = elcomm_1.xs[1] + h__ * g1[1];
    u[2] = elcomm_1.xs[2] + h__ * g1[2];

/*     SINUS OF GEOCENTRIC LATITZDE */

/* Computing 2nd power */
    d__1 = u[0];
/* Computing 2nd power */
    d__2 = u[1];
/* Computing 2nd power */
    d__3 = u[2];
    d1 = u[2] / sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);

/*     GEOGRAPHIC (GEODETIC) COORDINATES OF THE POINT */

    *rlat = atan(d1 * 1.0067391845079681f / sqrt(1. - d1 * d1));
    *rlon = atan2(u[1], u[0]);
    *ierr = 0;
    return 0;
} /* lpoint_ */

/* *********************************************************************** */
/* *********************************************************************** */
/* ** */
/* **   INTEGRAL SYSTEMS, INC. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   ROUTINE   : SNDELOC */
/* **   SOURCE    : F.SNDELOC */
/* **   LOAD NAME : ANY */
/* **   PROGRAMMER: IGOR LEVINE */
/* ** */
/* **   VER.    DATA    BY   COMMENT */
/* **   ----  --------  ---  --------------------------------------------- */
/* **   A     02/16/89  IL   INITIAL CREATION */
/* *********************************************************************** */
/* ** */
/* **   SNDELOC ACCEPTS THE MIRROR POSITION IN CYCLES AND INCREMENTS, */
/* **   SERVO ERROR VALUES, AND THE POSITIONAL OFFSETS FOR FOUR DETECTORS */
/* **   OF A SELECTED SOUNDER CHANNEL AND COMPUTES THE DETECTOR EARTH */
/* **   LOCATIONS IN LATITZDE/LONGITZDE COORDINATES. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   CALLED BY       : ANY */
/* **   COMMONS MODIFIED: NONE */
/* **   INPUTS          : NONE */
/* **   OUTPUTS         : NONE */
/* **   ROUTINES CALLED : LPOINT */
/* ** */
/* *********************************************************************** */
/* *********************************************************************** */
/* Subroutine */ int sndelo_(integer *cyew, integer *incew, integer *cyns, 
	integer *incns, doublereal *svew, doublereal *svns, doublereal *doff, 
	doublereal *geo)
{
    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal e, h__;
    static integer i__;
    static doublereal s, de, sc, ds, ev;
    static integer ier;
    static doublereal cose, sine;
    extern /* Subroutine */ int lpoint_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);


/*     CALLING PARAMETERS */

/*                         E-W CYCLES */
/*                         E-W INCREMENTS */
/*                         N-S CYCLES */
/*                         N-S INCREMENTS */
/*                         E-W SERVO ERROR IN RADIANS */
/*                         N-S SERVO ERROR IN RADIANS */
/*                         OFFSETS FOR 4 DETECTORS (RADIANS) */
/*                            DOFF(*,1) = E-W OFFSET */
/*                            DOFF(*,2) = N-S OFFSET */
/*                         GEOGRAPHIC COORDINATES RELATED TO 4 DETECTORS */
/*                            GEO(*,1) = LATITZDE IN RADIANS */
/*                            GEO(*,2) = LONGITZDE IN RADIANS */

/*     LOCAL VARIABLES */

/*     REAL*8 ALPHA,BETA */

/*     INCLUDE FILES */

/* *********************************************************************** */

/*     CONVERT THE MIRROR POSITION, GIVEN IN CYCLES AND INCREMENTS, TO */
/*     ELEVATION AND SCAN ANGLES */

/*     E=(CYNS*INCMAX(2)+INCNS)*ELVINC(2)-ELVMAX(2) */
/* *********************************************************************** */
/* *********************************************************************** */
/* ** */
/* **   INTEGRAL SYSTEMS, INC. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   NAME      : INSTCOMM */
/* **   TYPE      : DATA AREA */
/* **   SOURCE    : I.INSTCOMM */
/* ** */
/* **   VER.    DATA    BY                COMMENT */
/* **   ----  --------  ----------------  -------------------------------- */
/* **   A     02/16/89  I. LEVINE         INITIAL CREATION */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   DESCRIPTION */
/* **   COMMON AREA FOR INSTRUMENT-RELATED CONTROL PARAMETERS */
/* ** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     VARIABLES */
/*        CONSTANTS NEEDED TO PERFORM TRANSFORMATIONS BETWEEN THE */
/*        LATITUDE/LONGITUDE, LINE/PIXEL AND INSTRUMENT CYCLES/INCREMENTS */
/*        COORDINATES. */

/*                       NUMBER OF INCREMENTS PER CYCLE */
/*                       BOUNDS IN ELEVATION (RADIANS) */
/*                       BOUNDS IN SCAN ANGLE (RADIANS) */
/*                       CHANGE IN ELEVATION ANGLE PER INCREMENT (RAD) */
/*                       CHANGE IN SCAN ANGLE PER INCREMENT (RADIANS) */
/*                       ELEVATION ANGLE PER DETECTOR LINE (RADIANS) */
/*                       SCAN ANGLE PER PIXEL (RADIANS) */
    /* Parameter adjustments */
    geo -= 5;
    doff -= 5;

    /* Function Body */
    e = ((*cyns - 9) * instco_1.incmax[1] + *incns) * instco_1.elvinc[1] - 
	    instco_1.elvmax[1];
    s = (*cyew * instco_1.incmax[1] + *incew) * instco_1.scninc[1] - 
	    instco_1.scnmax[1];

/*     CORRECT ELEVATION AND SCAN ANGLES FOR SERVO ERRORS OBTAINING THE */
/*     TRUE MIRROR POINTING */

    e += *svns;
    s += *svew;
    sine = sin(e);
    cose = cos(e);
    h__ = instco_1.scnpx[1] * -2.;

/*     COMPUTE DETECTOR ROTATION OFFSETS FOR EACH DETECTOR */

/*      ALPHA = 0.643501D0 + E */
/*      BETA = 0.244979D0 - E */

/*      DOFF(1,1) = -0.064976D0 */
/*      DOFF(1,2) = 0.00042D0 */
/*      DOFF(2,1) = 0.00056D0 */
/*      DOFF(2,2) = 0.00014D0 */
/*      DOFF(3,1) = -0.064976D0 */
/*      DOFF(3,2) = -0.065396D0 */
/*      DOFF(4,1) = 0.00056D0 */
/*      DOFF(4,2) = -0.065116D0 */
/*      DOFF(1,1) = - 700.0D0*DCOS(ALPHA)*1.0D-6 */
/*      DOFF(1,2) =   700.0D0*DSIN(ALPHA)*1.0D-6 */
/*      DOFF(2,1) =   577.23479D0*DCOS(BETA)*1.D-6 */
/*      DOFF(2,2) =   577.23479D0*DSIN(BETA)*1.0D-6 */
/*      DOFF(3,1) = - 577.23479D0*DCOS(BETA)*1.0D-6 */
/*      DOFF(3,2) = - 577.23479D0*DSIN(BETA)*1.0D-6 */
/*      DOFF(4,1) =   700.0D0*DCOS(ALPHA)*1.0D-6 */
/*      DOFF(4,2) = - 700.0D0*DSIN(ALPHA)*1.0D-6 */

/*     COMPUTE EARTH LOCATIONS FOR FOUR DETECTORS */

    for (i__ = 1; i__ <= 4; ++i__) {
/*     COMPUTE POSITIONAL OFFSETS OF I-TH DETECTOR */
	de = (2.5f - i__) * instco_1.elvln[1] + doff[i__ + 8];
	ds = h__ + doff[i__ + 4];

/*     COMPUTE ELEVATION AND SCAN ANGLES RELATED TO I-TH DETECTOR */
/*     AND CORRECT THEM FOR THE DETECTOR POSITIONAL OFFSETS */

/*           EV=E+DOFF(I,2) */
/*           SC=S+DOFF(I,1) */

/*     CONVERT POSITIONAL OFFSETS TO ANGULAR OFFSETS AND */
/*     CORRECT ELEVATION AND SCAN ANGLES */
	ev = e + de * cose - ds * sine;
	sc = s + de * sine + ds * cose;
/*     TRANSFORM DETECTOR'S POINTING ANGLES TO GEOGRAPHIC COORDINATES */
/*     OF THE CORRESPONDING POINT ON THE EARTH SURFACE. */
/*     NOTE:  IF A DETECTOR LOOKS OFF THE EARTH, THE RELATED LATITZDE */
/*            AND LONGITZDE ARE SET TO 999999. */

	lpoint_(&ev, &sc, &geo[i__ + 4], &geo[i__ + 8], &ier);
	h__ = -h__;
/* L10: */
    }
    return 0;
} /* sndelo_ */

/* *********************************************************************** */
/* *********************************************************************** */
/* *********************************************************************** */
/* ** */
/* **   INTEGRAL SYSTEMS, INC. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   ROUTINE   : EVSC2LPF */
/* **   SOURCE    : F.EVSC2LPF */
/* **   LOAD NAME : ANY */
/* **   PROGRAMMER: IGOR LEVINE */
/* ** */
/* **   VER.    DATA    BY   COMMENT */
/* **   ----  --------  ---  --------------------------------------------- */
/* **   A     10/27/88  IL   INITIAL CREATION */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   THIS SUBROUTINE CONVERTS ELEVATION AND SCAN ANGLES */
/* **   TO THE FRACTIONAL LINE AND PIXEL NUMBERS. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   CALLED BY       : ANY */
/* **   COMMONS MODIFIED: NONE */
/* **   INPUTS          : NONE */
/* **   OUTPUTS         : NONE */
/* **   ROUTINES CALLED : NONE */
/* ** */
/* *********************************************************************** */
/* *********************************************************************** */
/* Subroutine */ int evsc2l_(integer *instr, doublereal *elev, doublereal *
	scan, doublereal *rl, doublereal *rp)
{

/*     CALLING PARAMETERS */

/*                       INSTRUMENT CODE (1-IMAGER, 2-SOUNDER) */
/*                       ELEVATION ANGLE IN RADIANS */
/*                       SCAN ANGLE IN RADIANS */
/*                       LINE NUMBER */
/*                       PIXEL NUMBER */

/*     LOCAL VARIABLES - NONE */


/*     INCLUDE FILES */

/* ************************************************************** */

/*     COMPUTE FRACTIONAL LINE NUMBER */

/* *********************************************************************** */
/* *********************************************************************** */
/* ** */
/* **   INTEGRAL SYSTEMS, INC. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   NAME      : INSTCOMM */
/* **   TYPE      : DATA AREA */
/* **   SOURCE    : I.INSTCOMM */
/* ** */
/* **   VER.    DATA    BY                COMMENT */
/* **   ----  --------  ----------------  -------------------------------- */
/* **   A     02/16/89  I. LEVINE         INITIAL CREATION */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   DESCRIPTION */
/* **   COMMON AREA FOR INSTRUMENT-RELATED CONTROL PARAMETERS */
/* ** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     VARIABLES */
/*        CONSTANTS NEEDED TO PERFORM TRANSFORMATIONS BETWEEN THE */
/*        LATITUDE/LONGITUDE, LINE/PIXEL AND INSTRUMENT CYCLES/INCREMENTS */
/*        COORDINATES. */

/*                       NUMBER OF INCREMENTS PER CYCLE */
/*                       BOUNDS IN ELEVATION (RADIANS) */
/*                       BOUNDS IN SCAN ANGLE (RADIANS) */
/*                       CHANGE IN ELEVATION ANGLE PER INCREMENT (RAD) */
/*                       CHANGE IN SCAN ANGLE PER INCREMENT (RADIANS) */
/*                       ELEVATION ANGLE PER DETECTOR LINE (RADIANS) */
/*                       SCAN ANGLE PER PIXEL (RADIANS) */
    *rl = (instco_1.elvmax[*instr - 1] - *elev) / instco_1.elvln[*instr - 1];
    if (*instr == 1) {
	*rl += 4.5;
    } else {
	*rl += 2.5;
    }

/*     COMPUTE FRACTIONAL PIXEL NUMBER */

    *rp = (instco_1.scnmax[*instr - 1] + *scan) / instco_1.scnpx[*instr - 1] 
	    + 1.;
    return 0;
} /* evsc2l_ */

/* *********************************************************************** */
/* ** */
/* **   INTEGRAL SYSTEMS, INC. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   ROUTINE   : EVLN */
/* **   SOURCE    : F.EVLN */
/* **   LOAD NAME : ANY */
/* **   PROGRAMMER: IGOR LEVINE */
/* ** */
/* **   VER.    DATA    BY   COMMENT */
/* **   ----  --------  ---  --------------------------------------------- */
/* **   A     10/27/88  IL   INITIAL CREATION */
/* *********************************************************************** */
/* ** */
/* **   THIS FUNCTION CONVERTS FRACTIONAL LINE NUMBER TO ELEVATION ANGLE */
/* **   IN RADIANS. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   CALLED BY       : ANY */
/* **   COMMONS MODIFIED: NONE */
/* **   INPUTS          : NONE */
/* **   OUTPUTS         : NONE */
/* **   ROUTINES CALLED : NONE */
/* ** */
/* *********************************************************************** */
/* *********************************************************************** */
doublereal evln_(integer *instr, doublereal *rline)
{
    /* System generated locals */
    doublereal ret_val;


/*     CALLING PARAMETERS */

/*                       INSTRUMENT CODE (1-IMAGER, 2-SOUNDER) */
/*                       FRACTIONAL LINE  NUMBER */

/*     LOCAL VARIABLES - NONE */


/*     INCLUDE FILES */

/* *********************************************************************** */
/* *********************************************************************** */
/* *********************************************************************** */
/* ** */
/* **   INTEGRAL SYSTEMS, INC. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   NAME      : INSTCOMM */
/* **   TYPE      : DATA AREA */
/* **   SOURCE    : I.INSTCOMM */
/* ** */
/* **   VER.    DATA    BY                COMMENT */
/* **   ----  --------  ----------------  -------------------------------- */
/* **   A     02/16/89  I. LEVINE         INITIAL CREATION */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   DESCRIPTION */
/* **   COMMON AREA FOR INSTRUMENT-RELATED CONTROL PARAMETERS */
/* ** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     VARIABLES */
/*        CONSTANTS NEEDED TO PERFORM TRANSFORMATIONS BETWEEN THE */
/*        LATITUDE/LONGITUDE, LINE/PIXEL AND INSTRUMENT CYCLES/INCREMENTS */
/*        COORDINATES. */

/*                       NUMBER OF INCREMENTS PER CYCLE */
/*                       BOUNDS IN ELEVATION (RADIANS) */
/*                       BOUNDS IN SCAN ANGLE (RADIANS) */
/*                       CHANGE IN ELEVATION ANGLE PER INCREMENT (RAD) */
/*                       CHANGE IN SCAN ANGLE PER INCREMENT (RADIANS) */
/*                       ELEVATION ANGLE PER DETECTOR LINE (RADIANS) */
/*                       SCAN ANGLE PER PIXEL (RADIANS) */
    if (*instr == 1) {
	ret_val = instco_1.elvmax[*instr - 1] * 1. - (*rline - 4.5f) * (
		instco_1.elvln[*instr - 1] * 1.);
    } else {
	ret_val = instco_1.elvmax[*instr - 1] * 1. - (*rline - 2.5) * (
		instco_1.elvln[*instr - 1] * 1.);
    }
    return ret_val;
} /* evln_ */

/* *********************************************************************** */
/* ** */
/* **   INTEGRAL SYSTEMS, INC. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   ROUTINE   : SCPX */
/* **   SOURCE    : F.SCPX */
/* **   LOAD NAME : ANY */
/* **   PROGRAMMER: IGOR LEVINE */
/* ** */
/* **   VER.    DATA    BY   COMMENT */
/* **   ----  --------  ---  --------------------------------------------- */
/* **   A     09/22/87  IL   INITIAL CREATION */
/* *********************************************************************** */
/* ** */
/* **   THIS FUNCTION CONVERTS FRACTIONAL PIXEL NUMBER TO SCAN ANGLE */
/* **   IN RADIANS. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   CALLED BY       : ANY */
/* **   COMMONS MODIFIED: NONE */
/* **   INPUTS          : NONE */
/* **   OUTPUTS         : NONE */
/* **   ROUTINES CALLED : NONE */
/* ** */
/* *********************************************************************** */
/* *********************************************************************** */
doublereal scpx_(integer *instr, doublereal *pix)
{
    /* System generated locals */
    doublereal ret_val;


/*     CALLING PARAMETERS */

/*                       INSTRUMENT CODE (1-IMAGER, 2-SOUNDER) */
/*                       FRACTIONAL PIXEL NUMBER */

/*     LOCAL VARIABLES */


/*     INCLUDE FILES */

/* *********************************************************************** */
/* *********************************************************************** */
/* *********************************************************************** */
/* ** */
/* **   INTEGRAL SYSTEMS, INC. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   NAME      : INSTCOMM */
/* **   TYPE      : DATA AREA */
/* **   SOURCE    : I.INSTCOMM */
/* ** */
/* **   VER.    DATA    BY                COMMENT */
/* **   ----  --------  ----------------  -------------------------------- */
/* **   A     02/16/89  I. LEVINE         INITIAL CREATION */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   DESCRIPTION */
/* **   COMMON AREA FOR INSTRUMENT-RELATED CONTROL PARAMETERS */
/* ** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     VARIABLES */
/*        CONSTANTS NEEDED TO PERFORM TRANSFORMATIONS BETWEEN THE */
/*        LATITUDE/LONGITUDE, LINE/PIXEL AND INSTRUMENT CYCLES/INCREMENTS */
/*        COORDINATES. */

/*                       NUMBER OF INCREMENTS PER CYCLE */
/*                       BOUNDS IN ELEVATION (RADIANS) */
/*                       BOUNDS IN SCAN ANGLE (RADIANS) */
/*                       CHANGE IN ELEVATION ANGLE PER INCREMENT (RAD) */
/*                       CHANGE IN SCAN ANGLE PER INCREMENT (RADIANS) */
/*                       ELEVATION ANGLE PER DETECTOR LINE (RADIANS) */
/*                       SCAN ANGLE PER PIXEL (RADIANS) */
    ret_val = (*pix * 1. - 1.) * (instco_1.scnpx[*instr - 1] * 1.) - 
	    instco_1.scnmax[*instr - 1] * 1.;
    return ret_val;
} /* scpx_ */

/* *********************************************************************** */
/* ** */
/* **   INTEGRAL SYSTEMS, INC. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT */
/* **   SYSTEM    : EARTH LOCATION USERS GUIDE */
/* **   ROUTINE   : GATT */
/* **   SOURCE    : F.GATT */
/* **   LOAD NAME : ANY */
/* **   PROGRAMMER: IGOR LEVINE */
/* ** */
/* **   VER.    DATA    BY   COMMENT */
/* **   ----  --------  ---  --------------------------------------------- */
/* **   A     12/01/88  IL   INITIAL CREATION */
/* ** */
/* *********************************************************************** */
/* ** */
/* **    THIS FUNCTION COMPUTES AN ATTITZDE/MISALIGNMENT ANGLE FROM A */
/* **    GIVEN SUBSET OF THE O&A PARAMETERS IN GVAR BLOK 0. */
/* **    ARGUMENT K0 INDICATES THE FIRST WORD OF THE SUBSET. */
/* ** */
/* *********************************************************************** */
/* ** */
/* **   CALLED BY       : LMODEL */
/* **   COMMONS MODIFIED: NONE */
/* **   INPUTS          : NONE */
/* **   OUTPUTS         : NONE */
/* **   ROUTINES CALLED : NONE */
/* ** */
/* *********************************************************************** */
/* *********************************************************************** */
doublereal gatt_(integer *k0, doublereal *rec, doublereal *wa, doublereal *te)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;
    static doublereal equiv_0[1], equiv_1[1];

    /* Builtin functions */
    double exp(doublereal), cos(doublereal), pow_di(doublereal *, integer *);

    /* Local variables */
    static integer i__;
#define j ((integer *)equiv_0)
    static integer k, l;
#define m ((integer *)equiv_1)
    static integer ll;
    static doublereal ir;
#define jr (equiv_0)
#define mr (equiv_1)
    static doublereal att;


/*     CALLING PARAMETERS */

/*                    STARTING POSITION OF A PARAMETER SUBSET IN THE */
/*                    O&A SET */
/*                    INPUT O&A PARAMETER SET */
/*                    INPUT SOLAR ORBIT ANGLE IN RADIANS */
/*                    INPUT EXPONENTIAL TIME DELAY FROM EPOCH (MINUTES) */

/*     LOCAL VARIABLES */


/*     INCLUDE FILES */

/* *********************************************************************** */

/*     CONSTANT COMPONENT */

    /* Parameter adjustments */
    --rec;

    /* Function Body */
    k = *k0;
    att = rec[k + 2];

/*     COMPUTES THE EXPONENTIAL TERM */

    if (*te >= 0. && rec[k + 1] > 0.) {
	att += rec[k] * exp(-(*te) / rec[k + 1]);
    }

/*     EXTRACTS THE NUMBER OF SINUSOIDS */

    ir = rec[k + 3];
    i__ = (integer) ir;

/*     CALCULATION OF SINUSOIDS */

    i__1 = i__;
    for (l = 1; l <= i__1; ++l) {
	att += rec[k + (l << 1) + 2] * cos(*wa * l + rec[k + (l << 1) + 3]);
/* L10: */
    }

/*     POINTER TO THE NUMBER OF MONOMIAL SINUSOIDS */

    k += 34;

/*     EXTACTS NUMBER OF MONOMIAL SINUSOIDS */

    ir = rec[k];
    i__ = (integer) ir;
/*     KKK=REC(K) */

/*     COMPUTES MONOMIAL SINUSOIDS */

    i__1 = i__;
    for (l = 1; l <= i__1; ++l) {
	ll = k + l * 5;

/*          ORDER OF SINUSOID */

	*jr = rec[ll - 4];

/*          ORDER OF MONOMIAL SINUSOID */

	*mr = rec[ll - 3];

	d__1 = *wa - rec[ll];
	att += rec[ll - 2] * pow_di(&d__1, m) * cos(*j * *wa + rec[ll - 1]);
/* L20: */
    }
    ret_val = att;
    return ret_val;
} /* gatt_ */

#undef mr
#undef jr
#undef m
#undef j


