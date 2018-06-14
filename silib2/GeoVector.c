/***************************************************************************
 *
 *       Module:    GeoVector.c
 *      Program:    SILib2
 *   Created by:    Darren Muff on 21/07/2013.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      C functions for handling geo-points. The functions use the SPVector class
 *      which is used for storing ECEF coords of the point. The functions use the
 *      WGS84(G873) ellipsoid and uses EGM96 for its GEOID. (The GEOID is the
 *      eqi-potential model for the Earth which is equivalent
 *      to Mean Sea Level.)
 *
 *      The calculations here make use of the references at the bottom which should
 *      be used to understand the maths.
 * 
 *      Additionally the file has two routines for converting between Universal 
 *      Transverse Mercator (UTM) and Latitide / Longitude.
 *      The calculations here are taken from [5] and adapted from code kindly provided by
 *      David Yip and Bob Harrison, Lockheed Martin, Denver
 *
 *
 *	References:
 *		[0] Burtch, R.; "A Comparision of Methods Used in Rectangular to Geodetic Coordinate
 *		Transformations"; ACSM Annual Conference and Lechnology Exhibition; Orlando, FL; Apr. 21-26, 2006; 21 pages.
 *
 *		[1] Lin, K-C. and J. Wang, 1995. “Transformation from Geocentric to Geodetic
 *		Coordinates using Newton’s Iteration”, Bulletin Géodésique, 69(4): 300-303.
 *
 *		[2] Bowring, B.R., 1976. “Transformation from Spatial to geographical
 *		Coordinates”, Survey Review, 23(181): 323-327.
 *
 *		[3] Gerdan, G.P., R.E. Deakin, 1999. “Transforming Cartesian Coordinates
 *		X, Y, Z to Geographical Coordinates”, The Australian Surveyor, 44(1): 55-63.
 *
 *		[4] http://en.wikipedia.org/wiki/Geodetic_system#From_geodetic_to_ECEF
 *
 *		[5] Hager, J. W., Behensky, J. F., Drew, B. W., "The Universal Grids:
 *			Universal Transverse Mercator (UTM) and Universal Polar Stereographic
 *			(UPS)", DMA Technical Manual, DMATM 8358.2, 18 September 1989.
 * 
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  13th September 2011
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 * THE SOFTWARE IN ITS ENTIRETY OR ANY SUBSTANTIAL PORTION SHOULD NOT BE
 * USED AS A WHOLE OR COMPONENT PART OF DERIVATIVE SOFTWARE THAT IS BEING
 * SOLD TO THE GOVERNMENT OF THE UNITED KINGDOM OF GREAT BRITAIN AND NORTHERN
 * IRELAND.
 *
 ***************************************************************************/
#include "SIlib2.h"
static const double correctionCoef1[65341];		///< 360*361/2 + 361
static const double correctionCoef2[65341];		///< 360*361/2 + 361
static const double potentialCoef1[65341];		///< 360*361/2 + 361
static const double potentialCoef2[65341];		///< 360*361/2 + 361

// Private function to calculate Legendre polynomial coefficients
//
void LEGFDN(int M,double THETA,double *RLEG, double *DLEG,int NMX, int *IR,double *RLNN,int IFLAG,
            double *DRTS, double *DIRT);
// Private function to sum Legendre coefficients
//
void DSCML(double RLON, int NMAX, double *SINML, double *COSML);
// private function to calculate the geoid undulation from the
// EGM96 potential coefficient model
//
void HUNDU(double *UNDU,int NMAX,double *P,const double *HC,const double *HS,
           double *SINML,double *COSML,double GR,double RE,const double *CCoef,
           const double *CS,double *HACO);

// Function to convert any angle to a correct latitude
// between -90 (south pole) and +90 (north pole)
//
double fixlat(double lat){
	lat = fmod(lat,360.);
	if(lat > 90. && lat <= 270.) lat = 180.-lat;
	if(lat > 270. && lat < 360.) lat = lat - 360.;
	if(lat < -90. && lat >= -270.) lat = -180.-lat;
	if(lat < -270. && lat > -360.) lat = lat + 360.;
	return(lat);
}

// function to convert any angle to a correct longitude
// between  -180 deg and +180 deg
//
double fixlon(double lon){
	lon = fmod(lon,360.);
	if(lon > 180.)lon=lon-360.;
	if(lon <= -180)lon=lon+360;
	return(lon);
	
}

// Function that returns the height difference between the geoid and ellipsoid
// for a given lat, lon
//
double GeoidHeightToEllipsoidHeight(double GeoidHeight, double lattitude_deg, double longitude_deg) {
    return ( GeoidHeight - undulation(lattitude_deg, longitude_deg ) ) ;
}

// Function that returns the height difference between the ellipsoid and geoid
// for a given lat, lon
//
double EllipsoidHeightToGeoidHeight(double ellipsoidHeight, double lattitude_deg, double longitude_deg) {
    return ( ellipsoidHeight + undulation(lattitude_deg, longitude_deg ) ) ;
}

// Function to convert a lat,lon.alt to an ECEF coordinate
// The alt can be specified as being ref to 'Ellipsoid' or 'Geoid' in
// the 4th parameter
//
SPStatus latLonAltToEcef(double latitude_deg, double longitude_deg, double alt, enum GeoAltReference GeoidOrEllipsoid, SPVector *ecef, SPStatus *status){
	
	// This function taken from [4]
    //
	
    CHECK_STATUS(status);
    double latitude_rad, longitude_rad;
    latitude_rad  = latitude_deg  * GeoConsts_DEGTORAD;
    longitude_rad = longitude_deg * GeoConsts_DEGTORAD;
	double e2 = (2/GeoConsts_Pff) - (1/(GeoConsts_Pff*GeoConsts_Pff));
	double chi = sqrt( 1-e2*sin(latitude_rad)*sin(latitude_rad));
	double altEllipsoid ;
    if(GeoidOrEllipsoid == Ellipsoid){
		altEllipsoid=alt;
	}else if(GeoidOrEllipsoid == Geoid){
        altEllipsoid = GeoidHeightToEllipsoidHeight(alt, latitude_deg, longitude_deg) ;
	}else{
		printf("ERROR: invalid reference for altitude ");
        status->status = UNKNOWN_GEOALT_REFERENCE ;
		return *status ;
	}
    
	ecef->x = (GeoConsts_Re/chi + altEllipsoid)*cos(latitude_rad)*cos(longitude_rad);
	ecef->y = (GeoConsts_Re/chi + altEllipsoid)*cos(latitude_rad)*sin(longitude_rad);
	ecef->z = (GeoConsts_Re*(1-e2)/chi + altEllipsoid)*sin(latitude_rad);
	return *status ;
}

// Function to convert ECEF parameters (X,Y,Z) to a lat,lon,alt SPVector
// The altitude in the return vector is referenced to the ellipsoid if the 4th parameter
// is set t be 'Ellipsoid' or referenced to the geoid if the 4th parameter is set
// to 'Geoid'
//
SPStatus ecefToLatLonAlt(double X, double Y, double Z, enum GeoAltReference GeoidOrEllipsoid, SPVector *latlonalt, SPStatus *status){
	
	// This method is the Lin and Wang Iterative method taken from [0]. It is detailed
    // in the attached documentation [0] and in [1]. It is selected as it has the accuracy
    // and convergence speed of Bowrings iterative method [2] but does not use
    // trigonometric or inverse functions which means its the most efficient in terms
    // of execution time (according to Gerden and Deakin [3]).
    //

    CHECK_STATUS(status);
    int cnt=0;

	// Key parameter definitions are:
	//
	// X,Y,Z are Cartesian geocentric coords of point in metres
    //
	
	double a=GeoConsts_Re;					// Semi-major axis
	double b=a*(1-(1/GeoConsts_Pff));       // Semi-minor axis
	double p = sqrt(X*X+Y*Y);				// Radius to equatorial axis
	double m;                               // m is a parameter that describes where along the normal, the point (x, y, z)
                                            // is located. For example,when m=0,x=X,y=Y,andz=Z.
	double fm;                              // function value for Newton-Raphson iteration
	double fmdash;                          // derivative of fm with respect to m
	double Ppdd;                            // Radius to equatorial axis of a point on the geoid
	double Zpdd;                            // Cartesian Z coord of the point on the geoid
                                            // whose normal to the geoid runs through
                                            // the ecef point (X,Y,Z)
	
	// To understand the following calculations please read [0]
    //
	
	// initial estimate for m
    //
	m = (a*b*pow(a*a*Z*Z+b*b*p*p,1.5)-a*a*b*b*(a*a*Z*Z+b*b*p*p))/(2*(a*a*a*a*Z*Z+b*b*b*b*p*p));
	fm = (p*p/((a+2*m/a)*(a+2*m/a)))+(Z*Z/((b+2*m/b)*(b+2*m/b)))-1;
	
	while(fabs(fm)>1.0e-6){			// one micron should be enough ;-)
		fmdash = -4*((p*p/(a*(a+2*m/a)*(a+2*m/a)*(a+2*m/a))) + (Z*Z/(b*(b+2*m/b)*(b+2*m/b)*(b+2*m/b))));
		m = m*fm/fmdash;
		fm = (p*p/((a+2*m/a)*(a+2*m/a)))+(Z*Z/((b+2*m/b)*(b+2*m/b)))-1;
        if(cnt++ >= 10){
            if(status->debug>=30){
                fprintf(stderr, "Warning poor conversion accurate for ecefToLatLonALt: File: %s, line: %d\n",__FILE__,__LINE__);
            }
            break ;
        }
	}
	
	double phi, lambda, h;
	
	Ppdd = fabs(p/(1+2*m/(a*a)));
	Zpdd = b*sqrt(1-((Ppdd*Ppdd)/(a*a)));
	phi = atan(a*a*Zpdd/(b*b*Ppdd));				// latitude
	lambda = atan(Y/X);								// longitude
	h = sqrt((p-Ppdd)*(p-Ppdd)+(Z-Zpdd)*(Z-Zpdd));	// height
	
	latlonalt->lat = phi    * GeoConsts_RADTODEG ;
	latlonalt->lon = lambda * GeoConsts_RADTODEG ;
    if(GeoidOrEllipsoid == Ellipsoid){
        latlonalt->z = h;
	}else if(GeoidOrEllipsoid == Geoid){
        latlonalt->z = EllipsoidHeightToGeoidHeight(h, latlonalt->lat, latlonalt->lon) ;
	}else{
		printf("ERROR: invalid reference for altitude ");
        status->status = UNKNOWN_GEOALT_REFERENCE ;
		return *status ;
	}

	return *status ;
}

// Undulation is the term refering to the difference in height in metres between
// the Earth's gravitional potential model and the Earth's ellipsoid model for a
// given lat and long.
// This function calculates the geoid undulation from the EGM96 potential coefficient model
// including the height anomaly to geoid undulation correction term and a
// correction term to have the undulations refer to the WGS84 ellipsoid. The
// geoid undulation unit is metres.
//
double undulation(double lattitude_deg, double longitude_deg){
    
	// Calculate the geocentric distance to the point, the geocentric
    // latitude, and an approximate value of the normal gravity at
    // the point using the constants of WGS84(G873)
    //
	double rlat = GeoConsts_RADTODEG*lattitude_deg;
	double rlon = GeoConsts_RADTODEG*longitude_deg;
	
	double T1 = sin(rlat)*sin(rlat);
	double e2 = (2/GeoConsts_Pff) - (1/(GeoConsts_Pff*GeoConsts_Pff));
	double N = GeoConsts_Re/sqrt(1.0 - e2*T1);
	double HT= 0.0; // initial value of height
	double T2 = (N+HT)*cos(rlat);
	double X = T2*cos(rlon);
	double Y = T2*sin(rlon);
	double Z = (N*(1.0-e2)+HT)*sin(rlat);
	double K = 0.00193185265246;
	// Compute the geocentric radius
    //
	double RE = sqrt(X*X + Y*Y + Z*Z);
	
    // Compute the geocentric latitude
    //
	double RLAT = atan(Z/sqrt(X*X + Y*Y));
	
    // Compute the normal gravity. Units are M/sec^2
    //
	double GR = GeoConsts_Geqt*(1.0+K*T1)/sqrt(1.0-e2*T1);
    
	// SETTING IFLAG=1 PREVENTS LEGENDRE FUNCTION DERIVATIVES BEING TAKEN
	// IN SUBROUTINE LEGFDN
    //
	int IFLAG=1;
	int IR=0;
	int k=361;
	
	RLAT = 1.5707963267948966-RLAT;
	double RLEG[361], DLEG[361], RLNN[361];
	double DRTS[1300], DIRT[1300];
	
	int LOC;
	double P[65341];
	int M;
	for(int J=1; J<=k; J++){
		M=J-1;
		LEGFDN(M,RLAT,RLEG,DLEG,360,&IR,RLNN,IFLAG,DRTS,DIRT);
		for(int I=J; I<=k; I++){
			N=I-1;
			LOC=(N*(N+1))/2+M+1;
			P[LOC-1] = RLEG[I-1];
		}
	}
	
	double SINML[361], COSML[361];
	double NMAX = 360;
	DSCML(rlon,NMAX,SINML,COSML);
	double U;
	double HACO;
	
	HUNDU(&U,NMAX,P, potentialCoef1,potentialCoef2,SINML,COSML,GR,RE,correctionCoef1,correctionCoef2,&HACO);
	
	return(U);
}

// convert a latitude and longitude to a Universal Transverse Mercator (UTM)
// easting and northing (and the associated UTMZone)
//
SPStatus LLtoUTM(double latitude_deg, double longitude_deg, double *UTMEasting, double *UTMNorthing, int *UTMZone, SPStatus *status){
	double latIn = latitude_deg;
	double lonIn = longitude_deg;
	double S1, utmz, lon_cm_deg, lon_cm, lat, lon, dlon, FN, v;
	double sin_lat, cos_lat, cos_lat2, cos_lat3,cos_lat4, cos_lat5, cos_lat6;
	double cos_lat7, cos_lat8, tan_lat, tan_lat2, tan_lat4, tan_lat6, T1, v_k0;
	double v_sinlat_k0,T2,temp3,T3,temp4,T4,temp5,T5,T6,T7,temp8,T8,temp9,T9;
	double dlon2,dlon3,dlon4,dlon5,dlon6,dlon7,dlon8;
	
    CHECK_STATUS(status) ;
    
	// Calc UTMzone - based on first value of input Longitude
    //
	utmz = floor((lonIn + 180.f)/6.f + 0.5f);
	if(latIn < 0) utmz = -utmz;
	*UTMZone = (int)utmz;
	
    
	lon_cm_deg = (abs(*UTMZone)*6.) - 180. - 3.;
	lon_cm = lon_cm_deg * GeoConsts_DEGTORAD;
    
	// Convert input degrees to radians
    //
    lat = latIn * GeoConsts_DEGTORAD;
    lon = lonIn * GeoConsts_DEGTORAD;
    
    // Calc Central Meridian for UTMzone
    //
    dlon = lon - lon_cm;
	FN = (lat >= 0) ? 0 : 10000000;
    
	// Calculate ellipsoid parameters that vary with lat
    //
	// p = a*(1-e2) / pow( (1 - e2 * pow(sin(lat), 2)), 1.5 );
	v = GeoConsts_a / pow( (1 - GeoConsts_e2* pow(sin(lat), 2)), 0.5);
	// Calc meridional arc (S)
    //
	S1 = GeoConsts_Ap*lat - GeoConsts_Bp*sin(2*lat) + GeoConsts_Cp*sin(4*lat) - GeoConsts_Dp*sin(6*lat) + GeoConsts_Ep*sin(8*lat);
	// Calculate General Equations
    //
	sin_lat  = sin(lat);
	cos_lat  = cos(lat);
	cos_lat2 = cos_lat * cos_lat;
	cos_lat3 = cos_lat2 * cos_lat;
	cos_lat4 = cos_lat3 * cos_lat;
	cos_lat5 = cos_lat4 * cos_lat;
	cos_lat6 = cos_lat5 * cos_lat;
	cos_lat7 = cos_lat6 * cos_lat;
	cos_lat8 = cos_lat7 * cos_lat;
	tan_lat  = tan(lat);
	tan_lat2 = tan_lat * tan_lat;
	tan_lat4 = tan_lat2 * tan_lat2;
	tan_lat6 = tan_lat4 * tan_lat2;
	T1 = S1 * GeoConsts_k0;
	v_k0 = v * GeoConsts_k0;
	v_sinlat_k0 = v_k0 * sin_lat;
	T2 = v_sinlat_k0/2.0 * cos_lat;
	temp3 = (5.0 - tan_lat2 + 9.0*GeoConsts_ep2*cos_lat2 + 4.0*GeoConsts_ep4*cos_lat4);
	T3 = v_sinlat_k0 * cos_lat3 /24.0 * temp3;
	temp4 = 61.0 - 58.0*tan_lat2 + tan_lat4 + 270.0*GeoConsts_ep2*cos_lat2
    - 330.0*tan_lat2*GeoConsts_ep2*cos_lat2 + 445.0*GeoConsts_ep4*cos_lat4 + 324.0*GeoConsts_ep6*cos_lat6
    - 680.0*tan_lat2*GeoConsts_ep4*cos_lat4 + 88.0*GeoConsts_ep8*cos_lat8
    - 600.0*tan_lat2*GeoConsts_ep6*cos_lat6 - 192.0*tan_lat2*GeoConsts_ep8*cos_lat8;
	T4 = v_sinlat_k0 * cos_lat5 / 720.0 * temp4;
	temp5 = 1385.0 - 3111.0*tan_lat2 + 543.0*tan_lat4 - tan_lat6;
	T5 = v_sinlat_k0 * cos_lat7 / 40320.0 * temp5;
	T6 = v_k0 * cos_lat;
	T7 = v_k0/6 * cos_lat3 * (1.0 - tan_lat2 + GeoConsts_ep2*cos_lat2);
	temp8 = 5.0 - 18.0*tan_lat2 + tan_lat4 + 14.0*GeoConsts_ep2*cos_lat2
    - 58.0*tan_lat2*GeoConsts_ep2*cos_lat2 + 1.0*GeoConsts_ep4*cos_lat4 + 4.0*GeoConsts_ep6*cos_lat6
    - 64.0*tan_lat2*GeoConsts_ep4*cos_lat4 - 24.0*tan_lat2*GeoConsts_ep6*cos_lat6;
	T8 = v_k0/120.0 * cos_lat5 * temp8;
	temp9 = 61.0 - 479.0*tan_lat2 + 179.0*tan_lat4 - tan_lat6;
	T9 = v_k0/5040.0 * cos_lat7 * temp9;
    
	dlon2 = dlon  * dlon;
	dlon3 = dlon2 * dlon;
	dlon4 = dlon3 * dlon;
	dlon5 = dlon4 * dlon;
	dlon6 = dlon5 * dlon;
	dlon7 = dlon6 * dlon;
	dlon8 = dlon7 * dlon;
    
	// NORTHING (in meters)
    //
	*UTMNorthing = (double)FN + (T1 + dlon2*T2 + dlon4*T3 + dlon6*T4 + dlon8*T5);
	// EASTING (in meters)
    //
	*UTMEasting  = (double)GeoConsts_FE + (dlon*T6 + dlon3*T7 + dlon5*T8 + dlon7*T9);
    
	return *status ;
}

// convert a Universal Transverse Mercator (UTM) easting and northing (and the associated UTMZone)
// to a latitude and longitude
//
SPStatus UTMtoLL(double UTMEasting, double UTMNorthing, int UTMZone, double *latitude_deg, double *longitude_deg, SPStatus *status){
	
	double lon_cm_deg,lon_cm,FN,dE,tmd,ftlat;
	double S1,p1=0,v,cos_f1,cos_f2,cos_f3,cos_f4,cos_f5,cos_f6,cos_f7,cos_f8;
	double tan_f1,tan_f2,tan_f4,tan_f6;
	double v2,v3,v4,v5,v6,v7;
	double T10,temp11,T11,temp12,T12,temp13,T13,T14,T15,T16,temp16,T17;
	double dE2,dE3,dE4,dE5,dE6,dE7,dE8,lat,lon;
	int j;
    
    CHECK_STATUS(status) ;
	
	// Calc Central Meridian for UTMzone
    //
	lon_cm_deg = (abs(UTMZone)*6.) - 180. - 3.;
	lon_cm = lon_cm_deg * GeoConsts_DEGTORAD;
	
	// Set False Northing & False Easting
    //
	FN = (UTMZone >= 0) ? 0 : 10000000;
    
	dE = UTMEasting - GeoConsts_FE;
	
    // Calc footpoint latitude
    //
	tmd = GeoConsts_S + ((UTMNorthing - FN)/GeoConsts_k0);
	ftlat = tmd / GeoConsts_p;
	for (j=1; j<=5; j++){
		S1 = GeoConsts_Ap*ftlat - GeoConsts_Bp*sin(2.0*ftlat) + GeoConsts_Cp*sin(4.0*ftlat) - GeoConsts_Dp*sin(6.0*ftlat) + GeoConsts_Ep*sin(8.0*ftlat);
		p1 = GeoConsts_a*(1-GeoConsts_e2) / pow( (1 - GeoConsts_e2 * pow(sin(ftlat), 2)), 1.5);
		ftlat = ftlat + (tmd-S1)/p1;
	}
	v = GeoConsts_a / pow( (1 - GeoConsts_e2 * pow(sin(ftlat), 2)), 0.5);
    
	// Calculate General Equations
    //
	cos_f1 = cos(ftlat);
	cos_f2 = cos_f1 * cos_f1;
	cos_f3 = cos_f2 * cos_f1;
	cos_f4 = cos_f3 * cos_f1;
	cos_f5 = cos_f4 * cos_f1;
	cos_f6 = cos_f5 * cos_f1;
	cos_f7 = cos_f6 * cos_f1;
	cos_f8 = cos_f7 * cos_f1;
	tan_f1 = tan(ftlat);
	tan_f2 = tan_f1 * tan_f1;
	tan_f4 = tan_f2 * tan_f2;
	tan_f6 = tan_f4 * tan_f2;
	v2 = v  * v;
	v3 = v2 * v;
	v4 = v3 * v;
	v5 = v4 * v;
	v6 = v5 * v;
	v7 = v6 * v;
	T10 = tan_f1 / (p1 * v * GeoConsts_k02_2);
	temp11 = 5.0 + 3.0*tan_f2 + GeoConsts_ep2*cos_f2 - GeoConsts_ep4_4*cos_f4 - tan_f2*GeoConsts_ep2_9*cos_f2;
	T11 = tan_f1 / (p1 * v3 * GeoConsts_k04_24) * temp11;
	temp12 = 61.0 + 90.0*tan_f2 + GeoConsts_ep2_46*cos_f2 + 45.0*tan_f4
    - tan_f2*GeoConsts_ep2_252*cos_f2 + GeoConsts_ep4_n3*cos_f4 + GeoConsts_ep6_100*cos_f6
    - tan_f2*GeoConsts_ep4_66*cos_f4 + tan_f4*GeoConsts_ep2_n90*cos_f2 + GeoConsts_ep8_88*cos_f8;
	T12 = tan_f1 / (p1 * v5 * GeoConsts_k06_720) * temp12;
	temp13 = 1385.0 + 3633.0*tan_f2 + 4095.0*tan_f4 + 1575.0*tan_f6;
	T13 = tan_f1 / (p1 * v7 * GeoConsts_k08_40320) * temp13;
	T14 = 1.0 / (v * cos_f1 * GeoConsts_k0);
	T15 = 1.0 / (v3 * cos_f1 * GeoConsts_k03_6) * (1.0 + 2.0*tan_f2 + GeoConsts_ep2*cos_f2);
	temp16 = 5.0 + GeoConsts_ep2_6*cos_f2 + 28.0*tan_f2 - GeoConsts_ep4_3*cos_f4
    + tan_f2*GeoConsts_ep2_8*cos_f2 + 24.0*tan_f4 - GeoConsts_ep6_4*cos_f6
    + tan_f2*GeoConsts_ep4_4*cos_f4 + tan_f2*GeoConsts_ep6_24*cos_f6;
	T16 = 1.0 / (v5*cos_f1*GeoConsts_k05_120) * temp16;
	T17 = 1.0 / (v7*cos_f1*GeoConsts_k07_5040)
    * (61.0 + 662.0*tan_f2 + 1320.0*tan_f4 +720.0*tan_f6);
	
	// Calc Lat & Lon
    //
	dE2 = dE  * dE;
	dE3 = dE2 * dE;
	dE4 = dE3 * dE;
	dE5 = dE4 * dE;
	dE6 = dE5 * dE;
	dE7 = dE6 * dE;
	dE8 = dE7 * dE;
	
    // Latitude (in degrees)
    //
	lat = ftlat - dE2*T10 + dE4*T11 - dE6*T12 + dE8*T13;
	*latitude_deg = lat * GeoConsts_RADTODEG;
	
    // Longitude (in degrees)
    //
	lon = lon_cm + dE*T14 - dE3*T15 + dE5*T16 - dE7*T17;
	*longitude_deg = lon * GeoConsts_RADTODEG;
	
    return *status ;
}

void DSCML(double RLON, int NMAX, double *SINML, double *COSML){
	
	double A = sin(RLON);
	double B = cos(RLON);
	SINML[0] = A;
	COSML[0] = B;
	SINML[1] = 2.0*B*A;
	COSML[1] = 2.0*B*B-1.0;
	int M;
	for (M=3; M<=NMAX; M++){
		SINML[M-1]=2.0*B*SINML[M-2]-SINML[M-3];
		COSML[M-1]=2.0*B*COSML[M-2]-COSML[M-3];
	}
	return;
}

void LEGFDN(int M,double THETA,double *RLEG, double *DLEG,int NMX, int *IR,double *RLNN,int IFLAG,
            double *DRTS, double *DIRT){
	
    // THIS SUBROUTINE COMPUTES  ALL NORMALIZED LEGENDRE FUNCTION
    // IN "RLEG" AND THEIR DERIVATIVES IN "DLEG". ORDER IS ALWAYS
    // M , AND COLATITUDE IS ALWAYS THETA  (RADIANS). MAXIMUM DEG
    // IS  NMX  . ALL CALCULATIONS IN DOUBLE PRECISION.
    // IR  MUST BE SET TO ZERO BEFORE THE FIRST CALL TO THIS SUB.
    // THE DIMENSIONS OF ARRAYS  RLEG, DLEG, AND RLNN  MUST BE
    // AT LEAST EQUAL TO  NMX+1  .
	//
    // THIS PROGRAM DOES NOT COMPUTE DERIVATIVES AT THE POLES .
    //
	// IF    IFLAG = 1  , ONLY THE LEGENDRE FUNCTIONS ARE
    // COMPUTED.
    //
    // ORIGINAL PROGRAMMER :OSCAR L. COLOMBO, DEPT. OF GEODETIC SCIENCE
    // THE OHIO STATE UNIVERSITY, AUGUST 1980 . ******************
    //
	
	int NMX1 = NMX+1;
	int NMX2P = 2*NMX+1;
	int M1 = M+1;
	int M2 = M+2;
	int M3 = M+3;
	
	if (*IR!=1){
		*IR = 1;
		for (int N=1; N<=NMX2P; N++){
			DRTS[N-1] = sqrt((double)N);
			DIRT[N-1] = 1.0/DRTS[N-1];
		}
	}
	double COTHET = cos(THETA);
	double SITHET = sin(THETA);
    double SITHI = 0;
	
	if(IFLAG != 1 && THETA != 0.0)SITHI = 1.0/SITHET;
	
	// Compute the Legendre Functions
    //
	
	RLNN[0] = 1.0;
	RLNN[1] = SITHET*DRTS[2];
	int N,N1,N2;
	for(N1=3; N1<=M1; N1++){
		N = N1-1;
		N2 = 2*N;
		RLNN[N1-1] = DRTS[N2]*DIRT[N2-1]*SITHET*RLNN[N1-2];
	}
	
	if(M<=1){
		if(M==0){
			RLEG[0] = 1.0;
			RLEG[1] = COTHET*DRTS[2];
		}else{
			RLEG[1] = RLNN[1];
			RLEG[2] = DRTS[4]*COTHET*RLEG[1];
		}
	}
	RLEG[M1-1] = RLNN[M1-1];
	if(M2 <= NMX1){
		RLEG[M2-1] = DRTS[M1*2]*COTHET*RLEG[M1-1];
		if (M3 <= NMX1){
			for (N1=M3; N1<=NMX1; N1++){
				N = N1-1;
				if( !((M==0 && N<2) || (M==1 && N<3)) ){
					N2 = 2*N;
					RLEG[N1-1] = DRTS[N2]*DIRT[N+M-1]*DIRT[N-M-1]*
					(DRTS[N2-1-1]*COTHET*RLEG[N1-2] - DRTS[N+M-2]*DRTS[N-M-2]*DIRT[N2-4]*RLEG[N1-3]);
				}
			}
		}
	}
	if (IFLAG == 1) return;
	if (SITHET == 0.0){
		printf("Error: LEGDFN Doesn't compute drivatives at the poles\n");
		return;
	}
	
	// Now compute all the derivatives of the Legendre Functions
    //
	
	RLNN[0] = 0.0;
	double RLN, RLN1;
	RLN = RLNN[1];
	RLNN[1] = DRTS[2]*COTHET;
	for (N1=3; N1<=M1; N1++){
		N = N1-1;
		N2 = 2*N;
		RLN1 = RLNN[N1-1];
		RLNN[N1-1] = DRTS[N2]*DIRT[N2-1]*(SITHET*RLNN[N-1]+COTHET*RLN);
		RLN = RLN1;
	}
	DLEG[M1-1] = RLNN[M1-1];
	if (M2>NMX1) return;
	for(N1=M2; N1<=NMX1; N1++){
		N=N1-1;
		N2 = N*2;
		DLEG[N1-1] = SITHI*(N*RLEG[N1-1]*COTHET-DRTS[N-M-1]*DRTS[N+M-1]*DRTS[N2]*DIRT[N2-2]*RLEG[N-1]);
	}
	return;
}

void HUNDU(double *UNDU,int NMAX,double *P,const double *HC,const double *HS,
           double *SINML,double *COSML,double GR,double RE,const double *CCoef,
           const double *CS,double *HACO){
	// U IS THE GEOID UNDULATION FROM THE EGM96 POTENTIAL COEFFICIENT MODEL
    // INCLUDING THE HEIGHT ANOMALY TO GEOID UNDULATION CORRECTION TERM
    // AND A CORRECTION TERM TO HAVE THE UNDULATIONS REFER TO THE
    // WGS84 ELLIPSOID. THE GEOID UNDULATION UNIT IS METERS.
    //
	
	const double AE=GeoConsts_Re;
	
	double AR  = AE/RE;
	double ARN = AR;
	double AC  = 0.0;
	double A   = 0.0;
	int K      = 3;
	int N,M;
	double SUM,SUMC,SUM2;
	double TEMP,TEMPC;
	
	for(N=2; N<=NMAX; N++){
		ARN = ARN*AR;
		K=K+1;
		SUM=P[K-1]*HC[K-1];
		SUMC=P[K-1]*CCoef[K-1];
		SUM2=0.0;
		for(M=1; M<=N; M++){
			K=K+1;
			TEMPC=CCoef[K-1]*COSML[M-1]+CS[K-1]*SINML[M-1];
			TEMP=HC[K-1]*COSML[M-1]+HS[K-1]*SINML[M-1];
			SUMC=SUMC+P[K-1]*TEMPC;
			SUM=SUM+P[K-1]*TEMP;
		}
		AC=AC+SUMC;
		A=A+SUM*ARN;
	}
	AC=AC+CCoef[0]+P[1]*CCoef[1]+P[2]*(CCoef[2]*COSML[0]+CS[2]*SINML[0]);
	*HACO = AC/100.0;
	*UNDU=A*GeoConsts_GM/(GR*RE);
	// ADD HACO TO CONVERT HEIGHT ANOMALY ON THE ELLIPSOID TO THE UNDULATION
    // ADD -0.53M TO MAKE UNDULATION REFER TO THE WGS84 ELLIPSOID.
    //
	*UNDU = *UNDU + *HACO - 0.53;
	return;
}

