/** @file********************************************************************
 *
 *       Module:    GeoVector.h
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

#ifndef SILib2_GeoVector_h
#define SILib2_GeoVector_h

enum GeoAltReference {Ellipsoid,Geoid};
enum GeoCoordsystems {WGS84, UTM, MGRS};

/// Function to convert any angle to a correct latitude
/// between -90 (south pole) and +90 (north pole)
///
double fixlat(double lat);

/// function to convert any angle to a correct longitude
/// between  -180 deg and +180 deg
///
double fixlon(double lon);

/// Function to convert a lat,lon.alt to an ECEF coordinate
/// The alt can be specified as being ref to 'Ellipsoid' or 'Geoid' in
/// the 4th parameter
///
SPStatus latLonAltToEcef(double lattitude_deg, double longitude_deg, double alt, enum GeoAltReference GeoidOrEllipsoid, SPVector *ecef, SPStatus *status);

/// Function to convert ECEF parameters (X,Y,Z) to a lat,lon,alt SPVector
/// The altitude in the return vector is referenced to the ellipsoid if the 4th parameter
/// is set t be 'Ellipsoid' or referenced to the geoid if the 4th parameter is set
/// to 'Geoid'
///
SPStatus ecefToLatLonAlt(double X, double Y, double Z, enum GeoAltReference GeoidOrEllipsoid, SPVector *latlonalt, SPStatus *status);

/// Undulation is the term refering to the difference in height in metres between
/// the Earth's gravitional potential model and the Earth's ellipsoid model for a
/// given lat and long.
/// This function calculates the geoid undulation from the EGM96 potential coefficient model
/// including the height anomaly to geoid undulation correction term and a
/// correction term to have the undulations refer to the WGS84 ellipsoid. The
/// geoid undulation unit is metres.
///
double undulation(double lattitude_deg, double longitude_deg);

/// convert a latitude and longitude to a Universal Transverse Mercator (UTM)
/// easting and northing (and the associated UTMZone)
///
SPStatus LLtoUTM(double latitude_deg, double longitude_deg, double *UTMEasting, double *UTMNorthing, int *UTMZone, SPStatus *status);

/// convert a Universal Transverse Mercator (UTM) easting and northing (and the associated UTMZone)
/// to a latitude and longitude
///
SPStatus UTMtoLL(double UTMEasting, double UTMNorthing, int UTMZone, double *latitude_deg, double *longitude_deg, SPStatus *status);

/// Function that returns the height difference between the geoid and ellipsoid
/// for a given lat, lon
///
double GeoidHeightToEllipsoidHeight(double GeoidHeight, double lattitude_deg, double longitude_deg) ;

/// Function that returns the height difference between the ellipsoid and geoid
/// for a given lat, lon
///
double EllipsoidHeightToGeoidHeight(double ellipsoidHeight, double lattitude_deg, double longitude_deg) ;

#endif
