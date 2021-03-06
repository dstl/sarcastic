/***************************************************************************
 * 
 *           Module :  latlon.h
 *          Program :  sarclib
 *       Created by :  Matt Nottingham on 04/01/2006.
 *                     Darren Muff on 21/07/2013. 
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      Functions for converting between ECEF & Lat/Longs.
 *
 *      Lat/Longs are stored in an SPVector type where:
 *          x = latitude
 *          y = longitude
 *          z = alt
 *
 * 
 *   (c) Crown Copyright 2018 Defence Science and Technology Laboratory
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software")
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
 ***************************************************************************/

#ifndef sarclib_LATLON_H__
#define sarclib_LATLON_H__

/// Convert ECEF ====> WGS84 Lat/long
///
SPStatus ecef2latlon(SPVector * ecef, SPVector * latlon, SPStatus * status); 

/// Convert WGS84 Lat/long ====> ECEF
///
SPStatus latlon2ecef(SPVector * latlon, SPVector * ecef, SPStatus * status);

/// Function that returns the height difference between the geoid and ellipsoid
/// for a given lat, lon
///
SPStatus geoid2ellipsoid(SPVector *latlongeoidalt, SPVector *latlonellipsoidalt, SPStatus *status);

/// Function that returns the height difference between the ellipsoid and geoid
/// for a given lat, lon
///
SPStatus ellipsoid2geoid(SPVector *latlonellipsoidalt, SPVector *latlongeoidalt, SPStatus *status);

#endif
