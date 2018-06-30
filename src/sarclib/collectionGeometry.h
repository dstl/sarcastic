/** @file********************************************************************
 *
 *       Module:    collectionGeometry.h
 *      Program:    sarclib
 *   Created by:    Darren on 24/07/2013.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      Function to calculate the collection geometry from a point within
 *      the synthetic aperture, defined by a pulse number, to a point on
 *      the ground.
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  14/03/2013
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

#ifndef sarclib_collectionGeometry_h
#define sarclib_collectionGeometry_h

#include "dataio_cphd.h"

enum lookDirection { LeftLook, RightLook }  ;

typedef struct {
    double range            ;    ///<  Range in metres
    double azimuthRad       ;    ///<  Azimuth angle in radians
    double grazingRad       ;    ///<  grazing angle in radians
    double squintRad        ;    ///<  squint angle in radians
    enum lookDirection look ;    ///<  look direction, left or right
    SPVector SRP            ;    ///<  SRP used for this collection
    SPVector grdPnt         ;    ///<  Ground point used for these calculations
    int pulse               ;    ///<  Pulse number used for these geometry calcs
    SPVector imgTL          ;    ///<  Image coordinates for top left in (lat,lon,Alt)
    SPVector imgTR          ;    ///<  Image coordinates for top right in (lat,lon,Alt)
    SPVector imgBR          ;    ///<  Image coordinates for bottom right in (lat,lon,Alt)
    SPVector imgBL          ;    ///<  Image coordinates for bottom left in (lat,lon,Alt)
    SPVector vel            ;    ///<  Velocity at synthetic aperture centre in m/2
} collectionGeom ;

/// Calculate the collection geometry between a location within the synthetic aperture and a point on the
/// ground. The location of the sensor is defined by the pulse index in 'pulse'. If 'pulse' is negative then
/// the function sets the sensor location to be the midpoint within the synthetic aperture. If 'pulse' is larger
/// than the number of pulses within the synthetic aperture then an error (INVALID_PULSE_INDEX) is thrown.
/// The ground location is specified by 'grdPnt'. If 'grdPnt' is (0,0,0) then the default stabilisation reference
/// point (SRP) for the pulse defined by 'pulse' is used for 'grdPnt'
///
int collectionGeometry(CPHDHeader *hdr, int pulse, SPVector grdPnt, collectionGeom *cGeom, SPStatus *status);

/// Print out a CPHD's information including its collection geometry
///
int printCPHDCollectionInfo(CPHDHeader *hdr1, SPStatus *status) ;

#endif
