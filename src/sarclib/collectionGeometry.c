/***************************************************************************
 * 
 *           Module :  collectionGeometry.c
 *          Program :  sarclib
 *       Created by :  Darren Muff on Sat Jun 30 10:34:25 2018 +0100
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      Function to calculate the collection geometry from a point within
 *      the synthetic aperture, defined by a pulse number, to a point on 
 *      the ground.
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

#include "collectionGeometry.h"

// calculate collection geometry for a cphd file. 'pulse' is the pulse to
// calculate the geometry for. If 'pulse' is negative then the centre of
// the synthetic aperture is used.
// grdPnt is the ground point to calculate the geometry for (in ECEF coords). If
// this is 0,0,0 then the scene SRP is used.
// hdr must be already loaded using cphd_read_header
//
int collectionGeometry(CPHDHeader *hdr, int pulse, SPVector grdPnt, collectionGeom *cGeom, SPStatus *status){
    CHECK_STATUS( status );
    const double PI = SIPC_pi ;
    const double TWOPI = 2 * PI ;
    SPVector SACentPnt, Rcent, SACentPntMn, SACentPntMx, vel, slantPlaneNorm, rVec;
    SPVector rCross, rGrnd, e, n, Z, rGrndNorm, rCrossNorm, tmp, grdPntNorm, SRPDir, ptRange ;
    double SACentTimeBefore, SACentTimeAfter ;
    int SACentPulse ;
    
    SACentPulse = hdr->num_azi / 2 ;
    
    // If pulse is negative then use the centre pulse to calculate the geometry
    //
    if (pulse < 0 ) {
        pulse = hdr->num_azi / 2 ;
    }
    
    // If pulse is greater than number of pulses in collection then fail with error
    //
    if (pulse >= hdr->num_azi ) {
        if ( status->debug >=10 ){
            printf("ERROR : requested pulse index greater than number of pulses in collection for 'collectionGeometry'\n");
            status->status = INVALID_PULSE_INDEX ;
            return ( status->status ) ;
        }
    }
    // Otherwise calculate the geometry for the pulse specified
    //
    cGeom->pulse = pulse ;
    
    // Check Groundpoint requested
    //
    if ( grdPnt.x == 0 && grdPnt.y == 0 && grdPnt.z == 0 ) {
        grdPnt = hdr->pulses[SACentPulse].srp ;
    }
    
    cGeom->SRP    = hdr->pulses[SACentPulse].srp ;
    cGeom->grdPnt = grdPnt ;
    
    // Slant plane is defined as being from the centre of the Synthetic Aperture
    //
    SACentPnt = hdr->pulses[SACentPulse].sat_ps_rx ;
    
    // Calculate range vector from centre position
    //
    VECT_SUB(cGeom->SRP, SACentPnt, Rcent) ;
    
    // Calculate the velocity vector
    //
    SACentPntMn      = hdr->pulses[(SACentPulse)-1].sat_ps_rx ;
    SACentPntMx      = hdr->pulses[(SACentPulse)+1].sat_ps_rx ;
    SACentTimeBefore = hdr->pulses[(SACentPulse)-1].sat_rx_time ;
    SACentTimeAfter  = hdr->pulses[(SACentPulse)+1].sat_rx_time ;
    VECT_SUB(SACentPntMx, SACentPntMn, vel);
    VECT_SCMULT(vel, (1.0/(SACentTimeAfter-SACentTimeBefore)), vel);
    cGeom->vel = vel;
    
    // Define the slant plane from the centre position of the synthetic aperture and the velocity vectors
    //
    VECT_NORM(Rcent, Rcent);
    VECT_CROSS(Rcent, vel, slantPlaneNorm);
    VECT_NORM(slantPlaneNorm, slantPlaneNorm);
    VECT_NORM(cGeom->SRP, SRPDir) ;
    
    // Calculate look direction
    //
    cGeom->look = ((VECT_DOT(slantPlaneNorm, grdPnt) >= 0 ) ? RightLook : LeftLook );
    
    // Now calculate ground point (and pulse) specific geometry
    //
    VECT_SUB(grdPnt, hdr->pulses[pulse].sat_ps_rx, ptRange) ;
    VECT_NORM(ptRange, rVec) ;
    VECT_NORM(grdPnt, grdPntNorm);
    VECT_CROSS(rVec, grdPntNorm, rCross) ;
    VECT_CROSS(grdPntNorm, rCross, rGrnd) ;
    
    VECT_CREATE(0, 0, 1, Z);
    VECT_CROSS(Z, grdPnt, e) ;
    VECT_NORM(e, e) ;
    
    // n is the local north vector - ie the tangent to the earths surface - the
    // direction a compass would point if you were standing at the grdPnt
    //
    VECT_CROSS(grdPntNorm, e, n) ;
    VECT_NORM(n, n) ;
    VECT_NORM(rGrnd, rGrndNorm) ;
    VECT_NORM(rCross, rCrossNorm) ;
    
    // Calculate azimuth angle - ie bearing from target to sensor
    //
    cGeom->azimuthRad = acos( VECT_DOT(n, rGrndNorm) ) ;
    if ( VECT_DOT(rGrndNorm, e) < 0 )   cGeom->azimuthRad *= -1.0 ;
    cGeom->azimuthRad += PI ;
    if (cGeom->azimuthRad < 0 )         cGeom->azimuthRad += TWOPI ;
    if (cGeom->azimuthRad > TWOPI )     cGeom->azimuthRad -= TWOPI ;
    cGeom->range = VECT_MAG(ptRange) ;
    cGeom->grazingRad = acos(VECT_DOT(rGrndNorm, rVec)) ;
    VECT_CROSS(slantPlaneNorm, grdPntNorm, tmp);
    
    SPVector vdir ;
    VECT_NORM(vel, vdir);
    
    SPVector velXInGroundPlane ; VECT_CROSS(SRPDir, vdir, velXInGroundPlane) ;
    SPVector velInGroundPlane  ; VECT_CROSS(velXInGroundPlane, SRPDir, velInGroundPlane) ;
    VECT_NORM(velInGroundPlane, velInGroundPlane) ;
    cGeom->squintRad = asin(VECT_DOT(rGrndNorm, velInGroundPlane)) ;
    
    // Calculate image scene coordinates based upon illuminated scene from beamWidth
    // and assuming increasing range (shadows) are down
    //
    if(hdr->antenna_width_az == 0 || hdr->antenna_width_el == 0){
        printf("Warning : Antenna dimension not set. Required for scene size\n");
        SPVector NULLVEC; VECT_CREATE(0, 0, 0, NULLVEC);
        cGeom->imgTL = cGeom->imgTR = cGeom->imgBR = cGeom->imgBL = NULLVEC ;
    }else{
        double groundIlluminated_az = hdr->antenna_width_az * cGeom->range ;
        double groundIlluminated_ra = hdr->antenna_width_el * cGeom->range / sin(cGeom->grazingRad) ;
        
        SPVector rghtExtent, leftExtent, nearExtent, farExtent ;
        VECT_SCMULT(rCrossNorm, groundIlluminated_az / 2, rghtExtent);
        VECT_SCMULT(rCrossNorm, -1 * groundIlluminated_az / 2, leftExtent);
        VECT_SCMULT(rGrndNorm, groundIlluminated_ra / 2, farExtent);
        VECT_SCMULT(rGrndNorm, -1 * groundIlluminated_ra / 2, nearExtent);
        
        VECT_ADD(cGeom->SRP, nearExtent, tmp) ;
        VECT_ADD(tmp, rghtExtent, tmp);
        ecefToLatLonAlt(tmp.x, tmp.y, tmp.z, Ellipsoid, &(cGeom->imgTL), status) ;
        
        VECT_ADD(cGeom->SRP, nearExtent, tmp) ;
        VECT_ADD(tmp, leftExtent, tmp);
        ecefToLatLonAlt(tmp.x, tmp.y, tmp.z, Ellipsoid, &(cGeom->imgTR), status) ;
        
        VECT_ADD(cGeom->SRP, farExtent, tmp) ;
        VECT_ADD(tmp, leftExtent, tmp);
        ecefToLatLonAlt(tmp.x, tmp.y, tmp.z, Ellipsoid, &(cGeom->imgBR), status) ;
        
        VECT_ADD(cGeom->SRP, farExtent, tmp) ;
        VECT_ADD(tmp, rghtExtent, tmp);
        ecefToLatLonAlt(tmp.x, tmp.y, tmp.z, Ellipsoid, &(cGeom->imgBL), status) ;
    }
    
    return ( status->status ) ;

}
