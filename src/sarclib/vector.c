/***************************************************************************
 * 
 *           Module :  vector.c
 *          Program :  sarclib
 *       Created by :  Darren Muff on 22/07/2013
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      vector calcultions that are not efficient to run as macros
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

#include "sarclib.h"
#define PROGEPSILON ((double)1.0e-8)

// Rotate a vector around another vector
//
void vectRotateAxis(SPVector inVect, SPVector axisVect, double angRads, SPVector *outVect){
    SPVector out ;
    if( fabs(axisVect.x) < PROGEPSILON && fabs(axisVect.y) < PROGEPSILON ){
        if(axisVect.z < 0) angRads = -angRads ;
        VECT_ROTATEZ(inVect, angRads, out) ;
    }else if(fabs(axisVect.x) < PROGEPSILON && fabs(axisVect.z) < PROGEPSILON){
        if(axisVect.y < 0) angRads = -angRads ;
        VECT_ROTATEY(inVect, angRads, out);
    }else if(fabs(axisVect.y) < PROGEPSILON && fabs(axisVect.z) < PROGEPSILON){
        if(axisVect.x < 0) angRads = -angRads ;
        VECT_ROTATEX(inVect, angRads, out);
    }else{
        SPVector ansa,ansb,ansc,ansd;
        double thetaz = atan2(axisVect.y,axisVect.x);
        double thetay = atan2( sqrt( (axisVect.x * axisVect.x) + (axisVect.y * axisVect.y)) , axisVect.z );
        
        // First rotate around the z axis
        //
        VECT_ROTATEZ(inVect, -thetaz, ansa);
        
        // Now rotate around y axis
        //
        VECT_ROTATEY(ansa, -thetay, ansb);
        
        // now rotate around z by the required angle theta
        //
        VECT_ROTATEZ(ansb, angRads, ansc);
        
        // Now add on the rotation axis around y and z to get back to the original reference frame
        //
        VECT_ROTATEY(ansc, thetay, ansd);
        VECT_ROTATEZ(ansd, thetaz, out);
    }
    outVect->x = out.x ;
    outVect->y = out.y ;
    outVect->z = out.z ;
}