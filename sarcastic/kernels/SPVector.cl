//***************************************************************************
//
//  SPVector.cl
//  SARCASTIC
//
//  Created by Darren Muff on 02/08/2012.
//  Copyright (c) 2012 [dstl]. All rights reserved.
//
//
// Description:
//
//
// CLASSIFICATION        :  UNCLASSIFIED
// Date of CLASSN        :  26/02/2014
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// THE SOFTWARE IN ITS ENTIRETY OR ANY SUBSTANTIAL PORTION SHOULD NOT BE
// USED AS A WHOLE OR COMPONENT PART OF DERIVATIVE SOFTWARE THAT IS BEING
// SOLD TO THE GOVERNMENT OF THE UNITED KINGDOM OF GREAT BRITAIN AND NORTHERN
// IRELAND.
//
//***************************************************************************

#ifndef SPVECTOR_CL
#define SPVECTOR_CL
#define PROGEPSILON ((double)1.0e-8)
/// defines the structure for a vector
///
typedef union {
    struct { double x, y, z; };
    struct { double r, g, b; };
    struct { double cell[3]; };
    struct { double R, theta, phi; };
    struct { double lat, lon, alt; };
    
} SPVector ;

/// Create a vector
///
#define VECT_CREATE(a,b,c,outVect) {    \
outVect.x = a;                      \
outVect.y = b;                      \
outVect.z = c;                      \
}

/// Multiply a vector by a constant
///
#define VECT_SCMULT(inVect,inScal,outVect)  {   \
outVect.x = (inVect.x)*inScal;              \
outVect.y = (inVect.y)*inScal;              \
outVect.z = (inVect.z)*inScal;              \
}

/// Subtract two vectors
///
#define VECT_SUB(aVect,bVect,outVect) {                 \
outVect.x = aVect.x - bVect.x;                      \
outVect.y = aVect.y - bVect.y;                      \
outVect.z = aVect.z - bVect.z;                      \
}

/// add two vectors
///
#define VECT_ADD(aVect,bVect,outVect) {                 \
outVect.x = aVect.x + bVect.x;                      \
outVect.y = aVect.y + bVect.y;                      \
outVect.z = aVect.z + bVect.z;                      \
}

/// Find the modulus (or magnitdue) of a vector
///
#define VECT_MOD(aVect) (                                       \
sqrt (aVect.x*aVect.x + aVect.y*aVect.y + aVect.z*aVect.z)  \
)
#define VECT_MAG(aVect)(                                        \
sqrt (aVect.x*aVect.x + aVect.y*aVect.y + aVect.z*aVect.z)  \
)

/// Negate a vector
///
#define VECT_MINUS(inVect,outVect) {                            \
VECT_CREATE(-inVect.x, -inVect.y, -inVect.z, outVect)       \
}

/// Find the dot product of two vectors
///
#define VECT_DOT(aVect,bVect) (aVect.x*bVect.x + aVect.y*bVect.y + aVect.z*bVect.z)

/// find a cross product of two vectors
///
#define VECT_CROSS(aVect,bVect,outVect) {                       \
SPVector ___tmp ;                                           \
___tmp.x = aVect.y*bVect.z - aVect.z*bVect.y;               \
___tmp.y = aVect.z*bVect.x - aVect.x*bVect.z;               \
___tmp.z = aVect.x*bVect.y - aVect.y*bVect.x;               \
outVect = ___tmp ;                                          \
}

/// Find the equivalent unit vector v
///
#define VECT_UNIT(aVect,outVect) {                              \
double vect_unit_tmp = 1.0 / VECT_MOD(aVect);               \
outVect.x = aVect.x*vect_unit_tmp;                          \
outVect.y = aVect.y*vect_unit_tmp;                          \
outVect.z = aVect.z*vect_unit_tmp;                          \
}
#define VECT_NORM(aVect,outVect) {                              \
double vect_unit_tmp = 1.0 / VECT_MOD(aVect);               \
outVect.x = aVect.x*vect_unit_tmp;                          \
outVect.y = aVect.y*vect_unit_tmp;                          \
outVect.z = aVect.z*vect_unit_tmp;                          \
}

/// project a vector along a normal
///
#define VECT_PROJ(vect, norm, outVect) {                        \
SPVector ___tmp;                                            \
double dot = VECT_DOT(vect, norm);                          \
___tmp.x = vect.x - dot*norm.x;                             \
___tmp.y = vect.y - dot*norm.y;                             \
___tmp.z = vect.z - dot*norm.z;                             \
VECT_UNIT(___tmp, outVect);                                 \
}

/// Rotate a vector around the X-axis
///
#define VECT_ROTATEX(inVect, angRads, outVect) {                                    \
double c   = cos(angRads) ;                                                     \
double s   = sin(angRads) ;                                                     \
double __y = inVect.y*c-inVect.z*s ;                                            \
double __z = inVect.y*s+inVect.z*c ;                                            \
VECT_CREATE(inVect.x, __y, __z, outVect);                                       \
}

/// Rotate a vector around the Y-axis
///
#define VECT_ROTATEY(inVect, angRads, outVect) {                                    \
double c   = cos(angRads) ;                                                     \
double s   = sin(angRads) ;                                                     \
double __x = inVect.x*c+inVect.z*s ;                                            \
double __z = inVect.z*c-inVect.x*s ;                                            \
VECT_CREATE(__x, inVect.y, __z, outVect);                                       \
}

/// Rotate a vector around the Z-axis
///
#define VECT_ROTATEZ(inVect, angRads, outVect) {                                    \
double c   = cos(angRads) ;                                                     \
double s   = sin(angRads) ;                                                     \
double __x =  inVect.x*c-inVect.y*s ;                                           \
double __y =  inVect.x*s+inVect.y*c ;                                           \
VECT_CREATE(__x, __y, inVect.z, outVect) ;                                      \
}

/// Rotate a vector around another vector
///
void vectRotateAxis(SPVector inVect, SPVector axisVect, double angRads, SPVector *outVect);

// Rotate a vector around another vector
//
void vectRotateAxis(SPVector inVect, SPVector axisVect, double angRads, SPVector *outVect){
    SPVector out ;
    if( fabs(axisVect.x) < PROGEPSILON && fabs(axisVect.y) < PROGEPSILON ){
        VECT_ROTATEZ(inVect, angRads, out) ;
    }else if(fabs(axisVect.x) < PROGEPSILON && fabs(axisVect.z) < PROGEPSILON){
        VECT_ROTATEY(inVect, angRads, out);
    }else if(fabs(axisVect.y) < PROGEPSILON && fabs(axisVect.z) < PROGEPSILON){
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

#endif
