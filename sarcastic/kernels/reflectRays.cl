//***************************************************************************
//
//  reflectRays.cl
//  SARCASTIC
//
//  Created by Darren Muff on 02/08/2012.
//  Copyright (c) 2012 [dstl]. All rights reserved.
//
//
// Description:
//  OpenCl kernel code to reflect an array of rays. Inputs are
//  the rays and an array of hitpoints (which contains the triangle id
//  for this hit) and an array of triangles
//
//
// CLASSIFICATION        :  UNCLASSIFIED
// Date of CLASSN        :  02/08/2012
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
#pragma OPENCL EXTENSION cl_khr_fp64: enable
//#include "materialProperties.h"
#define MATBYTES 128

typedef struct scatProps {
    char   matname[MATBYTES] ;
    float  corlen      ;
    float  roughness   ;
    float  resistivity ;
    float  specular    ;
    float  diffuse     ;
    float  shinyness   ;
} scatProps ;

#define NMATERIALS 9
static constant scatProps materialProperties[NMATERIALS] = {
    //   Name          corrLen      Roughness   Resistivity Specular    Diffuse     Shinyness
    {"MATERIAL",    100.0,      0.0,        0.0,        1.0,        0.0,        50.0        },
    {"ASPHALT",     0.5,        0.005,      1.0e18,     0.8,        0.2,        30.0        },
    {"BRICK",       0.1,        0.001,      1.0e18,     0.7,        0.3,        20.0        },
    {"CONCRETE",    0.2,        0.01,       120.0,      0.3,        0.7,        10.0        },
    {"METAL",       100.0,      0.0,        1.0e-8,     1.0,        0.0,        50.0        },
    {"ROOFING",     0.1,        0.1,        1.0e18,     0.6,        0.4,        40.0        },
    {"VEGETATION",  0.01,       0.1,        2000.0,     0.2,        0.8,        5.0         },
    {"WATER",       0.01,       0.1,        2.0e1,      1.0,        0.0,        50.0        },
    {"WOOD",        0.1,        0.001,      1.0e14,     0.6,        0.4,        10.0        }
} ;

// #include "SPVector.cl"
//+++++++++++++++++++++++++++ Start of SPVector.cl +++++++++++++++++++++++++++++++++++++++++++++++

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
//+++++++++++++++++++++++++++ End of SPVector.cl +++++++++++++++++++++++++++++++++++++++++++++++

//#include "structures.cl"
//+++++++++++++++++++++++++++ Start of structures.cl +++++++++++++++++++++++++++++++++++++++++++++++
typedef struct cplxf {
    float r;
    float i;
} cplxf;

typedef struct AABB {
    SPVector AA;
    SPVector BB;
} AABB;

typedef union {
    struct KdTreeLeaf {
        unsigned int flagDimAndOffset;
        // bits 0..1        : splitting dimension
        // bits 2..30       : offset bits
        // bits 31 (sign)   : flag whether node is a leaf
        float splitPosition;
        int Ropes[6];
        AABB aabb;
    } leaf;
    
    struct KdTreeBranch {
        unsigned int flagDimAndOffset;
        // bits 0..30       : offset to first child
        // bit 31 (sign)    : flag whether node is a leaf
        float splitPosition;
        int Ropes[6];
        AABB aabb;
    } branch;
    
} KdData ;
typedef struct Ray {
    SPVector org;    // Origin
    SPVector dir;    // Direction
    double   pow;    // Power for this ray
    double   len;    // Distance travelled to this ray's origin from transmission
    SPVector pol ;   // unit vector of direction of E field of ray
} Ray;

typedef struct rangeAndPower {
    double range ;
    double power ;
} rangeAndPower ;

typedef struct Hit {
    double dist;
    int trinum;
    double u;
    double v;
} Hit;

typedef struct ATS {
    int  triNum;    // Triangle ID
    double d;       // Constant of plane equation
    double nd_u;    // Normal.u / normal.k
    double nd_v;    // normal.v / normal.k
    int k;          // projection dimension
    double kbu;
    double kbv;
    double kbd;
    double kcu;
    double kcv;
    double kcd;
    int    matInd;  // Material Index
} ATS;

typedef struct TriCoords {
    SPVector A ;      // Cartesian coordinates of triangle
    SPVector B ;
    SPVector Cc ;
} TriCoords ;
//+++++++++++++++++++++++++++ End of structures.cl +++++++++++++++++++++++++++++++++++++++++++++++

__kernel void reflect(__global ATS * accelTriangles,
                      __global Ray *rays,
                      __global Hit *hits,
                      __global Ray *reflectedRays,
                      int nRays
                      )
{
    int ind, t, k, ku, kv ;
    unsigned int modulo[5];
    ATS T ;
    SPVector N, v, R, I, hp, perpol;
    float ks;
 
    modulo[0] = 0; modulo[1] = 1; modulo[2] = 2; modulo[3] = 0; modulo[4]=1;
    ind = get_global_id(0) ;
    
    if (ind >=0 && ind < nRays ) {
        
        t = hits[ind].trinum;
        T = accelTriangles[t];
        k = T.k;
        ku = modulo[k+1];
        kv = modulo[k+2];
        N.cell[k] = 1;
        N.cell[ku] = T.nd_u;
        N.cell[kv] = T.nd_v;
        VECT_NORM(N,N);
        
        I = rays[ind].dir;
        if(VECT_DOT(I,N)>0){
            VECT_MINUS(N,N);
        }
        VECT_SCMULT(N, 2.0 * VECT_DOT(I, N), v);
        VECT_SUB(I, v, R);
        
        VECT_SCMULT(rays[ind].dir, hits[ind].dist, hp);
        
        VECT_ADD(rays[ind].org, hp, reflectedRays[ind].org);
        reflectedRays[ind].dir = R ;
        reflectedRays[ind].len = rays[ind].len + hits[ind].dist ;
        
        // Now calculate forward scattered ray power by only considering the
        // specular component of the reflective surface texture.
        // (Diffuse and shinyness components are only taken into account on rays
        // returning to the sensor from a visible hit point
        // Note that this is because we use PO for last bounce but phong shding for all others
        //
        ks =  materialProperties[T.matInd].specular ;
        // power at reflected point is Pt*Gtx / 4 PI R^2
        //
        double islf = (4.0 * 3.1415926536 * hits[ind].dist * hits[ind].dist);
        islf = (islf < 1 ) ? 1 : islf ;
        reflectedRays[ind].pow = rays[ind].pow * ks ;
        
        // Calculate the polarisation of the reflected ray based upon the polarisation
        // of the incident ray
        //
        VECT_CROSS(rays[ind].pol, N, perpol);
        VECT_CROSS(perpol, reflectedRays[ind].dir, reflectedRays[ind].pol);
        VECT_NORM(reflectedRays[ind].pol,reflectedRays[ind].pol);

    }
    
    return ;
}

