//***************************************************************************
//
//  reflectPower.cl
//  GOSS
//
//  Created by Darren Muff on 02/08/2012.
//  Copyright (c) 2012 [dstl]. All rights reserved.
//
//
// Description:
//  OpenCl kernel code to calculate reflection power for an
//  array of rays back at the receiver
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
        VECT_ROTATEZ(inVect, angRads, out);
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

typedef struct Texture {
    float ka;    // Constant for ambient reflection
    float kd;    // Constant for diffuse scattering
    float ks;    // Constant for specular scattering
    float n;     // Shininess constant
} Texture;

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

typedef struct Triangle {
    int  triNum;    // Triangle ID
    double d;         // Constant of plane equation
    double nd_u;      // Normal.u / normal.k
    double nd_v;      // normal.v / normal.k
    int k;          // projection dimension
    double kbu;
    double kbv;
    double kbd;
    double kcu;
    double kcv;
    double kcd;
    int textureInd;
} Triangle;

typedef struct TriCoords {
    SPVector A ;      // Cartesian coordinates of triangle
    SPVector B ;
    SPVector Cc ;
} TriCoords ;
//+++++++++++++++++++++++++++ End of structures.cl +++++++++++++++++++++++++++++++++++++++++++++++


__kernel void reflectPower(__global Triangle * Triangles,       // Array of Triangles. Each triangle references a texture
                           __global Texture  * Textures,        // Array of Textures containing specular, diffuse adn shinyness constants
                           __global Hit *hits,                  // Array of hit locations to x-ref with triangles (and then Textures) for material props
                           const SPVector RxPos,                // Location of Receiver in x,y,z
                           const double GrxOverFourPi,          // Receiver antenna gain / 4Pi.
                           int nRays,
                           __global Ray *rays,                  // unit vector rays arriving at hitpoint
                           __global Ray *Rrays,                 // Reflection unit vector rays
                           __global Ray *Vrays,                 // Viewpoint unit vector. The vector direction to calculate power for (usually to receiver)
                           __global double *ranges,             // Range to receiver for each shadow ray (precalculated in shadowRay generation)
                           __global rangeAndPower *rnp           // Output array of ray power at, and range to reciever
                           )
{
    int ind, k, ku, kv ;
    Triangle T ;
    Texture tex ;
    SPVector N, R, L, V;
    double is, kd, ks, n, LdotN, RdotV, pow, rng;
    /* double id; */
    unsigned int modulo[5];
    
    ind = get_global_id(0) ;
    
    if (ind >=0 && ind < nRays ) {
        modulo[0] = 0; modulo[1] = 1; modulo[2] = 2; modulo[3] = 0; modulo[4]=1;

        T   = Triangles[hits[ind].trinum];
        tex = Textures[T.textureInd];
        
        // Create surface normal for taxture associated with this hitpoint
        //
        k  = T.k;
        ku = modulo[k+1];
        kv = modulo[k+2];
        
        N.cell[k]  = 1;
        N.cell[ku] = T.nd_u;
        N.cell[kv] = T.nd_v;
        VECT_NORM(N,N);
        
        // Set up Phong shading parameters for this texture
        //
        /*ia = 0 ;  */          // Ambient RF power. Assumed 0 for RF frequencies
        /*id = rays[ind].pow;*/ // Diffuse component incident energy on hitpoint
        is = rays[ind].pow ;    // Specular component of incident energy on hitpoint
        kd = tex.kd ;           // Diffuse constant of hitpoint's texture
        ks = tex.ks ;           // Specular constant of hitpoint's texture
        n  = tex.n ;            // Shinyness constant of hitpoints texture
//        printf("ray hit props are : is %f, kd %f, ks %f n %f\n",is,kd,ks,n);
        
        VECT_MINUS(rays[ind].dir, L); // L is the unit vector illumination ray (ray point to source of illumination
        R = Rrays[ind].dir;
        V = Vrays[ind].dir;
        
        LdotN = VECT_DOT(L,N);
        RdotV = VECT_DOT(R, V);
        
        LdotN = (LdotN < 0) ? -LdotN : LdotN;   // removes 'one-way mirror' effect
        RdotV = (RdotV < 0) ? 0 : RdotV;        // if ray bouncing away from radar then there is no specular component
        
        // This is the equation for the power but as there is no ambient Rf energy and is = id then we can refine it a bit
        //   pow = (ia * ka ) + (kd * LdotN * id) + ( ks * pow( RdotV , n) * is );
        //
        pow = is * ((kd * LdotN) + ( ks * pow( RdotV, n) )) ;
        
        // Now calculate the amount of power recieved at the receiver, taking into account
        // the inverse square law for return path and receiver antenna gain
        //
        rng = ranges[ind] ;
        rnp[ind].power = pow * GrxOverFourPi / (rng * rng) ;
        rnp[ind].range = (rng + Vrays[ind].len) * 0.5 ;
    }
    return ;
}