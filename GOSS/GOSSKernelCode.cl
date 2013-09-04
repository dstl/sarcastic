//***************************************************************************
//
//  kernelCode.c
//  GOSS
//
//  Created by Darren Muff on 02/08/2012.
//  Copyright (c) 2012 [dstl]. All rights reserved.
//
//
// Description:
//     This file needs to be self contained as it will become the
//     OpenCL / Cuda kernel
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

#if defined(cl_khr_fp64)  // Khronos extension available?
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)  // AMD extension available?
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

#define MAXBOUNCES 8
#define PI 3.1415927
#define CC 299792458.0
#define NOINTERSECTION -1
#define true 1
#define false 0
#define KDT_ISLEAF(n)    (n->leaf.flagDimAndOffset & (unsigned int)(1<<31))
#define KDT_DIMENSION(n) (n->leaf.flagDimAndOffset & 0x3)
#define KDT_OFFSET(n)    ((n->leaf.flagDimAndOffset & (0x7FFFFFFC))>>2)
#define NILROPE ((int) -666666 )
#define EPSILON ((double) 0.000001)
#define RAD2DEG ((double)(57.2957795))  // Radians to Degrees conversion
#define TRUE 1
#define FALSE 0
#define HIDDENSCATTERERS 0
#define RXLNA 1000        // Receiver Low-Noise amplifier in dB (30dB)
#define RXImp 50          // Receiver Impedence

typedef int OutCode;

#define INSIDEAABB  0   // 000000
#define LEFTAABB    1   // 000001
#define RIGHTAABB   2   // 000010
#define NEARAABB    4   // 000100
#define FARAABB     8   // 001000
#define TOPAABB     16  // 010000
#define BOTTOMAABB  32  // 100000
#define NOINTERSECTION -1

void ClipToRect(float *x0, float *y0, float *x1, float *y1, float xmin, float xmax, float ymin, float ymax, int *status);

// defines the structure for a vector
//
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
    outVect.x = aVect.y*bVect.z - aVect.z*bVect.y;              \
    outVect.y = aVect.z*bVect.x - aVect.x*bVect.z;              \
    outVect.z = aVect.x*bVect.y - aVect.y*bVect.x;              \
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
    SPVector __tmp;                                             \
    double __dot = VECT_DOT(vect, norm);                        \
    __tmp.x = vect.x - __dot*norm.x;                            \
    __tmp.y = vect.y - __dot*norm.y;                            \
    __tmp.z = vect.z - __dot*norm.z;                            \
    VECT_UNIT(tmp, outVect);                                    \
}

/// Rotate a vector around the X-axis
///
#define VECT_ROTATEX(inVect, angRads, outVect) {                                    \
    double c = cos(angRads) ;                                                       \
    double s = sin(angRads) ;                                                       \
    VECT_CREATE(inVect.x, inVect.y*c-inVect.z*s, inVect.y*s+inVect.z*c, outVect)    \
}

/// Rotate a vector around the Y-axis
///
#define VECT_ROTATEY(inVect, angRads, outVect) {                                    \
    double c = cos(angRads) ;                                                       \
    double s = sin(angRads) ;                                                       \
    VECT_CREATE(inVect.x*c+inVect.z*s, inVect.y, inVect.z*c-inVect.x*s, outVect)    \
}

/// Rotate a vector around the Z-axis
///
#define VECT_ROTATEZ(inVect, angRads, outVect) {                                    \
    double c = cos(angRads) ;                                                       \
    double s = sin(angRads) ;                                                       \
    VECT_CREATE(inVect.x*c-inVect.y*s, inVect.x*s+inVect.y*c, inVect.z, outVect)    \
}

/// Rotate a vector around another vector
///
void vectRotateAxis(SPVector inVect, SPVector axisVect, double angRads, SPVector *outVect, int debug);

typedef struct Ray {
    SPVector org;    // Origin
    SPVector dir;    // Direction
    float   san;    // Solid Angle for Ray
} Ray;

typedef struct Hit {
    double dist;
    int trinum;
    double u;
    double v;
} Hit;

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

typedef struct rangeAndPower {
    double range ;
    double power ;
} rangeAndPower ;

typedef struct beamParams {
    int nAzBeam ;                // Number of azimuth slices in beam
    int nElBeam ;                // Number of elevation slices in beam
    SPVector RxPos ;             // Receiver position for this pulse
    SPVector TxPos ;             // Transmitter position for this pulse
    double dAz ;                 // Azimuth slice in radians
    double dEl ;                 // Elevation slice in radians
    double raySolidAng ;         // Solid angle of a single ray
    double TxPowPerRay ;         // Transmitter power per ray
    AABB SceneBoundingBox ;      // Scene bounding box
    double Aeff ;                // The effective area of the Receive Antenna
    int bounceToShow;            // Which bounce to output
} beamParams ;

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


// forward declaration of functions
//
int occluded(SPVector point, SPVector dir,AABB SceneBoundingBox, __global KdData * KdTree, __global int *triangleListData, __global int * triangleListPtrs, __global Triangle * Triangles);
Hit StacklessTraverse(Ray ray, AABB SceneBoundingBox, __global KdData * KdTree, __global int *triangleListData, __global int * triangleListPtrs, __global Triangle * Triangles, int debug);
double reflectPower(double rayPower, SPVector L,    SPVector R, SPVector V, Hit h, __global Triangle * Triangles, __global Texture * textureData);
Ray reflect(Ray ray, Hit h, __global Triangle * Triangles);
int clipToAABB(AABB boundingBox, SPVector *lineStart, SPVector *lineEnd);
void ClipToBox(SPVector *p0, SPVector *p1, SPVector min, SPVector max, int *status);
void Intersect(__global Triangle *tri, Ray *ray, Hit *hit);
OutCode ComputeOutCode(SPVector p, SPVector min, SPVector max);
SPVector hitPoint (Hit h, Ray r);


__kernel void rayTraceBeam (const int nAzBeam,              // Number of azimuth slices in beam
                            const int nElBeam,              // Number of elevation slices in beam
                            const SPVector RxPos,           // Receiver position for this pulse
                            const SPVector TxPos,           // Transmitter position for this pulse
                            const double dAz,               // Azimuth slice in radians
                            const double dEl,               // Elevation slice in radians
                            const double raySolidAng,       // Solid angle of a single ray
                            const double TxPowPerRay,       // Transmitter power per ray
                            const AABB SceneBoundingBox,    // Scene bounding box
                            const double Aeff,              // The effective area of the Receive Antenna
                            const int bounceToShow,         // Which bounce to print out
                            __global Triangle * Triangles,  // array - triangle data
                            __global Texture * textureData, // array - texture data
                           __global KdData * KdTree,        // array containing KdTree
                           __global int *triangleListData,  // array of triangle indices into Triangles
                           __global int *triangleListPtrs,  // array of pointers into triangleListData for each node
                           __global rangeAndPower *rnp,      // Array storing ranges and powers (dim: nAzbeam x nElBeam x MAXBOUNCES )
                            const double testdbl
                            )
{
    int xId = get_global_id(0);
    int yId = get_global_id(1);
    
    SPVector aimDir, elAxis, azAxis, zHat, azRayDir, elRayDir, v_tmp ;
    double rayPow, totalDist=0, RCSArea, beamWidth, beamHeight, thetaAz, thetaEl ;
    double range, refPow, refPowatRx, rangeOfreturn ;
    int bounces=0,debug=0,ind;
    Ray r;
    Hit h;
        
    if (xId >=0 && xId < nAzBeam && yId >=0 && yId < nElBeam) {

        // DEBUG a single thread by setting values below and uncommenting
        //
        if(xId==0 && yId == 0){debug = 10;}else{debug=0;}
        
        if(debug){
            printf("%d,%d : Device parameters :\n",xId,yId);
            printf("%d,%d : nAzBeam : %d\n",xId,yId,nAzBeam);
            printf("%d,%d : nElBeam : %d\n",xId,yId,nElBeam);
            printf("%d,%d : RxPos : %1.9e,%1.9e,%1.9e\n",xId,yId,RxPos.x, RxPos.y,RxPos.z);
            printf("%d,%d : TxPos : %1.9e,%1.9e,%1.9e\n",xId,yId,TxPos.x, TxPos.x,TxPos.z);
            printf("%d,%d : dAz : %1.9e\n",xId,yId,dAz);
            printf("%d,%d : dEl : %1.9e\n",xId,yId,dEl);
            printf("%d,%d : raySolidAng : %1.9e\n",xId,yId,raySolidAng);
            printf("%d,%d : TxPowPerRay : %e\n",xId,yId,TxPowPerRay);
            printf("%d,%d : SceneBoundingBox : %f,%f,%f-%f,%f,%f\n",xId,yId,SceneBoundingBox.AA.x,SceneBoundingBox.AA.y,SceneBoundingBox.AA.z,SceneBoundingBox.BB.x,SceneBoundingBox.BB.y,SceneBoundingBox.BB.z);
            printf("%d,%d : Aeff : %e\n",xId,yId,Aeff);
            printf("%d,%d : bounceToShow : %d\n",xId,yId,bounceToShow);
            printf("%d,%d : testdbl is %1.9e\n",xId,yId,testdbl);
        }
        
        beamWidth  = nAzBeam * dAz;
        beamHeight = nElBeam * dEl;
        
        if(debug>=10)printf("rxp : %f,%f,%f\n",RxPos.x,RxPos.y,RxPos.z);

        VECT_MINUS(RxPos, aimDir) ;
        if(debug>=10)printf("aimDir : %f,%f,%f\n",aimDir.x,aimDir.y,aimDir.z);

        VECT_CREATE(0,0,1.,zHat) ;
        VECT_CROSS(aimDir,zHat,elAxis) ;
        VECT_CROSS(elAxis,aimDir,azAxis) ;
        
        // Angle of this azimuth slice from the aimDir in radians
        //
        thetaAz = (yId * dAz) - (beamWidth/2) + (dAz/2);
        
        // Rotate aimDir around azimuth beam axis of rotation to get the direction of this Ray
        //
        if(debug>=10)printf("inVect at call: %f,%f,%f\n",aimDir.x,aimDir.y,aimDir.z);
        vectRotateAxis(aimDir, azAxis, thetaAz, &azRayDir, 0);

        // Angle of this elevation slice from the aimDir in radians
        //
        thetaEl = (xId * dEl) - (beamHeight/2) + (dEl/2);
        
        // For every elevation slice within pulse rotate the ray in the elevation direction.
        //
        vectRotateAxis(azRayDir, elAxis, thetaEl, &elRayDir,debug);
        if(debug>=10)printf("azRayDir : %f,%f,%f\n",azRayDir.x,azRayDir.y,azRayDir.z);

        
        r.org  = TxPos ;
        r.san  = raySolidAng ;
        rayPow = TxPowPerRay ;
        if(debug>=10)printf("elRayDir : %f,%f,%f\n",elRayDir.x,elRayDir.y,elRayDir.z);

        VECT_NORM(elRayDir, r.dir) ;
        if(debug>=10)printf("r.dir : %f,%f,%f\n",r.dir.x,r.dir.y,r.dir.z);

        
        while (bounces < MAXBOUNCES){
            if(debug>=10)printf("================start of stacklesstraverse===============\n");
            h = StacklessTraverse(r, SceneBoundingBox, KdTree, triangleListData, triangleListPtrs, Triangles,debug);
            if(debug>=10)printf("================= end of stacklesstraverse===============\n");
            
            if(debug>=10)printf("h: trinum = %d, dist = %f\n",h.trinum,h.dist);
            if( h.trinum == NOINTERSECTION){
                break ;
            }
            // if here then we hit something
            //
            
            // For power calculations we need the effective area of the triangle being hit
            //
            totalDist = totalDist + h.dist; // total dist to hitpoint
            RCSArea = r.san * totalDist * totalDist ;
            rayPow = RCSArea * TxPowPerRay / (4 * PI * totalDist * totalDist);  // pow at hitpoint
            // rayPow = rayPow * 1000000000000.0; // to cover gaps between rays
            // Above value should be the 'collander fill factor' ie the difference in solid angles between
            // the sum of all the rays and the part of teh beam being used.
            //
            Ray reflected = reflect(r, h, Triangles);
            
            SPVector returnVect, returnVectDir ;
            VECT_SUB(RxPos, reflected.org, returnVect);
            VECT_NORM(returnVect, returnVectDir) ;
            
            if (! occluded(reflected.org, returnVectDir, SceneBoundingBox, KdTree, triangleListData, triangleListPtrs, Triangles) )
                // If not occluded - ie there is a path from the found intersection point back to the SAR receiver
                // Calculate power and range and return the values
                //
            {
                if(bounces == bounceToShow){
                    printf("%f,%f,%f\n",reflected.org.x,reflected.org.y,reflected.org.z);
                }
                if(debug>=10)printf("%f,%f,%f\n",reflected.org.x,reflected.org.y,reflected.org.z);
                
                range = VECT_MAG( returnVect );
                VECT_MINUS(r.dir, v_tmp) ;
                refPow = reflectPower(rayPow, v_tmp, reflected.dir, returnVectDir, h,Triangles, textureData);
                
                refPowatRx = refPow * Aeff / (4.0 * PI * range * range);
                rangeOfreturn = (totalDist + range)*0.5;
                
                ind = (bounces * nAzBeam * nElBeam) + (yId * nAzBeam) + xId ;
                rnp[ind].range = rangeOfreturn ;
                rnp[ind].power = refPowatRx ;
                
            }
            r = reflected;
            bounces++;
        }
        
    }
}

Hit StacklessTraverse(Ray ray, AABB SceneBoundingBox, __global KdData * KdTree, __global int *triangleListData, __global int *triangleListPtrs,
                      __global Triangle * Triangles, int debug){
    
    int i,trisInLeaf;
    float t_entry = 0;
    float t_exit  = VECT_MAG( ray.org ) + 100;
    SPVector volumeEntry, volumeExit, PEntry, v_tmp, dirInverse;
;
    Hit hit;
    float t1,t2,xinv,yinv,zinv;
    
    hit.dist = 10e6;
    hit.trinum = NOINTERSECTION;
    
    // Create an infinite line from the ray
    //
    volumeEntry  = ray.org;
    VECT_SCMULT(ray.dir, t_exit, v_tmp);
    VECT_ADD(ray.org, v_tmp, volumeExit);

    if(debug>=10){
        printf("volumeEntry : %f,%f,%f\n",volumeEntry.x,volumeEntry.y,volumeEntry.z);
        printf("volumeExit : %f,%f,%f\n",volumeExit.x,volumeExit.y,volumeExit.z);
        printf("ray.dir: %f,%f,%f\n",ray.dir.x,ray.dir.y,ray.dir.z);
        printf("t_exit: %f\n",t_exit);
        SPVector at;
        at.x = ray.dir.x * t_exit;
        at.y = ray.dir.y * t_exit;
        at.z = ray.dir.z * t_exit;
        printf("at: %f,%f,%f\n",at.x,at.y,at.z);
        printf("v_tmp: %f,%f,%f\n",v_tmp.x,v_tmp.y,v_tmp.z);
    }
    
    // Calculate the ray segment within the scene volume
    // Do this early to reduce calcs for ray misses
    //
    if(!clipToAABB(SceneBoundingBox, &volumeEntry, &volumeExit)){
        return hit; // ray misses bounding volume
    }
    
    // Calc inverse direction so only multiplies (not divides) in loop (cheaper)
    //
    xinv = (ray.dir.x == 0) ? 0 : 1./ray.dir.x;
    yinv = (ray.dir.y == 0) ? 0 : 1./ray.dir.y;
    zinv = (ray.dir.z == 0) ? 0 : 1./ray.dir.z;
    VECT_CREATE(xinv, yinv, zinv, dirInverse);
    
    // initialise the root of the tree
    //
    __global KdData * node  = &(KdTree[0]);
    
    int dimToUse=0;
    while (dirInverse.cell[dimToUse] == 0)dimToUse++;
    t_entry = 0.0;
    t_exit  = (volumeExit.cell[dimToUse]   - volumeEntry.cell[dimToUse])* dirInverse.cell[dimToUse];
    
    while (t_entry <= t_exit ) {
        
        VECT_SCMULT(ray.dir, t_entry, v_tmp);
        VECT_ADD(volumeEntry, v_tmp, PEntry);
        
        if(debug>=10){
            printf("PEntry:%f,%f,%f\n",PEntry.x,PEntry.y,PEntry.z);
            printf("volumeEntry : %f,%f,%f\n",volumeEntry.x,volumeEntry.y,volumeEntry.z);
            printf("ray.dir: %f,%f,%f\n",ray.dir.x,ray.dir.y,ray.dir.z);
            printf("t_entry : %f\n",t_entry);
        }
        
        while (!KDT_ISLEAF(node)){  // Branch node
            
            if (PEntry.cell[KDT_DIMENSION(node)] < (node->branch.splitPosition-EPSILON)) {
                if(debug>=10)printf("[%f < %f (k:%d)], next node is %d\n",PEntry.cell[KDT_DIMENSION(node)],node->branch.splitPosition-EPSILON,KDT_DIMENSION(node),KDT_OFFSET(node));
                node = &(KdTree[KDT_OFFSET(node)]);
            }else if (PEntry.cell[KDT_DIMENSION(node)] > (node->branch.splitPosition+EPSILON)) {
                if(debug>=10)printf("[%f > %f (k:%d)], next node is %d\n",PEntry.cell[KDT_DIMENSION(node)],node->branch.splitPosition+EPSILON,KDT_DIMENSION(node),KDT_OFFSET(node)+1);
                node = &(KdTree[KDT_OFFSET(node)+1]);
            }else{
                // PEntry is on splitposition
                // Two different situations here:
                // 1. PEntry is on the split plane having arrived from neg side.
                // 2. PEntry is on the split plane having arrived from pos side.
                // in Case 1. next node is node
                // in Case 2. next node is node+1
                if (ray.org.cell[KDT_DIMENSION(node)] < node->branch.splitPosition-EPSILON){
                    if(debug>=10)printf("o[%f < %f (k:%d)] next node is %d\n",ray.org.cell[KDT_DIMENSION(node)],node->branch.splitPosition,KDT_DIMENSION(node),KDT_OFFSET(node));
                    node = &(KdTree[KDT_OFFSET(node)]);
                } else if (ray.org.cell[KDT_DIMENSION(node)] > node->branch.splitPosition+EPSILON){
                    if(debug>=10)printf("o[%f > %f (k:%d)] next node is %d\n",ray.org.cell[KDT_DIMENSION(node)],node->branch.splitPosition,KDT_DIMENSION(node),KDT_OFFSET(node)+1);
                    node = &(KdTree[KDT_OFFSET(node)+1]);
                } else {
                    // ray origin on split plane. Determine next node by direction of ray
                    //
                    if(ray.dir.cell[KDT_DIMENSION(node)] > 0 ){
                        if(debug>=10)printf("|[] next node is %d, dir is %f\n",KDT_OFFSET(node)+1,ray.dir.cell[KDT_DIMENSION(node)]);
                        node = &(KdTree[KDT_OFFSET(node)+1]);
                    }else {
                        // Includes situation where origin is on split position and ray is travelling parallel to split position
                        //
                        if(debug>=10)printf("|[] next node is %d dir is %f\n",KDT_OFFSET(node),ray.dir.cell[KDT_DIMENSION(node)]);
                        node = &(KdTree[KDT_OFFSET(node)]);
                    }
                }
            }
        }
        // Have a leaf now
        //
        trisInLeaf = triangleListData[triangleListPtrs[KDT_OFFSET(node)]];
        if(debug>=10)printf("Leaf found : %d triangles\n",trisInLeaf);
        for (i=0; i<trisInLeaf; i++){
            __global Triangle * tri = &(Triangles[triangleListData[triangleListPtrs[KDT_OFFSET(node)]+i+1]]);
            Intersect(tri, &ray, &hit);
        }
        if((hit.trinum != NOINTERSECTION) && (hit.dist > 0.01)){
            return hit;
        }
        
        if(debug>=10)printf("No Intersection found; Chasing rope...\n");
        // If ray doesnt intersect triangle in this leaf then propagate the
        // ray to an adjacent node using this leaf's Ropes. (rather than popping
        // back up the tree)
        // Set t_entry to be the ray distance to the adjacent AABB
        //
        
        float t_max = 10e6;
        float tpos;
        int ropeInd, ropeIndSide, ropeIndOff;
        for(int i=0; i<3; i++){
            if (dirInverse.cell[i] != 0){
                t1 = (node->branch.aabb.AA.cell[i] - volumeEntry.cell[i]) * dirInverse.cell[i];
                t2 = (node->branch.aabb.BB.cell[i] - volumeEntry.cell[i]) * dirInverse.cell[i];
                
                if(t2 > t1 && t2 >= 0){
                    tpos = t2;
                    ropeIndSide = 1;
                }else if( t1 > t2 && t1 >= 0){
                    tpos = t1;
                    ropeIndSide = 0;
                }else if (t2 == t1){
                    // AABB is planar so select rope based upon direction of ray
                    //
                    tpos = t1;
                    ropeIndSide = (dirInverse.cell[i] < 0 ) ? 0 : 1;
                }
                if(tpos < t_max){
                    t_max = tpos;
                    ropeInd = i;
                    ropeIndOff = ropeIndSide;
                }
            }
        }
        t_entry = t_max;
        
        // Follow the rope for this side to jump across to adjacent node
        //
        
        if ( node->branch.Ropes[2*ropeInd+ropeIndOff] == NILROPE ) {
            hit.trinum = NOINTERSECTION;
            return hit ;
        }
        node = &(KdTree[node->branch.Ropes[2*ropeInd+ropeIndOff]]);
        
    }
    return hit ;
}

int clipToAABB(AABB boundingBox, SPVector *lineStart, SPVector *lineEnd){
    
    int status=0;
    ClipToBox(lineStart, lineEnd, boundingBox.AA, boundingBox.BB, &status);
    if (status == NOINTERSECTION) return 0; // No intersect in any dim means no intersection with volume
    
    return 1;
    
}

// Cohen–Sutherland clipping algorithm clips a line from
// P0 = (x0, y0) to P1 = (x1, y1) against a rectangle with
// diagonal from (xmin, ymin) to (xmax, ymax).
//
void ClipToBox(SPVector *p0, SPVector *p1, SPVector min, SPVector max, int *status)
{
    // compute outcodes for P0, P1, and whatever point lies outside the clip rectangle
    //
    OutCode outcode0 = ComputeOutCode(*p0, min, max);
    OutCode outcode1 = ComputeOutCode(*p1, min, max);
    int accept = FALSE;
    
    while (TRUE) {
        if (!(outcode0 | outcode1)) { // Bitwise OR is 0. Trivially accept and get out of loop
            accept = TRUE;
            break;
        } else if (outcode0 & outcode1) { // Bitwise AND is not 0. Trivially reject and get out of loop
            break;
        } else {
            // failed both tests, so calculate the line segment to clip
            // from an outside point to an intersection with clip edge
            //
            double x, y, z, m;
            
            // At least one endpoint is outside the clip rectangle; pick it.
            //
            OutCode outcodeOut = outcode0 ? outcode0 : outcode1;
            
            // Now find the intersection point;
            //
            if (outcodeOut & FARAABB) {             // Point is beyond the box (y-axis)
                m = (max.y - p0->y) / (p1->y - p0->y);
                x = p0->x + m * (p1->x - p0->x) ;
                y = max.y;
                z = p0->z + m * (p1->z - p0->z) ;
            }else if (outcodeOut & NEARAABB){     // Point is before the box (y-axis)
                m = (min.y - p0->y) / (p1->y - p0->y);
                x = p0->x + m * (p1->x - p0->x) ;
                y = min.y;
                z = p0->z + m * (p1->z - p0->z) ;
            }else if (outcodeOut & LEFTAABB){       // Point is left of box (x-axis)
                m = (min.x - p0->x) / (p1->x - p0->x);
                x = min.x;
                y = p0->y + m * (p1->y - p0->y);
                z = p0->z + m * (p1->z - p0->z);
            }else if (outcodeOut & RIGHTAABB){      // Point is right of box (x-axis)
                m = (max.x - p0->x) / (p1->x - p0->x);
                x = max.x;
                y = p0->y + m * (p1->y - p0->y);
                z = p0->z + m * (p1->z - p0->z);
            }else if (outcodeOut & TOPAABB){        // Point is above box (z-axis)
                m = (max.z - p0->z) / (p1->z - p0->z);
                x = p0->x + m * (p1->x - p0->x);
                y = p0->y + m * (p1->y - p0->y);
                z = max.z;
            }else{                              // Point is below box (z-axis)
                m = (min.z - p0->z) / (p1->z - p0->z);
                x = p0->x + m * (p1->x - p0->x);
                y = p0->y + m * (p1->y - p0->y);
                z = min.z;
            }
            
            // Now we move outside point to intersection point to clip
            // and get ready for next pass.
            //
            if (outcodeOut == outcode0) {
                p0->x = x;
                p0->y = y;
                p0->z = z;
                outcode0 = ComputeOutCode(*p0, min, max);
            } else {
                p1->x = x;
                p1->y = y;
                p1->z = z;
                outcode1 = ComputeOutCode(*p1, min, max);
            }
        }
    }
    if (accept) {
        return ;
    }else{
        *status = NOINTERSECTION ;
        return ;
    }
}

// Compute the bit code for a point (x, y) using the clip rectangle
// bounded diagonally by (xmin, ymin), and (xmax, ymax)
//
OutCode ComputeOutCode(SPVector p, SPVector min, SPVector max)
{
    OutCode code;
    
    code = INSIDEAABB;          // initialised as being inside of clip window
    
    if (p.x < min.x)            // to the left of box
        code |= LEFTAABB;
    else if (p.x > max.x)       // to the right of box
        code |= RIGHTAABB;
    if (p.y < min.y)            // nearer that the box
        code |= NEARAABB;
    else if (p.y > max.y)       // farther than the box
        code |= FARAABB;
    if (p.z < min.z)            // Below bottom of box
        code |= BOTTOMAABB;
    else if(p.z > max.z)        // Above top of box
        code |= TOPAABB;
    
    return code;
}


void Intersect(__global Triangle *tri, Ray *ray, Hit *hit){
    unsigned int modulo[5];
    modulo[0] = 0; modulo[1] = 1; modulo[2] = 2; modulo[3] = 0; modulo[4]=1;
    
    int ku = modulo[tri->k+1];
    int kv = modulo[tri->k+2];
    
    const double nd = 1.0/(ray->dir.cell[tri->k] + ((tri->nd_u) * ray->dir.cell[ ku ]) + ((tri->nd_v) * ray->dir.cell[ kv ]) );
    const double thit = (tri->d - ray->org.cell[tri->k] - tri->nd_u * ray->org.cell[ ku ] - tri->nd_v * ray->org.cell[ kv ]) * nd;
    
    // check for valid distance.
    //
    if ( !(hit->dist > thit && thit >  EPSILON  ) ) return;
    
    // compute hitpoint positions on uv plane
    //
    const double hu = (ray->org.cell[ku] + thit * ray->dir.cell[ ku ]);
    const double hv = (ray->org.cell[kv] + thit * ray->dir.cell[ kv ]);
    
    // check first barycentric coordinate
    //
    const double beta = (hu * tri->kbu + hv * tri->kbv + tri->kbd);
    if (beta < 0.0f) return ;
    
    // check second barycentric coordinate￼
    //
    const double gamma = (hu * tri->kcu + hv * tri->kcv + tri->kcd);
    if (gamma < 0.0f) return;
    
    // check third barycentric coordinate
    //
    if (beta+gamma > 1.0f) return ;
    
    // have a valid hitpoint here. store it.
    //
    hit->dist = thit;
    hit->trinum = tri->triNum;
    hit->u = beta;
    hit->v = gamma;
    return ;
}


Ray reflect(Ray ray, Hit h, __global Triangle * Triangles){
    unsigned int modulo[5];
    modulo[0] = 0; modulo[1] = 1; modulo[2] = 2; modulo[3] = 0; modulo[4]=1;
    int t = h.trinum;
    Triangle T = Triangles[t];
    int k = T.k;
    int ku = modulo[k+1];
    int kv = modulo[k+2];
    SPVector N;
    N.cell[k] = 1;
    N.cell[ku] = T.nd_u;
    N.cell[kv] = T.nd_v;
    
    SPVector R, v_tmp;
    SPVector I = ray.dir;
    
    VECT_SCMULT(N, (2.0 * (VECT_DOT(I,N))), v_tmp );
    VECT_SUB(I,v_tmp, R) ;
    Ray reflected;
    reflected.org = hitPoint(h, ray);
    reflected.dir = R;
    // For now just set the reflected solida angle to be the same as that of
    // the inbound ray. Later this can be adjusted to calculate the SA based upon
    // the dimensions of the reflecting triangle
    //
    reflected.san = ray.san ;
    
    return reflected;
}

SPVector hitPoint (Hit h, Ray r){
    SPVector ans;
    double x = r.org.x + h.dist * r.dir.x;
    double y = r.org.y + h.dist * r.dir.y;
    double z = r.org.z + h.dist * r.dir.z;
    VECT_CREATE(x,y,z, ans);
    return ans;
}

int occluded(SPVector point, SPVector dir, AABB SceneBoundingBox,
             __global KdData * KdTree,
             __global int * triangleListData,
             __global int * triangleListPtrs,
             __global Triangle * Triangles){
    Hit h;
    Ray r;
    r.org = point;
    r.dir = dir;
    r.san = 0;
    h = StacklessTraverse(r,SceneBoundingBox,KdTree,triangleListData,triangleListPtrs,Triangles,0);
    if (h.trinum != NOINTERSECTION) return 1;
    return 0;
}

double reflectPower(double rayPower, SPVector L, SPVector R, SPVector V, Hit h, __global Triangle * Triangles, __global Texture * textureData){
    // Phong equation
	// Ip = ka*ia + kd(L.N)id + ks((R.V)^n)is
	// Ip is the power at the intersection point reflected in direction V
    // ka = ambient reflection constant of material
    // ia = incident ambient power
	// kd = diffuse reflection constant of the material
    // L = unit vector pointing to light source (ie opposite direction of srcray)
	// N is the surface normal unit vector
    // id = incident diffuse power reflected
    // is = incident specular power
	// ks is the specular reflection constant of the material
	// R is the Reflection unit vector
	// V is the viewpoint unit vector (ie vector pointing to RXposition)
	// n is the shininess constant.
	// Note - it is assumed that there is no ambient power at these frequencies hence ia is zero
    
    Triangle T = Triangles[h.trinum];
    Texture texture = textureData[T.textureInd];
    int k = T.k;
    unsigned int modulo[5];
    modulo[0] = 0; modulo[1] = 1; modulo[2] = 2; modulo[3] = 0; modulo[4]=1;
    int ku = modulo[k+1];
    int kv = modulo[k+2];
    SPVector N;
    N.cell[k] = 1;
    N.cell[ku] = T.nd_u;
    N.cell[kv] = T.nd_v;
    double id = rayPower;
    double is = rayPower;
    double ia = 0;
    double ka = 0.0;
    double kd = 0.0 ;  //texture.kd;
    double ks = 1.0 ;  //texture.ks;
    double n  = 50.0 ; // texture.n;
    double LdotN = VECT_DOT(L, N);
    VECT_NORM(R,R) ;
    VECT_NORM(V,V) ;
    double RdotV = VECT_DOT(R, V);
    
    LdotN = (LdotN < 0) ? -LdotN : LdotN;   // removes 'one-way mirror' effect
    RdotV = (RdotV < 0) ? 0 : RdotV;        // if ray bouncing away from radar then there is no specular component
    double ans = (ia * ka ) +  (kd * LdotN * id) + ( ks * pow( RdotV , n) * is );    
    return ans;
}

// Rotate a vector around another vector
//
void vectRotateAxis(SPVector inVect, SPVector axisVect, double angRads, SPVector *outVect, int debug){
    if(debug>=10)printf("inVect on entry: %f,%f,%f\n", inVect.x,inVect.y,inVect.z);
    SPVector out;
    if(axisVect.x == 0 && axisVect.y == 0 ){
        VECT_ROTATEZ(inVect, angRads, out) ;
    }else if(axisVect.x == 0 && axisVect.z == 0){
        VECT_ROTATEY(inVect, angRads, out);
    }else if(axisVect.y == 0 && axisVect.z == 0){
        VECT_ROTATEZ(inVect, angRads, out);
    }else{
        SPVector ansa,ansb,ansc,ansd;
        double thetaz = atan(axisVect.y/axisVect.x);
        double thetay = atan( sqrt( (axisVect.x * axisVect.x) + (axisVect.y * axisVect.y)) / axisVect.z );
        if(debug>=10){
            printf("thetaz : %f\n",thetaz);
            printf("thetay : %f\n",thetaz);
            printf("inVect : %f,%f,%f \n",inVect.x,inVect.y,inVect.z);
            printf("axisVect : %f,%f,%f\n",axisVect.x,axisVect.y,axisVect.z);
            printf("angRads : %f\n", angRads);
        }
        
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
    
    outVect->x = out.x;
    outVect->y = out.y;
    outVect->z = out.z;
    
}