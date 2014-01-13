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
//     If compiled with DEBUG defined then various levels of debug information
//      are displayed.
//          DEBUG >= 10 top level output from rayTraceBeam
//          DEBUG >= 15 output from stacklessTraverse() function
//          DEBUG >= 20 output from intersect() function
//          DEBUg >= 25 output from reflectPower() function
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
#define RXLNA 1000000000000        // Receiver Low-Noise amplifier in dB (30dB)
#define MAXTRAVERSAL 1000          // Maximum kdTree traversal steps before we blow up
typedef int OutCode;

#define INSIDEAABB  0   // 000000
#define LEFTAABB    1   // 000001
#define RIGHTAABB   2   // 000010
#define NEARAABB    4   // 000100
#define FARAABB     8   // 001000
#define TOPAABB     16  // 010000
#define BOTTOMAABB  32  // 100000
#define NOINTERSECTION -1

#ifdef DEBUG
#ifndef DBGX
#define DBGX 0
#endif
#ifndef DBGY
#define DBGY 0
#endif
#endif

void ClipToRect(float *x0, float *y0, float *x1, float *y1, float xmin, float xmax, float ymin, float ymax, int *status);

typedef union {
    struct { double x, y, z; };
    struct { double r, g, b; };
    struct { double cell[3]; };
    struct { double R, theta, phi; };
    struct { double lat, lon, alt; };
    
} VectorH ;

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
    VectorH AA;
    VectorH BB;
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
    VectorH org;    // Origin
    VectorH dir;    // Direction
    float   san;    // Solid Angle for Ray
} Ray;

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
    VectorH A ;      // Cartesian coordinates of triangle
    VectorH B ;
    VectorH Cc ;
} TriCoords ;

typedef struct RadarParams {
    int nPulses;                    // radar parameter - pulses in output pulses
    int nSamps;                     // radar parameter - samples in output pulses
    float oneOverLambda;            // radar parameter - one over wavelength
    double TxPowPerRay;             // radar parameter - transmitter power per ray in beam
    double  TxRaySolidAngle;         // radar parameter - Solid angle of a single ray cast by transmitter
    float chirpRate;                // radar parameter - chirp rate Hz/Sec
    float ADRate;                   // radar parameter - ad sample rate
    float dAz;                      // radar parameter - az slice in radians
    float dEl;                      // radar parameter - el slice in radians
    float StartFrequency;           // radar parameter - startFrequency of transmitted chirp (Fcent - Tp*gamma/2)
    int nAzBeam;                    // num of rays in azimuth in radar beam
    int nElBeam;                    // num of rays in elevation in radar beam
    int bounceToShow ;              // Which bounce to display -1= none, 0=1st bounce etc...
    double Aeff;                     // Radar Parameter - The effective area of the Receive Antenna
} RadarParams;

typedef struct rangeAndPower {
    double range ;
    double power ;
} rangeAndPower ;

VectorH hitPoint (Hit h, Ray r);
VectorH vectCreate(double a, double b, double c);
VectorH vectSub(VectorH a, VectorH b);
VectorH vectAdd(VectorH a, VectorH b);
VectorH vectMult(VectorH v, double sc);
VectorH vectMinus(VectorH v);
double vectMag(VectorH v);
double vectDot(VectorH a, VectorH b);
VectorH vectCross(VectorH a, VectorH b);
VectorH vectRotatex(VectorH v,double a);
VectorH vectRotatey(VectorH v,double a);
VectorH vectRotatez(VectorH v,double a);
VectorH vectRotateAxis(VectorH v, VectorH axis, double theta);
VectorH vectNorm(VectorH v);

int occluded(VectorH point, VectorH dir,AABB SceneBoundingBox, __global KdData * KdTree, __global int *triangleListData, __global int * triangleListPtrs, __global Triangle * Triangles);
Hit StacklessTraverse(Ray ray, AABB SceneBoundingBox, __global KdData * KdTree, __global int *triangleListData, __global int * triangleListPtrs, __global Triangle * Triangles, int debug, int rayx, int rayy, int pulseIndex);
double reflectPower(double rayPower, VectorH L, VectorH R, VectorH V, Hit h, __global Triangle * Triangles, __global Texture * textureData, int debug);
Ray reflect(Ray ray, Hit h, __global Triangle * Triangles);
int clipToAABB(AABB boundingBox, VectorH *lineStart, VectorH *lineEnd);
void ClipToBox(VectorH *p0, VectorH *p1, VectorH min, VectorH max, int *status);
void Intersect(__global Triangle *tri, Ray *ray, Hit *hit, int debug);
OutCode ComputeOutCode(VectorH p, VectorH min, VectorH max);

VectorH vectCreate(double a, double b, double c){
    VectorH v;
    v.x = a; v.y = b; v.z = c;
    return(v);
}
VectorH vectSub(VectorH a, VectorH b){ return vectCreate(a.x-b.x, a.y-b.y, a.z-b.z) ;}
VectorH vectAdd(VectorH a, VectorH b){ return vectCreate(a.x+b.x, a.y+b.y, a.z+b.z) ;}
VectorH vectMult(VectorH v, double sc){ return vectCreate(v.x*sc, v.y*sc, v.z*sc) ;}
double vectDot(VectorH a, VectorH b){return(a.x*b.x + a.y*b.y + a.z*b.z);}
VectorH vectMinus(VectorH v){ return vectCreate(-v.x, -v.y, -v.z) ; }
double vectMag(VectorH v){ return sqrt(v.x*v.x+v.y*v.y+v.z*v.z) ; }
VectorH vectCross(VectorH a, VectorH b){ return vectCreate(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x); }
VectorH vectNorm(VectorH v){ double l = 1.0f / vectMag(v); return vectCreate(v.x*l,v.y*l,v.z*l); }
VectorH vectRotatex(VectorH v,double a){
    double c = cos(a); double s = sin(a);
    return vectCreate(v.x, v.y*c-v.z*s, s*v.y + c*v.z);
}
VectorH vectRotatey(VectorH v,double a){
    double c = cos(a); double s = sin(a);
    return vectCreate(c*v.x+s*v.z, v.y, c*v.z - v.x*s);
}
VectorH vectRotatez(VectorH v,double angle){
    double c = cos(angle); double s = sin(angle);
    return vectCreate(c*v.x-s*v.y, s*v.x+c*v.y, v.z);
}
VectorH vectRotateAxis(VectorH v, VectorH axis, double theta){
    VectorH ans;
    if(axis.x == 0 && axis.y == 0 ){
        ans = vectRotatez(v,theta);
    }else if(axis.x == 0 && axis.z == 0){
        ans = vectRotatey(v,theta);
    }else if(axis.y == 0 && axis.z == 0){
        ans = vectRotatex(v,theta);
    }else{
        double thetaz = atan(axis.y/axis.x);
        double thetay = atan( sqrt( (axis.x * axis.x) + (axis.y * axis.y)) / axis.z );
                
        // First rotate around the z axis
        ans = vectRotatez(v,-thetaz);
        
        // Now rotate around y axis
        ans = vectRotatey(ans,-thetay);
        
        // now rotate around z by the required angle theta
        ans = vectRotatez(ans,theta);
        
        // Now add on the rotation axis around y and z to get back to the original reference frame
        ans = vectRotatey(ans,thetay);
        ans = vectRotatez(ans,thetaz);
    }
    return(ans);
}

__kernel void rayTraceBeam (const int nAzBeam,              // Number of azimuth slices in beam
                            const int nElBeam,              // Number of elevation slices in beam
                            const VectorH RxPos,            // Receiver position for this pulse
                            const VectorH TxPos,            // Transmitter position for this pulse
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
                           __global rangeAndPower *rnp,     // Array storing ranges and powers (dim: nAzbeam x nElBeam x MAXBOUNCES )
                            const int pulseIndex            // For debugging - Index of this pulse
                            )
{
    int xId = get_global_id(0);
    int yId = get_global_id(1);
    
    double beamHeight, beamWidth, thetaAz, thetaEl, rayPow, totalDist, RCSArea;
    VectorH aimDir, elAxis, azAxis, azRayDir, elRayDir, rxp;
    VectorH srcPt, dstPt ;
    double dist;
    int debug = 0, bounces;
    Ray r;
    Hit h;

    if (xId >=0 && xId < nAzBeam && yId >=0 && yId < nElBeam) {

#ifdef DEBUG
        if( xId == DBGX && yId == DBGY ){debug = DEBUG;}else{debug=0;}
        if(debug>=10)printf("+++++++++++++++++++++++++++++++++++++++\n");
#endif
        beamWidth  = nAzBeam * dAz;
        beamHeight = nElBeam * dEl;
        rxp        = RxPos;
        aimDir     = vectMinus(rxp);
        elAxis     = vectCross(aimDir,vectCreate(0, 0, 1));
        azAxis     = vectCross(elAxis,aimDir);
        thetaAz    = (xId * dAz) - (beamWidth/2) + (dAz/2);
        azRayDir   = vectRotateAxis(aimDir, azAxis, thetaAz);
        thetaEl    = (yId * dEl) - (beamHeight/2) + (dEl/2);
        elRayDir   = vectRotateAxis(azRayDir, elAxis, thetaEl);
        r.org      = TxPos;
        r.dir      = vectNorm (  elRayDir );
        r.san      = raySolidAng ;
        rayPow     = TxPowPerRay;
        totalDist  = 0.0;
        bounces    = 0;
                
        while (bounces < MAXBOUNCES){
#ifdef DEBUG
            if(debug>=10)printf("================start of stacklesstraverse===============\n");
#endif
            h = StacklessTraverse(r, SceneBoundingBox, KdTree, triangleListData, triangleListPtrs, Triangles,debug,xId,yId, pulseIndex);
#ifdef DEBUG
            if(debug>=10)printf("================= end of stacklesstraverse===============\n");
            if(debug>=10)printf("h: trinum = %d, dist = %f\n",h.trinum,h.dist);
#endif
            if( h.trinum == NOINTERSECTION){
                break ;
            }
            // if here then we hit something
            //
            
            // For power calculations we need the effective area of the triangle being hit
            //
            Ray reflected = reflect(r, h, Triangles);
            srcPt = r.org;
            dstPt = reflected.org;
            dist  =  vectMag( vectSub(dstPt,srcPt) );
#ifdef DEBUG
            if(debug>=10)printf("dist = %f\n",h.trinum,dist);
#endif
            totalDist = totalDist + dist ;
            RCSArea = r.san * totalDist * totalDist ;
            rayPow = RCSArea * TxPowPerRay / (4 * PI * totalDist * totalDist);  // pow at hitpoint
            
            VectorH returnVect = vectSub(rxp, reflected.org );
            VectorH returnVectDir = vectNorm(returnVect);

#ifdef DEBUG
            if(debug>=10){
                printf("%f,%f,%f\n",reflected.org.x,reflected.org.y,reflected.org.z);
                printf("ray dir: %f,%f,%f\n", reflected.dir.x, reflected.dir.y, reflected.dir.z);
            }
#endif
            if(bounces == bounceToShow){
                printf("%f,%f,%f\n",reflected.org.x,reflected.org.y,reflected.org.z);
            }
            
            if (! occluded(reflected.org, returnVectDir, SceneBoundingBox, KdTree, triangleListData, triangleListPtrs, Triangles) )
                // If not occluded - ie there is a path from the found intersection point back to the SAR receiver
                // Calculate power and range and mix the signal back into the receiver
                //
            {
                
                double range = vectMag( returnVect );
                double refPow = reflectPower(rayPow, vectMinus(r.dir), reflected.dir, returnVectDir, h,Triangles, textureData, debug);
                double refPowatRx = RXLNA*refPow * Aeff / (4.0 * PI * range * range);
                double rangeOfreturn = (totalDist + range)*0.5;
                
                rnp[((yId * nAzBeam) + xId)+(nAzBeam*nElBeam*bounces)].range = rangeOfreturn ;
                rnp[((yId * nAzBeam) + xId)+(nAzBeam*nElBeam*bounces)].power = refPowatRx ;
            }
            r = reflected;
            bounces++;
        }
    }
}


Hit StacklessTraverse(Ray ray, AABB SceneBoundingBox, __global KdData * KdTree, __global int *triangleListData, __global int *triangleListPtrs,
                      __global Triangle * Triangles, int debug, int rayx, int rayy, int pulseIndex){
    
    int i,trisInLeaf;
    float t_entry = 0;
    float t_exit  = vectMag( ray.org ) + 100;
    VectorH volumeEntry, volumeExit, PEntry, hp;
    Hit hit;
    float t1,t2,xinv,yinv,zinv;
    
    hit.dist = 10e6;
    hit.trinum = NOINTERSECTION;
    
    // Create an infinite line from the ray
    //
    volumeEntry  = ray.org;
    volumeExit   = vectAdd(ray.org, vectMult(ray.dir, t_exit ));
    
#ifdef DEBUG
    if(debug>=15){
        printf("volumeEntry : %f,%f,%f\n",volumeEntry.x,volumeEntry.y,volumeEntry.z);
        printf("volumeExit  : %f,%f,%f\n",volumeExit.x,volumeExit.y,volumeExit.z);
        printf("ray.dir     : %f,%f,%f\n",ray.dir.x,ray.dir.y,ray.dir.z);
    }
#endif
    
    // Calculate the ray segment within the scene volume
    // Do this early to reduce calcs for ray misses
    if(!clipToAABB(SceneBoundingBox, &volumeEntry, &volumeExit)){
        return hit; // ray misses bounding volume
    }
    
    // Calc inverse direction so only multiplies (not divides) in loop (cheaper)
    xinv = (ray.dir.x == 0) ? 0 : 1./ray.dir.x;
    yinv = (ray.dir.y == 0) ? 0 : 1./ray.dir.y;
    zinv = (ray.dir.z == 0) ? 0 : 1./ray.dir.z;
    VectorH dirInverse = vectCreate(xinv, yinv, zinv);
    
    // initialise the root of the tree
    __global KdData * node  = &(KdTree[0]);
    
    int dimToUse=0;
    while (dirInverse.cell[dimToUse] == 0)dimToUse++;
    t_entry = 0.0;
    t_exit  = (volumeExit.cell[dimToUse]   - volumeEntry.cell[dimToUse])* dirInverse.cell[dimToUse];

#ifdef DEBUG
    if(debug>=15){
        printf("t_exit : %f\n",t_exit);
    }
#endif
    
    int cnt =0;
    while (t_entry <= t_exit ) {
        cnt++;
        if(cnt > MAXTRAVERSAL){
            printf("Exception : max stacktraversal exceeded for ray [%d,%d] in pulse %d\n",rayx,rayy,pulseIndex);
            break;
        }
        PEntry = vectAdd(volumeEntry, vectMult(ray.dir, t_entry));
        
#ifdef DEBUG
        if(debug>=15){
            printf("volumeEntry : %f,%f,%f\n",volumeEntry.x,volumeEntry.y,volumeEntry.z);
            printf("PEntry      : %f,%f,%f\n",PEntry.x,PEntry.y,PEntry.z);
            printf("ray.org     : %f,%f,%f\n",ray.org.x,ray.org.y,ray.org.z);
            printf("ray.dir     : %f,%f,%f\n",ray.dir.x,ray.dir.y,ray.dir.z);
            printf("t_entry     : %f\n",t_entry);
        }
#endif
        
        while (!KDT_ISLEAF(node)){  // Branch node
            hp.x = hp.y = hp.z = -666;

            if (PEntry.cell[KDT_DIMENSION(node)] < (node->branch.splitPosition-EPSILON)) {
#ifdef DEBUG
                if(debug>15)printf("[%f < %f (k:%d)], next node is %d\n",PEntry.cell[KDT_DIMENSION(node)],node->branch.splitPosition-EPSILON,KDT_DIMENSION(node),KDT_OFFSET(node));
#endif
                node = &(KdTree[KDT_OFFSET(node)]);
            }else if (PEntry.cell[KDT_DIMENSION(node)] > (node->branch.splitPosition+EPSILON)) {
#ifdef DEBUG
                if(debug>15)printf("[%f > %f (k:%d)], next node is %d\n",PEntry.cell[KDT_DIMENSION(node)],node->branch.splitPosition+EPSILON,KDT_DIMENSION(node),KDT_OFFSET(node)+1);
#endif
                node = &(KdTree[KDT_OFFSET(node)+1]);
            }else{
                // PEntry is on splitposition
                // Two different situations here:
                // 1. PEntry is on the split plane having arrived from neg side.
                // 2. PEntry is on the split plane having arrived from pos side.
                // in Case 1. next node is node
                // in Case 2. next node is node+1
                if (ray.org.cell[KDT_DIMENSION(node)] < node->branch.splitPosition-EPSILON){
#ifdef DEBUG
                    if(debug>15)printf("o[%f < %f (k:%d)] next node is %d\n",ray.org.cell[KDT_DIMENSION(node)],node->branch.splitPosition,KDT_DIMENSION(node),KDT_OFFSET(node));
#endif
                    node = &(KdTree[KDT_OFFSET(node)]);
                } else if (ray.org.cell[KDT_DIMENSION(node)] > node->branch.splitPosition+EPSILON){
#ifdef DEBUG
                    if(debug>15)printf("o[%f > %f (k:%d)] next node is %d\n",ray.org.cell[KDT_DIMENSION(node)],node->branch.splitPosition,KDT_DIMENSION(node),KDT_OFFSET(node)+1);
#endif
                    node = &(KdTree[KDT_OFFSET(node)+1]);
                } else {
                    // ray origin on split plane. Determine next node by direction of ray
                    //
                    if(ray.dir.cell[KDT_DIMENSION(node)] > 0 ){
#ifdef DEBUG
                        if(debug>15)printf("|[] next node is %d, dir is %f\n",KDT_OFFSET(node)+1,ray.dir.cell[KDT_DIMENSION(node)]);
#endif
                        node = &(KdTree[KDT_OFFSET(node)+1]);
                    }else {
                        // Includes situation where origin is on split position and ray is travelling parallel to split position
                        //
#ifdef DEBUG
                        if(debug>15)printf("|[] next node is %d dir is %f\n",KDT_OFFSET(node),ray.dir.cell[KDT_DIMENSION(node)]);
#endif
                        node = &(KdTree[KDT_OFFSET(node)]);
                    }
                }
            }
        }
        // Have a leaf now
        trisInLeaf = triangleListData[triangleListPtrs[KDT_OFFSET(node)]];
        hit.dist = 10e6;
#ifdef DEBUG
        if(debug>15)printf("Leaf found : %d triangles (hit.dist : %f)\n",trisInLeaf,hit.dist);
#endif
        for (i=0; i<trisInLeaf; i++){
            __global Triangle * tri = &(Triangles[triangleListData[triangleListPtrs[KDT_OFFSET(node)]+i+1]]);
            Intersect(tri, &ray, &hit, debug);
        }
        if((hit.trinum != NOINTERSECTION) && (hit.dist > 0.001)){
            // hitpoint may be outside box defining node and there may therefore be another
            // node beyond this node that has a nearer hitpoint. If so then follow rope to next
            // node
            hp = hitPoint(hit,ray);
            if(   hp.x <= (node->leaf.aabb.BB.x+EPSILON) && hp.x >= (node->leaf.aabb.AA.x-EPSILON)
               && hp.y <= (node->leaf.aabb.BB.y+EPSILON) && hp.y >= (node->leaf.aabb.AA.y-EPSILON)
               && hp.z <= (node->leaf.aabb.BB.z+EPSILON) && hp.z >= (node->leaf.aabb.AA.z-EPSILON))
            return hit;
        }
        
#ifdef DEBUG
        if(debug>15){
            printf("No Intersection found; Chasing rope...(hit.dist : %f, Hit: %2.8f,%2.8f,%2.8f)\n",hit.dist, hp.x,hp.y,hp.z);
            if( !(hp.x <= (node->leaf.aabb.BB.x+EPSILON) && hp.x >= (node->leaf.aabb.AA.x-EPSILON)) )printf("X hitpoint component not in AABB (%2.8f - %2.8f)\n",node->leaf.aabb.AA.x,node->leaf.aabb.BB.x);
            if( !(hp.y <= (node->leaf.aabb.BB.y+EPSILON) && hp.y >= (node->leaf.aabb.AA.y)-EPSILON) )printf("Y hitpoint component not in AABB (%2.8f - %2.8f)\n",node->leaf.aabb.AA.y,node->leaf.aabb.BB.y);
            if( !(hp.z <= (node->leaf.aabb.BB.z+EPSILON) && hp.z >= (node->leaf.aabb.AA.z-EPSILON)) )printf("Z hitpoint component not in AABB (%2.8f - %2.8f)\n",node->leaf.aabb.AA.z,node->leaf.aabb.BB.z);
        }
#endif

        // If ray doesnt intersect triangle in this leaf then propagate the
        // ray to an adjacent node using this leaf's Ropes. (rather than popping
        // back up the tree)
        // Set t_entry to be the ray distance to the adjacent AABB
        
        float t_max = 10e6;
        float tpos;
        int ropeInd, ropeIndSide, ropeIndOff;
        for(int i=0; i<3; i++){
            if (dirInverse.cell[i] != 0){
                t1 = (node->branch.aabb.AA.cell[i] - volumeEntry.cell[i]) * dirInverse.cell[i];
                t2 = (node->branch.aabb.BB.cell[i] - volumeEntry.cell[i]) * dirInverse.cell[i];
                if(t2-t1 > EPSILON && t2 >= 0){
#ifdef DEBUG
                    if(debug>15)printf("t2>t1 option\n");
#endif
                    tpos = t2;
                    ropeIndSide = 1;
                }else if( t1-t2 > EPSILON && t1 >= 0){
#ifdef DEBUG
                    if(debug>15)printf("t1>t2 option\n");
#endif
                    tpos = t1;
                    ropeIndSide = 0;
                }else{
#ifdef DEBUG
                    if(debug>15)printf("t1==t2 option\n");
#endif
                    // AABB is planar so select rope based upon direction of ray
                    tpos = t1;
                    ropeIndSide = (dirInverse.cell[i] < 0 ) ? 0 : 1;
                }
                if(tpos < t_max){
                    t_max = tpos;
                    ropeInd = i;
                    ropeIndOff = ropeIndSide;
                }
#ifdef DEBUG
                if(debug>=10)printf("Which rope? t1(%d)=%2.8f, t2(%d)=%2.8f, dirInverse(%d)=%f; tmax=%f ropeInd=%d, ropeIndOff=%d\n"
                                   ,i,t1,i,t2,i,dirInverse.cell[i],t_max,ropeInd,ropeIndOff);
#endif

            }
        }
        t_entry = t_max;
        
        // Follow the rope for this side to jump across to adjacent node
        
        if ( node->branch.Ropes[2*ropeInd+ropeIndOff] == NILROPE ) {
            hit.trinum = NOINTERSECTION;
            return hit ;
        }
#ifdef DEBUG
        if(debug>10)printf("Chased rope to node %d\n",node->branch.Ropes[2*ropeInd+ropeIndOff]);
#endif
        node = &(KdTree[node->branch.Ropes[2*ropeInd+ropeIndOff]]);
        
    }
    return hit ;
}

int clipToAABB(AABB boundingBox, VectorH *lineStart, VectorH *lineEnd){
    
    int status=0;
    ClipToBox(lineStart, lineEnd, boundingBox.AA, boundingBox.BB, &status);
    if (status == NOINTERSECTION) return 0;     // No intersect in any dim means no intersection with volume
    
    return 1;
    
}

// Cohen–Sutherland clipping algorithm clips a line from
// P0 = (x0, y0) to P1 = (x1, y1) against a rectangle with
// diagonal from (xmin, ymin) to (xmax, ymax).
void ClipToBox(VectorH *p0, VectorH *p1, VectorH min, VectorH max, int *status)
{
    // compute outcodes for P0, P1, and whatever point lies outside the clip rectangle
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
            double x, y, z, m;
            
            // At least one endpoint is outside the clip rectangle; pick it.
            OutCode outcodeOut = outcode0 ? outcode0 : outcode1;
            
            // Now find the intersection point;
            
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

OutCode ComputeOutCode(VectorH p, VectorH min, VectorH max)
{
    OutCode code;
    
    code = INSIDEAABB;          // initialised as being inside of clip window
    
    if (p.x < min.x)        // to the left of box
        code |= LEFTAABB;
    else if (p.x > max.x)   // to the right of box
        code |= RIGHTAABB;
    if (p.y < min.y)        // nearer that the box
        code |= NEARAABB;
    else if (p.y > max.y)   // farther than the box
        code |= FARAABB;
    if (p.z < min.z)        // Below bottom of box
        code |= BOTTOMAABB;
    else if(p.z > max.z)    // Above top of box
        code |= TOPAABB;
    
    return code;
}


void Intersect(__global Triangle *tri, Ray *ray, Hit *hit, int debug){
    unsigned int modulo[5];
    modulo[0] = 0; modulo[1] = 1; modulo[2] = 2; modulo[3] = 0; modulo[4]=1;
    
    int ku = modulo[tri->k+1];
    int kv = modulo[tri->k+2];
    
    const double nd = 1.0/(ray->dir.cell[tri->k] + ((tri->nd_u) * ray->dir.cell[ ku ]) + ((tri->nd_v) * ray->dir.cell[ kv ]) );
    const double thit = (tri->d - ray->org.cell[tri->k] - tri->nd_u * ray->org.cell[ ku ] - tri->nd_v * ray->org.cell[ kv ]) * nd;
#ifdef DEBUG
    if(debug>=20)printf("thit: %f hit->dist : %f\n",thit,hit->dist);
#endif
    // check for valid distance.
    if ( !(hit->dist > thit && thit >  EPSILON  ) ) return;
    
    // compute hitpoint positions on uv plane
    const double hu = (ray->org.cell[ku] + thit * ray->dir.cell[ ku ]);
    const double hv = (ray->org.cell[kv] + thit * ray->dir.cell[ kv ]);
    
    // check first barycentric coordinate
    const double beta = (hu * tri->kbu + hv * tri->kbv + tri->kbd);
    if (beta < 0.0f) return ;
    
    // check second barycentric coordinate￼
    const double gamma = (hu * tri->kcu + hv * tri->kcv + tri->kcd);
    if (gamma < 0.0f) return;
    
    // check third barycentric coordinate
    if (beta+gamma > 1.0f) return ;
    
    // have a valid hitpoint here. store it.
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
    VectorH N;
    N.cell[k] = 1;
    N.cell[ku] = T.nd_u;
    N.cell[kv] = T.nd_v;
    N = vectNorm(N);
    
    VectorH R;
    VectorH I = ray.dir;
    R = vectSub(I, vectMult(N, 2.0 * vectDot(I, N)) );
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

VectorH hitPoint (Hit h, Ray r){
    double x = r.org.x + h.dist * r.dir.x;
    double y = r.org.y + h.dist * r.dir.y;
    double z = r.org.z + h.dist * r.dir.z;
    return vectCreate(x, y, z);
}

int occluded(VectorH point, VectorH dir, AABB SceneBoundingBox,
             __global KdData * KdTree,
             __global int * triangleListData,
             __global int * triangleListPtrs,
             __global Triangle * Triangles){
    Hit h;
    Ray r;
    r.org = point;
    r.dir = dir;
    r.san = 0;
    h = StacklessTraverse(r,SceneBoundingBox,KdTree,triangleListData,triangleListPtrs,Triangles,0,0,0,0);
    if (h.trinum != NOINTERSECTION) return 1;
    return 0;
}

double reflectPower(double rayPower, VectorH L, VectorH R, VectorH V, Hit h, __global Triangle * Triangles, __global Texture * textureData, int debug){
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
    //    Texture texture = textureData[T.textureInd];
    int k = T.k;
    unsigned int modulo[5];
    modulo[0] = 0; modulo[1] = 1; modulo[2] = 2; modulo[3] = 0; modulo[4]=1;
    int ku = modulo[k+1];
    int kv = modulo[k+2];
    VectorH N;
    N.cell[k] = 1;
    N.cell[ku] = T.nd_u;
    N.cell[kv] = T.nd_v;
    N = vectNorm(N);

    double id = rayPower;
    double is = rayPower;
    double ia = 0;
    double ka = 0.0;
    double kd = 0.0 ;  //texture.kd;
    double ks = 1.0 ;  //texture.ks;
    double n  = 50.0 ; // texture.n;
    double LdotN = vectDot(L, N);
    double RdotV = vectDot(vectNorm(R), vectNorm(V));
    
    LdotN = (LdotN < 0) ? -LdotN : LdotN;   // removes 'one-way mirror' effect
    RdotV = (RdotV < 0) ? 0 : RdotV;        // if ray bouncing away from radar then there is no specular component
    double ans = (ia * ka ) +  (kd * LdotN * id) + ( ks * pow( RdotV , n) * is );

#ifdef DEBUG
    if(debug >= 25){
        VectorH vtmp,rtmp,ltmp;
        vtmp = V ; //vectNorm(V);
        rtmp = R ; //vectNorm(R);
        ltmp = L ; //vectNorm(L);
        printf("\n");
        printf("               Ray Scattering Properties: \n");
        printf(" ------------------------------------------------------\n");
        printf(" Incident specular power (is)               : %f\n",is);
        printf(" Incident diffuse power (id)                : %f\n",id);
        printf(" Ambient incident power (ia)                : %f\n",ia);
        printf(" Material ambient reflection constant (ka)  : %f\n",ka);
        printf(" Material diffuse reflection constant (kd)  : %f\n",kd);
        printf(" Material specular reflection constant (ks) : %f\n",ks);
        printf(" Material shininess constant (n)            : %f\n",n);
        printf(" Reflection vector (R)                      : %f, %f, %f\n",rtmp.x,rtmp.y,rtmp.z);
        printf(" Viewpoint vector (V)                       : %f, %f, %f\n",vtmp.x,vtmp.y,vtmp.z);
        printf(" Illumination vector (L)                    : %f, %f, %f\n",ltmp.x,ltmp.y,ltmp.z);
        printf(" L . N                                      : %f\n",LdotN);
        printf(" R . V                                      : %f\n",RdotV);
        printf(" Ambient power component                    : %f\n",ia*ka);
        printf(" Specular power component                   : %f \n", ( ks * pow( RdotV , n) * is));
        printf(" Diffuse power component                    : %f \n",(kd * LdotN * id));
        printf("Total power for this ray                    : %f\n",ans);
    }
#endif
    
    return ans;
}
