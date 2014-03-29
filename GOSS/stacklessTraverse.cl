//***************************************************************************
//
//  stacklessTraverse.cl
//  GOSS
//
//  Created by Darren Muff on 02/08/2012.
//  Copyright (c) 2012 [dstl]. All rights reserved.
//
//
// Description:
//  OpenCl kernel code to perform a stackless traversal of a ray through
//  a KdTree. Inputs are the tree (and associated triangle primitives) and
//  an array of rays. Outputs are
//  ray intersections and ray directions
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

//#include "clipToAABB.cl"
//+++++++++++++++++++++++++++ Start of clipToAABB.cl +++++++++++++++++++++++++++++++++++++++++++++++

typedef int OutCode;

#define INSIDEAABB  0   // 000000
#define LEFTAABB    1   // 000001
#define RIGHTAABB   2   // 000010
#define NEARAABB    4   // 000100
#define FARAABB     8   // 001000
#define TOPAABB     16  // 010000
#define BOTTOMAABB  32  // 100000
#define NOINTERSECTION -1
#define TRUE 1
#define FALSE 0

int     clipToAABB    (AABB boundingBox, SPVector *lineStart, SPVector *lineEnd);
void    ClipToBox     (SPVector *p0,     SPVector *p1,        SPVector min, SPVector max, int *status);
OutCode ComputeOutCode(SPVector p,       SPVector min,        SPVector max);

int clipToAABB(AABB boundingBox, SPVector *lineStart, SPVector *lineEnd){
    
    int status=0;
    ClipToBox(lineStart, lineEnd, boundingBox.AA, boundingBox.BB, &status);
    if (status == NOINTERSECTION) return 0;     // No intersect in any dim means no intersection with volume
    
    return 1;
    
}

// Cohen–Sutherland clipping algorithm clips a line from
// P0 = (x0, y0) to P1 = (x1, y1) against a rectangle with
// diagonal from (xmin, ymin) to (xmax, ymax).
void ClipToBox(SPVector *p0, SPVector *p1, SPVector min, SPVector max, int *status)
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

OutCode ComputeOutCode(SPVector p, SPVector min, SPVector max)
{
    OutCode code;
    
    code = INSIDEAABB;      // initialised as being inside of clip window
    
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
//+++++++++++++++++++++++++++ End of clipToAABB.cl +++++++++++++++++++++++++++++++++++++++++++++++

//#include "intersect.cl"
//+++++++++++++++++++++++++++ Start of intersect.cl +++++++++++++++++++++++++++++++++++++++++++++++

#define EPSILON          ((double) 0.000001)

void Intersect(__global Triangle *tri, __global Ray *ray, __global Hit *hit);

void Intersect(__global Triangle *tri, __global Ray *ray, __global Hit *hit){
    unsigned int modulo[5];
    modulo[0] = 0; modulo[1] = 1; modulo[2] = 2; modulo[3] = 0; modulo[4]=1;
    
    int ku = modulo[tri->k+1];
    int kv = modulo[tri->k+2];
    
    const double nd = 1.0/(ray->dir.cell[tri->k] + ((tri->nd_u) * ray->dir.cell[ ku ]) + ((tri->nd_v) * ray->dir.cell[ kv ]) );
    const double thit = (tri->d - ray->org.cell[tri->k] - tri->nd_u * ray->org.cell[ ku ] - tri->nd_v * ray->org.cell[ kv ]) * nd;
    
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
//+++++++++++++++++++++++++++ End of intersect.cl +++++++++++++++++++++++++++++++++++++++++++++++

#define KDT_ISLEAF(n)    (n->leaf.flagDimAndOffset & (unsigned int)(1<<31))
#define KDT_DIMENSION(n) (n->leaf.flagDimAndOffset & 0x3)
#define KDT_OFFSET(n)    ((n->leaf.flagDimAndOffset & (0x7FFFFFFC))>>2)
#define NILROPE          ((int) -666666 )
#define NOINTERSECTION -1
#define MAXTRAVERSAL 1000           // Maximum kdTree traversal steps before we blow up

__kernel void stacklessTraverse(__global KdData * KdTree,
                                __global int * triangleListData,
                                __global int * triangleListPtrs,
                                __global Triangle * Triangles,
                                const    AABB SceneBoundingBox,
                                const    int nRays,                 // Number of rays
                                __global Ray * rays,                // array of rays to process.
                                __global Hit *hits                  // Location of ray hits
                                
                                
){
    int ind = get_global_id(0) ;
    int dimToUse ;
    int cnt, i ;
    int trisInLeaf ;
    float t1,t2,xinv,yinv,zinv;
    float t_entry, t_exit ;
    
    SPVector volumeEntry, volumeExit, PEntry, hp, dirInverse, v ;
    
    __global KdData * node ;
    
    if (ind >=0 && ind < nRays ) {
        t_entry = 0;
        t_exit  = VECT_MAG(rays[ind].org) + 1000 ;
        
        hits[ind].dist   = 10e6;
        hits[ind].trinum = NOINTERSECTION ;
        
        // Create an infinite line from the ray
        //
        volumeEntry  = rays[ind].org;
        VECT_SCMULT(rays[ind].dir, t_exit, v);
        VECT_ADD(rays[ind].org, v, volumeExit);
        
        // Calculate the ray segment within the scene volume
        // Do this early to reduce calcs for ray misses
        //
        if(!clipToAABB(SceneBoundingBox, &volumeEntry, &volumeExit)) return ;
        
        // Calc inverse direction so only multiplies (not divides) in loop (cheaper)
        //
        xinv = (rays[ind].dir.x == 0) ? 0 : 1./rays[ind].dir.x;
        yinv = (rays[ind].dir.y == 0) ? 0 : 1./rays[ind].dir.y;
        zinv = (rays[ind].dir.z == 0) ? 0 : 1./rays[ind].dir.z;
        VECT_CREATE(xinv, yinv, zinv, dirInverse);

        // initialise the root of the tree
        //
        node  = &(KdTree[0]);
        
        // get first non-zero dimension of direction vector
        //
        dimToUse = 0 ;
        while (dirInverse.cell[dimToUse] == 0)dimToUse++;
        
        // Set the extent of the ray 't' from the origin of the ray to the point where it exits
        // the volume
        //
        t_entry = 0.0 ;
        t_exit  = (volumeExit.cell[dimToUse]   - volumeEntry.cell[dimToUse]) * dirInverse.cell[dimToUse];
        
        cnt = 0 ;
        while ((t_entry <= t_exit) && (cnt++ <= MAXTRAVERSAL)) {
            VECT_SCMULT(rays[ind].dir, t_entry, v);
            VECT_ADD(volumeEntry, v, PEntry);
            
            while (!KDT_ISLEAF(node)){  // Branch node
                hp.x = hp.y = hp.z = -666;
                
                if (PEntry.cell[KDT_DIMENSION(node)] < (node->branch.splitPosition-EPSILON)) {
                    node = &(KdTree[KDT_OFFSET(node)]);
                }else if (PEntry.cell[KDT_DIMENSION(node)] > (node->branch.splitPosition+EPSILON)) {
                    node = &(KdTree[KDT_OFFSET(node)+1]);
                }else{
                    // PEntry is on splitposition
                    // Two different situations here:
                    // 1. PEntry is on the split plane having arrived from neg side.
                    // 2. PEntry is on the split plane having arrived from pos side.
                    // in Case 1. next node is node
                    // in Case 2. next node is node+1
                    if (rays[ind].org.cell[KDT_DIMENSION(node)] < node->branch.splitPosition-EPSILON){
                        node = &(KdTree[KDT_OFFSET(node)]);
                    } else if (rays[ind].org.cell[KDT_DIMENSION(node)] > node->branch.splitPosition+EPSILON){
                        node = &(KdTree[KDT_OFFSET(node)+1]);
                    } else {
                        // ray origin on split plane. Determine next node by direction of ray
                        //
                        if(rays[ind].dir.cell[KDT_DIMENSION(node)] > 0 ){
                            
                            node = &(KdTree[KDT_OFFSET(node)+1]);
                        }else {
                            // Includes situation where origin is on split position and ray is travelling parallel to split position
                            //
                            node = &(KdTree[KDT_OFFSET(node)]);
                        }
                    }
                }
            }
            
            // Have a leaf now
            //
            trisInLeaf = triangleListData[triangleListPtrs[KDT_OFFSET(node)]];
            hits[ind].dist = 10e6;
            
            for (i=0; i<trisInLeaf; i++){
                __global Triangle * tri = &(Triangles[triangleListData[triangleListPtrs[KDT_OFFSET(node)]+i+1]]);
                Intersect(tri, &(rays[ind]), &(hits[ind]));
            }
            
            if((hits[ind].trinum != NOINTERSECTION) && (hits[ind].dist > 0.001)){
                // hitpoint may be outside box defining node and there may therefore be another
                // node beyond this node that has a nearer hitpoint. If so then follow rope to next
                // node
                //
                VECT_SCMULT(rays[ind].dir, hits[ind].dist,hp);
                VECT_ADD(rays[ind].org, hp, hp);
                if(   hp.x <= (node->leaf.aabb.BB.x+EPSILON) && hp.x >= (node->leaf.aabb.AA.x-EPSILON)
                   && hp.y <= (node->leaf.aabb.BB.y+EPSILON) && hp.y >= (node->leaf.aabb.AA.y-EPSILON)
                   && hp.z <= (node->leaf.aabb.BB.z+EPSILON) && hp.z >= (node->leaf.aabb.AA.z-EPSILON))
                    return ;
            }
            
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
                    if(t2-t1 > EPSILON && t2 >= 0){
                        tpos = t2;
                        ropeIndSide = 1;
                    }else if( t1-t2 > EPSILON && t1 >= 0){
                        tpos = t1;
                        ropeIndSide = 0;
                    }else{
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
            
            t_entry = t_max ;
            
            // Follow the rope for this side to jump across to adjacent node
            //
            if ( node->branch.Ropes[2*ropeInd+ropeIndOff] == NILROPE ) {
                hits[ind].trinum = NOINTERSECTION;
                return ;
            }
            
            node = &(KdTree[node->branch.Ropes[2*ropeInd+ropeIndOff]]);
            
        }
        return ;
    }
    return ;
}

