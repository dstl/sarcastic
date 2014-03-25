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

#include "SPVector.cl"
#include "structures.cl"
#include "clipToAABB.cl"
#include "intersect.cl"

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

