//
//  rayTrace.cpp
//  sarcastic
//
//  Created by Darren Muff on 05/03/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#include "rayTrace.hpp"
#include <SIlib2/SIlib2.h>


void  rayTrace(TriangleMesh *mesh, kdTree::KdData **kdTree, int *numNodesInTree) {
    
    // generate rays
    
    // traverse tree for intersection and return coords
    
}

void stacklessTraverse(kdTree::KdData * KdTree,
                                 int * triangleListData,
                                 int * triangleListPtrs,
                                 ATS * accelTriangles,
                                 const    AABB SceneBoundingBox,
                                 const    int nRays,                 // Number of rays
                                 Ray * rays,                // array of rays to process.
                                 Hit *hits                  // Location of ray hits


){
//    int ind = get_global_id(0) ;
    int dimToUse ;
    int cnt, i ;
    int trisInLeaf ;
    float t1,t2,xinv,yinv,zinv;
    float t_entry, t_exit ;
    
    SPVector volumeEntry, volumeExit, PEntry, hp, dirInverse, v ;
    
     kdTree::KdData * node ;
    
    for(int ind =0; i<nRays; ++ind) {
//    if (ind >=0 && ind < nRays ) {
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
                
                if (PEntry.cell[KDT_DIMENSION(node)] < (node->brch.splitPosition-EPSILON)) {
                    node = &(KdTree[KDT_LEFTCHILD(node)]);
                }else if (PEntry.cell[KDT_DIMENSION(node)] > (node->brch.splitPosition+EPSILON)) {
                    node = &(KdTree[KDT_RGHTCHILD(node)]);
                }else{
                    // PEntry is on splitposition
                    // Two different situations here:
                    // 1. PEntry is on the split plane having arrived from neg side.
                    // 2. PEntry is on the split plane having arrived from pos side.
                    // in Case 1. next node is node
                    // in Case 2. next node is node+1
                    if (rays[ind].org.cell[KDT_DIMENSION(node)] < node->brch.splitPosition-EPSILON){
                        node = &(KdTree[KDT_LEFTCHILD(node)]);
                    } else if (rays[ind].org.cell[KDT_DIMENSION(node)] > node->brch.splitPosition+EPSILON){
                        node = &(KdTree[KDT_RGHTCHILD(node)]);
                    } else {
                        // ray origin on split plane. Determine next node by direction of ray
                        //
                        if(rays[ind].dir.cell[KDT_DIMENSION(node)] > 0 ){
                            node = &(KdTree[KDT_RGHTCHILD(node)]);
                        }else {
                            // Includes situation where origin is on split position and ray is travelling parallel to split position
                            //
                            node = &(KdTree[KDT_LEFTCHILD(node)]);
                        }
                    }
                }
            }
            
            // Have a leaf now
            //
            trisInLeaf = KDT_NUMTRIS(node) ;
            hits[ind].dist   = 10e6;
            hits[ind].trinum = NOINTERSECTION ;
            
            for (i=0; i<trisInLeaf; i++){
                
                ATS * tri = &(accelTriangles[KDT_INDEX(node)+1]) ;
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
            int ropeInd=0, ropeIndSide, ropeIndOff;
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

int clipToAABB(AABB boundingBox, SPVector *lineStart, SPVector *lineEnd){
    
    int status=0;
    ClipToBox(lineStart, lineEnd, boundingBox.AA, boundingBox.BB, &status);
    if (status == NOINTERSECTION) return 0;     // No intersect in any dim means no intersection with volume
    
    return 1;
    
}

void Intersect(ATS *tri, Ray *ray, Hit *hit){
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
