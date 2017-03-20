//
//  rayTrace.cpp
//  sarcastic
//
//  Created by Darren Muff on 05/03/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#include "rayTrace.hpp"
#include <SIlib2/SIlib2.h>


void buildRays(Ray **rayArray, int nAzRays, int nElRays, SPVector TxPos, AABB SceneBoundingBox) ;

void  rayTrace(TriangleMesh *mesh, kdTree::KdData *kdTree, int *numNodesInTree) {
    
    kdTree::KdData *node;
    
    node = &(kdTree[0]) ;
    
    // build Ropes and Boxes
    //
    int Ropes[6] ;
    AABB sceneAABB = kdTree[0].brch.aabb ;

    for(int i=0; i<6; i++) Ropes[i] = NILROPE;
    BuildRopesAndBoxes(node, Ropes, sceneAABB, kdTree);
    
    // generate rays
    //
    Ray *rays;
    int nAzRays = 10;
    int nElRays = 10;
    int nRays = nAzRays * nElRays ;
    SPVector TxPos;
    VECT_CREATE(-100.0, 0.0, 100.0, TxPos) ;

    buildRays(&rays, nAzRays, nElRays, TxPos, sceneAABB) ;
    
    // generate triangles
    //
    ATS *accelTriangles;
    accelerateTriangles(mesh,&accelTriangles) ;
    Hit * hits = new Hit [nRays] ;
    
    Timer rayTimer ;
    SPStatus status ;
    im_init_status(status, 0);

    startTimer(&rayTimer, &status) ;
    
    shootRay(kdTree, accelTriangles, nRays, rays, hits);
    
   
    endTimer(&rayTimer, &status);
    
    printf("%d rays traced in  %f seconds (%f rays / second) \n",nRays,timeElapsedInSeconds(&rayTimer, &status), nRays / timeElapsedInSeconds(&rayTimer, &status)) ;
    
    SPVector hp;
    for (int i=0; i<nRays; ++i) {
        if(hits[i].trinum != -1){
            VECT_SCMULT(rays[i].dir, hits[i].dist, hp) ;
            VECT_ADD(rays[i].org, hp, hp) ;
            printf("%f, %f, %f\n",hp.x,hp.y,hp.z);
        }
    }
    delete [] accelTriangles ;
    return ;
    
}

void shootRay(kdTree::KdData * KdTree,ATS * accelTriangles,const int nRays, Ray * rays, Hit *hits){
    for (int ind=0; ind<nRays; ++ind) {
        stacklessTraverse(ind, KdTree, accelTriangles, nRays, rays, hits);
    }
}


void stacklessTraverse(const int ind,           // Index of ray to trace in rays array
                       kdTree::KdData * KdTree, // The Kd-Tree - prebuilt and passed in here
                       ATS * accelTriangles,    // An array of triangles in accelerated format
                       const int nRays,         // Number of rays - size of rays array and hits array
                       Ray * rays,              // array of rays to process.
                       Hit *hits                // Location of ray hits
){
    int dimToUse ;
    int cnt, i ;
    int trisInLeaf ;
    float t1,t2,xinv,yinv,zinv;
    float t_entry, t_exit ;
    
    SPVector volumeEntry, volumeExit, PEntry, hp, dirInverse, v ;
    AABB sceneBoundingBox ;
    
    kdTree::KdData * node ;
    
    sceneBoundingBox = KdTree[0].brch.aabb ;
    
    if (ind >=0 && ind < nRays ) {
        if (ind == 6) {
            printf("here index is 6\n");
        }
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
        if(!clipToAABB(sceneBoundingBox, &volumeEntry, &volumeExit)) return ;
        
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
                
                int triInx = KdTree[KDT_INDEX(node)+1+i].leaf.triangleIndex ;
                ATS * tri = &(accelTriangles[triInx]) ;
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
                   && hp.z <= (node->leaf.aabb.BB.z+EPSILON) && hp.z >= (node->leaf.aabb.AA.z-EPSILON)){
//                    printf("%f,%f,%f\n",hp.x,hp.y,hp.z);
                    if (hp.z == 0.0) {
                        printf("Here\n");
                    }
                    return ;
                }
            }
            
            // If ray doesnt intersect triangle in this leaf then propagate the
            // ray to an adjacent node using this leaf's Ropes. (rather than popping
            // back up the tree)
            // Set t_entry to be the ray distance to the adjacent AABB
            //
            float t_max = 10e6;
            float tpos;
            int ropeInd=0, ropeIndSide, ropeIndOff=0;
            for(int i=0; i<3; i++){
                if (dirInverse.cell[i] != 0){
                    t1 = (node->brch.aabb.AA.cell[i] - volumeEntry.cell[i]) * dirInverse.cell[i];
                    t2 = (node->brch.aabb.BB.cell[i] - volumeEntry.cell[i]) * dirInverse.cell[i];
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
            if ( node->brch.Ropes[2*ropeInd+ropeIndOff] == NILROPE ) {
                hits[ind].trinum = NOINTERSECTION;
                return ;
            }
            
            node = &(KdTree[node->brch.Ropes[2*ropeInd+ropeIndOff]]);
            
        }
//        return ;
    }
    return ;
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

//void BuildRopesAndBoxes(kdTree::KdData * Node, int *RS, AABB aabb, kdTree::KdData * KdTree){
//    
//    int Sl, Sr; //SplitPosLeft, SPlitPosRight;
//    float V;
//    int RSLeft[6], RSRight[6]; // Ropes for left and right side
//    AABB aabbLeft, aabbRight;
//    kdTree::KdData * Nr, * Nl;
//    
//    if( KDT_ISLEAF(Node) ){
//        for(int i=0; i<3; i++){
//            Node->leaf.aabb.AA.cell[i] = aabb.AA.cell[i];
//            Node->leaf.aabb.BB.cell[i] = aabb.BB.cell[i];
//            Node->leaf.Ropes[2*i] = RS[2*i];
//            Node->leaf.Ropes[2*i+1] = RS[2*i+1];
//        }
//    } else {
//        for(int i=0; i<3; i++){
//            Node->brch.aabb.AA.cell[i] = aabb.AA.cell[i];
//            Node->brch.aabb.BB.cell[i] = aabb.BB.cell[i];
//        }
//        
//        Sl = KDT_DIMENSION(Node) * 2;
//        Sr = KDT_DIMENSION(Node) * 2 + 1;
//        V  = Node->brch.splitPosition;
//        for(int i=0; i<6; i++)RSLeft[i]=RSRight[i]=RS[i];
//        
//        RSLeft[Sr] = KDT_RGHTCHILD(Node); // Right child for Node
//        aabbLeft = aabb;
//        aabbLeft.BB.cell[KDT_DIMENSION(Node)] = V;
//        Nl = &(KdTree[KDT_LEFTCHILD(Node)]);
//        
//        BuildRopesAndBoxes(Nl, RSLeft, aabbLeft, KdTree);
//        
//        RSRight[Sl] = KDT_LEFTCHILD(Node); // Left child for Node
//        aabbRight = aabb;
//        aabbRight.AA.cell[KDT_DIMENSION(Node)] = V;
//        Nr = &(KdTree[KDT_RGHTCHILD(Node)]);
//        
//        BuildRopesAndBoxes(Nr, RSRight, aabbRight, KdTree);
//        
//    }
//    
//    return ;
//}

void buildRays(Ray **rayArray, int nAzRays, int nElRays, SPVector TxPos, AABB SceneBoundingBox){
    int nAzBeam = nAzRays;
    int nElBeam = nElRays;
    
    double PowPerRay = 1.0;
    int nRays = nAzBeam*nElBeam ;
    
    SPVector rVect,zHat,unitBeamAz,unitBeamEl;
    VECT_MINUS( TxPos, rVect ) ;
    VECT_CREATE(0, 0, 1., zHat) ;
    VECT_CROSS(rVect, zHat, unitBeamAz);
    VECT_NORM(unitBeamAz, unitBeamAz) ;
    VECT_CROSS(unitBeamAz, rVect, unitBeamEl) ;
    VECT_NORM(unitBeamEl, unitBeamEl) ;
    
    double maxEl,maxAz, minEl, minAz;
    maxEl = maxAz = minEl = minAz = 0.0 ;
    
    SPVector boxPts[8];
    VECT_CREATE(SceneBoundingBox.AA.x, SceneBoundingBox.AA.y, SceneBoundingBox.AA.z, boxPts[0]);
    VECT_CREATE(SceneBoundingBox.AA.x, SceneBoundingBox.BB.y, SceneBoundingBox.AA.z, boxPts[1]);
    VECT_CREATE(SceneBoundingBox.BB.x, SceneBoundingBox.BB.y, SceneBoundingBox.AA.z, boxPts[2]);
    VECT_CREATE(SceneBoundingBox.BB.x, SceneBoundingBox.AA.y, SceneBoundingBox.AA.z, boxPts[3]);
    VECT_CREATE(SceneBoundingBox.AA.x, SceneBoundingBox.AA.y, SceneBoundingBox.BB.z, boxPts[4]);
    VECT_CREATE(SceneBoundingBox.AA.x, SceneBoundingBox.BB.y, SceneBoundingBox.BB.z, boxPts[5]);
    VECT_CREATE(SceneBoundingBox.BB.x, SceneBoundingBox.BB.y, SceneBoundingBox.BB.z, boxPts[6]);
    VECT_CREATE(SceneBoundingBox.BB.x, SceneBoundingBox.AA.y, SceneBoundingBox.BB.z, boxPts[7]);
    
    for( int k=0; k<8; k++){
        double El = VECT_DOT(boxPts[k], unitBeamEl) ;
        double Az = VECT_DOT(boxPts[k], unitBeamAz) ;
        maxEl = ( maxEl < El ) ? El : maxEl ;
        maxAz = ( maxAz < Az ) ? Az : maxAz ;
        minEl = ( minEl > El ) ? El : minEl ;
        minAz = ( minAz > Az ) ? Az : minAz ;
    }
    
    SPVector aimpoint;
    VECT_CREATE(SceneBoundingBox.AA.x+((SceneBoundingBox.BB.x - SceneBoundingBox.AA.x)/2),
                SceneBoundingBox.AA.y+((SceneBoundingBox.BB.y - SceneBoundingBox.AA.y)/2),
                SceneBoundingBox.AA.z+((SceneBoundingBox.BB.z - SceneBoundingBox.AA.z)/2), aimpoint);
    
    *rayArray = (Ray *)sp_malloc(nRays * sizeof(Ray));
    SPVector elVect, azVect, aimpnt,Opnt,Hdir,Vdir;
    
    for( int i=0; i < nRays; i++){
        double el, az;
        el = box_muller(minEl+((maxEl-minEl)/2), (maxEl-minEl)/2 );
        az = box_muller(minAz+((maxAz-minAz)/2), (maxAz-minAz)/2 );
        VECT_SCMULT(unitBeamEl, el, elVect);
        VECT_SCMULT(unitBeamAz, az, azVect);
        VECT_ADD(elVect, azVect, aimpnt);
        Opnt = TxPos ;
        VECT_ADD(Opnt, elVect, Opnt);
        VECT_ADD(Opnt, azVect, Opnt);
        VECT_SUB(aimpnt, Opnt, (*rayArray)[i].dir );
        VECT_NORM((*rayArray)[i].dir, (*rayArray)[i].dir) ;
        (*rayArray)[i].org = Opnt ;
        (*rayArray)[i].pow = PowPerRay ;
        (*rayArray)[i].len = 0 ;
        VECT_CROSS((*rayArray)[i].dir, zHat, Hdir);
        VECT_CROSS(Hdir, (*rayArray)[i].dir, Vdir);
        VECT_NORM(Vdir, (*rayArray)[i].pol) ;
    }

    return ;
}

//void accelerateTriangles(TriangleMesh *mesh, ATS **accelTriangles) {
//
//    SPVector tmp,b,c;
//    
//    *accelTriangles = new ATS [mesh->triangles.size()] ;
//    if (*accelTriangles == NULL) {
//        printf("Error allocating memory for accelTriangles \n");
//        exit(1);
//    }
//    
//    for(int i=0; i<mesh->triangles.size(); ++i){
//        SPVector Aa = mesh->vertAforTri(i) ;
//        SPVector Bb = mesh->vertBforTri(i) ;
//        SPVector Cc = mesh->vertCforTri(i) ;
//        SPVector N  = mesh->triangles[i].N.asSPVector() ;
//        
//        // calculate the projection dimension
//        //
//        float nx = N.x ;
//        float ny = N.y ;
//        float nz = N.z ;
//        float anx = fabs(nx);
//        float any = fabs(ny);
//        float anz = fabs(nz);
//        
//        if( anx > any ){
//            if (anx > anz){
//                (*accelTriangles)[i].k=0;           /* X */
//            }else{
//                (*accelTriangles)[i].k=2;           /* Z */
//            }
//        }else{
//            if ( any > anz){
//                (*accelTriangles)[i].k=1;           /* Y */
//            }else{
//                (*accelTriangles)[i].k=2;           /* Z */
//            }
//        }
//    
//        int u = quickmodulo[(*accelTriangles)[i].k+1];
//        int v = quickmodulo[(*accelTriangles)[i].k+2];
//    
//        (*accelTriangles)[i].nd_u = N.cell[u] / N.cell[(*accelTriangles)[i].k] ;
//        (*accelTriangles)[i].nd_v = N.cell[v] / N.cell[(*accelTriangles)[i].k] ;
//        VECT_SCMULT(N, 1/N.cell[(*accelTriangles)[i].k], tmp) ;
//        (*accelTriangles)[i].d = VECT_DOT(Aa, tmp) ;
//        
//        VECT_SUB(Cc, Aa, b) ;
//        VECT_SUB(Bb, Aa, c) ;
//    
//        float denom = 1.0/((b.cell[u]*c.cell[v]) - (b.cell[v]*c.cell[u]));
//        
//        (*accelTriangles)[i].kbu        = -b.cell[v] * denom ;
//        (*accelTriangles)[i].kbv        =  b.cell[u] * denom ;
//        (*accelTriangles)[i].kbd        = ((b.cell[v] * Aa.cell[u]) - (b.cell[u] * Aa.cell[v])) * denom;
//        (*accelTriangles)[i].kcu        =  c.cell[v] * denom;
//        (*accelTriangles)[i].kcv        = -c.cell[u] * denom;
//        (*accelTriangles)[i].kcd        =  ((c.cell[u] * Aa.cell[v]) - (c.cell[v] * Aa.cell[u])) * denom;
//        (*accelTriangles)[i].textureInd = mesh->triangles[i].mat ;
//        (*accelTriangles)[i].triNum     = i;
//        
//        mesh->triangles[i].buildTransMats() ;
//        
//    }
//
//    return ;
//}
//
