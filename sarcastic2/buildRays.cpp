//
//  buildRays.cpp
//  sarcastic
//
//  Created by Darren Muff on 08/05/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#include "buildRays.hpp"


void buildRays(Ray **rayArray, int *nRays, int nAzRays, int nElRays, TriangleMesh *mesh, SPVector TxPos,
               double PowPerRay, AABB SceneBoundingBox,SPVector **rayAimPoints)
{
    
    int METHOD  = 1;
    
    // 1 - each ray aimed at triangle centre
    // 2 - random rays on each call across scene
    // 3 - random rays created first time but the same hitpoints used for each subsequent call
    // 4 - like 2 (random rays on each call across the scene) but rays are parallel from Tx
    
    if(METHOD == 1){
        *nRays = (int)mesh->triangles.size() ;
        *rayArray = (Ray *)sp_malloc(*nRays * sizeof(Ray));
        
        for(int i=0;i<mesh->triangles.size(); i++){
            SPVector mean;
            SPVector a = mesh->vertAforTri(i) ;
            SPVector b = mesh->vertBforTri(i) ;
            SPVector c = mesh->vertCforTri(i) ;
            
            for(int j=0; j<3; j++) mean.cell[j] = (a.cell[j] + b.cell[j] + c.cell[j]) / 3.0 ;
            
            Ray r;
            r.org = TxPos ;
            r.pow = PowPerRay ;
            r.len = 0;
            SPVector aimdir ;
            VECT_SUB(mean, TxPos, aimdir);
            VECT_NORM(aimdir, r.dir);
            SPVector zHat, Hdir, Vdir;
            VECT_CREATE(0, 0, 1, zHat);
            VECT_CROSS(r.dir, zHat, Hdir);
            VECT_CROSS(Hdir, r.dir, Vdir);
            VECT_NORM(Vdir, r.pol);
            (*rayArray)[i] = r;
        }
        return ;
        
    }else if(METHOD == 2){
        
        int nAzBeam = nAzRays;
        int nElBeam = nElRays;
        *nRays = nAzBeam*nElBeam ;
        
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
        
        *rayArray = (Ray *)sp_malloc(*nRays * sizeof(Ray));
        
        for( int i=0; i < *nRays; i++){
            SPVector elVect, azVect;
            double el, az;
            el = box_muller(minEl+((maxEl-minEl)/2), (maxEl-minEl)/2 );
            az = box_muller(minAz+((maxAz-minAz)/2), (maxAz-minAz)/2 );
            VECT_SCMULT(unitBeamEl, el, elVect);
            VECT_SCMULT(unitBeamAz, az, azVect);
            VECT_ADD(elVect, azVect, (*rayAimPoints)[i]) ;
        }
        
        SPVector Hdir,Vdir;
        
        for(int i=0; i< *nRays; i++){
            VECT_SUB((*rayAimPoints)[i], TxPos, (*rayArray)[i].dir );
            VECT_NORM((*rayArray)[i].dir, (*rayArray)[i].dir) ;
            (*rayArray)[i].org = TxPos ;
            (*rayArray)[i].pow = PowPerRay ;
            (*rayArray)[i].len = 0 ;
            VECT_CROSS((*rayArray)[i].dir, zHat, Hdir);
            VECT_CROSS(Hdir, (*rayArray)[i].dir, Vdir);
            VECT_NORM(Vdir, (*rayArray)[i].pol) ;
        }
        
        return;
        
    }else if(METHOD == 3){
        
        int nAzBeam = nAzRays;
        int nElBeam = nElRays;
        *nRays = nAzBeam*nElBeam ;
        
        if(*rayAimPoints == NULL){
            *rayAimPoints = (SPVector *)sp_malloc(sizeof(SPVector) * *nRays);
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
            
            for( int i=0; i < *nRays; i++){
                SPVector elVect, azVect;
                double el, az;
                el = box_muller(minEl+((maxEl-minEl)/2), (maxEl-minEl)/2 );
                az = box_muller(minAz+((maxAz-minAz)/2), (maxAz-minAz)/2 );
                VECT_SCMULT(unitBeamEl, el, elVect);
                VECT_SCMULT(unitBeamAz, az, azVect);
                VECT_ADD(elVect, azVect, (*rayAimPoints)[i]) ;
            }
        }
        
        *rayArray = (Ray *)sp_malloc(*nRays * sizeof(Ray));
        SPVector zHat,Hdir,Vdir;
        VECT_CREATE(0, 0, 1, zHat);
        
        for(int i=0; i< *nRays; i++){
            VECT_SUB((*rayAimPoints)[i], TxPos, (*rayArray)[i].dir );
            VECT_NORM((*rayArray)[i].dir, (*rayArray)[i].dir) ;
            (*rayArray)[i].org = TxPos ;
            (*rayArray)[i].pow = PowPerRay ;
            (*rayArray)[i].len = 0 ;
            VECT_CROSS((*rayArray)[i].dir, zHat, Hdir);
            VECT_CROSS(Hdir, (*rayArray)[i].dir, Vdir);
            VECT_NORM(Vdir, (*rayArray)[i].pol) ;
        }
        return ;
        
    }else if(METHOD == 4){      // 4 - like 2 (random rays on each call across the scene) but rays are parallel from Tx
        
        int nAzBeam = nAzRays;
        int nElBeam = nElRays;
        *nRays = nAzBeam*nElBeam ;
        
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
        
        *rayArray = (Ray *)sp_malloc(*nRays * sizeof(Ray));
        SPVector elVect, azVect, aimpnt,Opnt,Hdir,Vdir;
        
        for( int i=0; i < *nRays; i++){
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
        
        return;
    }
    
}
