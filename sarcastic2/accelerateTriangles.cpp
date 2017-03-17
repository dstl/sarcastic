//
//  accelerateTriangles.cpp
//  sarcastic
//
//  Created by Darren Muff on 17/03/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#include "accelerateTriangles.hpp"
static const unsigned int quickmodulo[] = {0,1,2,0,1};

void accelerateTriangles(TriangleMesh *mesh, ATS **accelTriangles) {
    
    SPVector tmp,b,c;
    
    *accelTriangles = new ATS [mesh->triangles.size()] ;
    if (*accelTriangles == NULL) {
        printf("Error allocating memory for accelTriangles \n");
        exit(1);
    }
    
    for(int i=0; i<mesh->triangles.size(); ++i){
        SPVector Aa = mesh->vertices[mesh->triangles[i].a].asSPVector() ;
        SPVector Bb = mesh->vertices[mesh->triangles[i].b].asSPVector() ;
        SPVector Cc = mesh->vertices[mesh->triangles[i].c].asSPVector() ;
        SPVector N  = mesh->triangles[i].N.asSPVector() ;
        
        // calculate the projection dimension
        //
        float nx = N.x ;
        float ny = N.y ;
        float nz = N.z ;
        float anx = fabs(nx);
        float any = fabs(ny);
        float anz = fabs(nz);
        
        if( anx > any ){
            if (anx > anz){
                (*accelTriangles)[i].k=0;           /* X */
            }else{
                (*accelTriangles)[i].k=2;           /* Z */
            }
        }else{
            if ( any > anz){
                (*accelTriangles)[i].k=1;           /* Y */
            }else{
                (*accelTriangles)[i].k=2;           /* Z */
            }
        }
        
        int u = quickmodulo[(*accelTriangles)[i].k+1];
        int v = quickmodulo[(*accelTriangles)[i].k+2];
        
        (*accelTriangles)[i].nd_u = N.cell[u] / N.cell[(*accelTriangles)[i].k] ;
        (*accelTriangles)[i].nd_v = N.cell[v] / N.cell[(*accelTriangles)[i].k] ;
        VECT_SCMULT(N, 1/N.cell[(*accelTriangles)[i].k], tmp) ;
        (*accelTriangles)[i].d = VECT_DOT(Aa, tmp) ;
        
        VECT_SUB(Cc, Aa, b) ;
        VECT_SUB(Bb, Aa, c) ;
        
        float denom = 1.0/((b.cell[u]*c.cell[v]) - (b.cell[v]*c.cell[u]));
        
        (*accelTriangles)[i].kbu        = -b.cell[v] * denom ;
        (*accelTriangles)[i].kbv        =  b.cell[u] * denom ;
        (*accelTriangles)[i].kbd        = ((b.cell[v] * Aa.cell[u]) - (b.cell[u] * Aa.cell[v])) * denom;
        (*accelTriangles)[i].kcu        =  c.cell[v] * denom;
        (*accelTriangles)[i].kcv        = -c.cell[u] * denom;
        (*accelTriangles)[i].kcd        =  ((c.cell[u] * Aa.cell[v]) - (c.cell[v] * Aa.cell[u])) * denom;
        (*accelTriangles)[i].textureInd = mesh->triangles[i].mat ;
        (*accelTriangles)[i].triNum     = i;        
    }
    
    return ;
}
