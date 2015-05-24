//
//  main.c
//  sarcastic
//
//  Created by Darren Muff on 03/10/2014.
//  Copyright (c) 2014 Dstl. All rights reserved.
//

#include <stdio.h>
#include <SIlib/SIlib.h>
#include "POTriangle.h"
#include "matrixMultiplication.h"
#define TXPOL "V"
#define RXPOL "-"
#define LAMBDA      ((double)0.3)       // Wavelegth in metres
#define RAYPOW      ((double)1.0e8)     // Transmit power used for incident E field = Pt * Gt
#define NPHIS       ((int)360)          // Number of observation points to calculate in azimuth
#define NTHETAS     ((int)90)           // Number of observation points to calculate in elevation
#define STARTPHI    ((double)0)         // Start azimuth angle (deg)
#define ENDPHI      ((double)360)       // End azimuth angle (deg)
#define STARTTHETA  ((double)0)         // Start angle of incidence (deg)
#define ENDTHETA    ((double)90)        // End angle of incidence (deg)
#define OBSDIST     ((double)20000.0)   // Observation distance in metres
#define OUTPUTDBS   1                   // Boolean - should output be in DBs?
#define ILLRANGE    ((double)20000.0)   // Range in metres from origin of illumination source
#define ILLAZ       ((double)180.0)     // Azimuth angke in degrees of source of illumination
#define ILLINC      ((double)30.0)      // Incidence ange in degrees of source of illumination
#define PRINTTRIS   0                   // Boolean - just print out triangles and quit
#define RCSOUTPUT   1
#define USEOPENCL   0
#define NRAYSX      10
#define NRAYSY      10
#define GOONLY      0                   // Perform calculations using geometrical optics rather than physical optics (PO)

void buildTriangle(SPVector AA, SPVector BB, SPVector CC, triangle * tri) ;
void readTriFile(triangle **tris, int *ntris, SPVector illuminationDir, const char *fname) ;
void oclPOField(triangle *tris, int ntris, Ray *rays, int nrays, SPVector *hitPoints, SPVector RxPnt, double k, SPVector RXVdir, SPVector RXHdir, SPCmplx *EsV, SPCmplx *EsH);
int PointInTriangle(SPVector p, SPVector p0, SPVector p1, SPVector p2);

int main(int argc, const char * argv[])
{
    
    // Read in triangles from triFile
    //
    triangle *tris ;
    int ntris ;
    SPVector zhat, Hpol, Vpol ;
    VECT_CREATE(0, 0, 1, zhat);
    
    // Define the illumination Origin in terms of range and azimuth, and incidence angle
    //
    SPVector illOrigin, illDir;
    double illRange, illAz, illInc, Ei_pow ;
    illRange = ILLRANGE ;
    illAz    = DEG2RAD(ILLAZ) ;
    illInc   = DEG2RAD(ILLINC) ;
    VECT_CREATE(illRange*sin(illInc)*cos(illAz), illRange*sin(illInc)*sin(illAz), illRange*cos(illInc), illOrigin) ;
    VECT_NORM(illOrigin, illDir) ;
    
//    readTriFile(&tris, &ntris, illDir ,"/Users/Darren/Development/DATA/triangles.tri") ;
    
    /*
     For debugging - allows you to set your own triangles without using a triangle file
     */
//    ntris = 1 ;
//    tris = sp_malloc(sizeof(triangle) * ntris);
//    SPVector AA,BB,CC ;
//    VECT_CREATE(   0, -1, 0.0, AA);
//    VECT_CREATE(   1, 1, 0.0, BB);
//    VECT_CREATE(   -1, 1, 0.0, CC);
//    double dx,dy;
//    dx=  0.0;
//    dy=  0.33333333333333331;
//    VECT_CREATE(   0+dx,  1+dy, 0.0, AA);
//    VECT_CREATE(  -1+dx, -1+dy, 0.0, BB);
//    VECT_CREATE(   1+dx, -1+dy, 0.0, CC);
//    buildTriangle(AA, BB, CC, &(tris[0])) ;
    
    ntris = 2 ;
    tris = sp_malloc(sizeof(triangle) * ntris);
    SPVector AA,BB,CC ;
    VECT_CREATE(  0.5,  0.5, 0.0, AA);
    VECT_CREATE( -0.5,  0.5, 0.0, BB);
    VECT_CREATE( -0.5, -0.5, 0.0, CC);
    buildTriangle(AA, BB, CC, &(tris[0])) ;
    
    VECT_CREATE(  -0.5, -0.5, 0.0, AA);
    VECT_CREATE(   0.5, -0.5, 0.0, BB);
    VECT_CREATE(   0.5,  0.5, 0.0, CC);
    buildTriangle(AA, BB, CC, &(tris[1])) ;
    
    
//    ntris = 4 ;
//    tris = sp_malloc(sizeof(triangle) * ntris);
//    SPVector AA,BB,CC ;
//    VECT_CREATE(  0.5, 0.5, 0.0, AA);
//    VECT_CREATE( -0.5, 0.5, 0.0, BB);
//    VECT_CREATE(  0.0, 0.0, 0.0, CC);
//    buildTriangle(AA, BB, CC, &(tris[0])) ;
//    
//    VECT_CREATE( -0.5,  0.5, 0.0, AA);
//    VECT_CREATE( -0.5, -0.5, 0.0, BB);
//    VECT_CREATE(  0.0,  0.0, 0.0, CC);
//    buildTriangle(AA, BB, CC, &(tris[1])) ;
//    
//    VECT_CREATE( -0.5, -0.5, 0.0, AA);
//    VECT_CREATE(  0.5, -0.5, 0.0, BB);
//    VECT_CREATE(  0.0,  0.0, 0.0, CC);
//    buildTriangle(AA, BB, CC, &(tris[2])) ;
//    
//    VECT_CREATE(  0.5, -0.5, 0.0, AA);
//    VECT_CREATE(  0.5,  0.5, 0.0, BB);
//    VECT_CREATE(  0.0,  0.0, 0.0,  CC);
//    buildTriangle(AA, BB, CC, &(tris[3])) ;
    
//    ntris = 8 ;
//    tris = sp_malloc(sizeof(triangle) * ntris);
//    SPVector AA,BB,CC ;
//    VECT_CREATE(  -1,  -1, 0.0, AA);
//    VECT_CREATE( 0,  -1, 0.0, BB);
//    VECT_CREATE(   0 ,0, 0.0, CC);
//    buildTriangle(AA, BB, CC, &(tris[0])) ;
//    
//    VECT_CREATE(   0, -1, 0.0, AA);
//    VECT_CREATE(  1, -1, 0.0, BB);
//    VECT_CREATE(   0 ,0, 0.0, CC);
//    buildTriangle(AA, BB, CC, &(tris[1])) ;
//    
//    VECT_CREATE(   -1,  0, 0.0, AA);
//    VECT_CREATE(   -1, -1, 0.0, BB);
//    VECT_CREATE(   0 ,0, 0.0, CC);
//    buildTriangle(AA, BB, CC, &(tris[2])) ;
//
//    VECT_CREATE(   1,  -1, 0.0, AA);
//    VECT_CREATE(   1, 0, 0.0, BB);
//    VECT_CREATE(   0 ,0, 0.0, CC);
//    buildTriangle(AA, BB, CC, &(tris[3])) ;
//    
//    VECT_CREATE(   -1,  1, 0.0, AA);
//    VECT_CREATE(   -1, 0, 0.0, BB);
//    VECT_CREATE(   0 ,0, 0.0, CC);
//    buildTriangle(AA, BB, CC, &(tris[4])) ;
//    
//    VECT_CREATE(   1, 0, 0.0, AA);
//    VECT_CREATE(   1, 1, 0.0, BB);
//    VECT_CREATE(   0 ,0, 0.0, CC);
//    buildTriangle(AA, BB, CC, &(tris[5])) ;
//    
//    VECT_CREATE(   0, 1, 0.0, AA);
//    VECT_CREATE(  -1, 1, 0.0, BB);
//    VECT_CREATE(   0 ,0, 0.0, CC);
//    buildTriangle(AA, BB, CC, &(tris[6])) ;
//    
//    VECT_CREATE(  1,  1, 0.0, AA);
//    VECT_CREATE( 0,  1, 0.0, BB);
//    VECT_CREATE( 0, 0, 0.0, CC);
//    buildTriangle(AA, BB, CC, &(tris[7])) ;
    
    if(PRINTTRIS){
        exit(0);
    }
    

    // We need parallel illumination so that we don't have to worry about transmit power or
    // phase effects from curved wavefronts
    //
    
    int nRays = NRAYSX * NRAYSY;
    Ray *rays       =      (Ray *)sp_malloc(sizeof(Ray) * nRays) ;
    SPVector * hits = (SPVector *)sp_malloc(sizeof(SPVector) * nRays);
    int * triForHit =      (int *)sp_malloc(sizeof(int) * nRays);
    int * hitsForTri=      (int *)sp_malloc(sizeof(int) * ntris);
    
    double xmin,xmax,ymin,ymax,x,y,delx,dely;
    int ix, iy, index;
    xmin = -0.5;
    xmax = +0.5;
    ymin = -0.5;
    ymax = +0.5;
    delx = (xmax-xmin) / NRAYSX ;
    printf("Hits\n");
    for(iy=0; iy<NRAYSY; iy++){
        for(ix=0; ix<NRAYSX; ix++){
            x = xmin + ix*delx ;
            y = ymin + iy*dely ;
            index = iy*NRAYSX+ix;
            VECT_CREATE( x, y, 0, hits[index]);
            printf("%f,%f\n",hits[index].x,hits[index].y);
            SPVector dir;
            VECT_SUB(hits[index], illOrigin, dir);
            VECT_NORM(dir, rays[index].dir);
            rays[index].org = illOrigin ;
            rays[index].pow = RAYPOW ;
            rays[index].len = 0.0;
            if(fabs(VECT_DOT(rays[index].dir, zhat)) >= 0.99 ){
                VECT_CREATE(1, 0, 0, Hpol);
            }else{
                VECT_CROSS(rays[index].dir, zhat, Hpol) ;
                VECT_NORM(Hpol, Hpol) ;
            }
            VECT_CROSS(Hpol, rays[index].dir, Vpol) ;
            if(TXPOL == "V"){
                rays[index].pol = Vpol ;
            }else{
                rays[index].pol = Hpol ;
            }
        }
    }
    
    // Now trace each ray through the volume to find the hits
    // The purpose here is to build a set of hit points 
    //
    
    
    if (RCSOUTPUT) {
        Ei_pow = RAYPOW / (4 * SIPC_pi * ILLRANGE * ILLRANGE) ; // power at facet in watts/m^2
    }
//    
//    for (int i=0; i<nRays; i++) {
//        rays[i].org = illOrigin ;
//        rays[i].pow = RAYPOW ;
//        VECT_SUB(hits[i], rays[i].org, rays[i].dir);
//        rays[i].len = VECT_MAG(rays[i].dir);
//        VECT_NORM(rays[i].dir, rays[i].dir);
//        rays[i].pow = RAYPOW / (4.0 * SIPC_pi*rays[i].len*rays[i].len);
//        if(fabs(VECT_DOT(rays[i].dir, zhat)) >= 0.99 ){
//            VECT_CREATE(1, 0, 0, Hpol);
//        }else{
//            VECT_CROSS(rays[i].dir, zhat, Hpol) ;
//            VECT_NORM(Hpol, Hpol) ;
//        }
//        VECT_CROSS(Hpol, rays[i].dir, Vpol) ;
//        if(TXPOL == "V"){
//            rays[i].pol = Vpol ;
//        }else{
//            rays[i].pol = Hpol ;
//        }
//        
//        // find triangle for this ray
//        //
//        for(int t=0; t<ntris; t++){
//            if( PointInTriangle(hits[i], tris[t].AA, tris[t].BB, tris[t].CC) ){
//                triForHit[i] = t;
//                hitsForTri[t]++ ;
////                printf("Ray %d hits triangle %d at (%f,%f)\n",i,triForHit[i],hits[i].x,hits[i].y);
//                t=ntris;
//            }else{
//                triForHit[i] = -1;
//            }
//
//        }
//    }
    
        
    /*SPVector vhatpar,vhatper, parvec, pervec ;
    double par,per ;

    if(fabs(VECT_DOT(illDir, zhat)) == 1.0 ){
        VECT_CREATE(1, 0, 0, vhatper);
    }else{
        VECT_CROSS(illDir, zhat, vhatper) ;
        VECT_NORM(vhatper, vhatper) ;
    }
    VECT_CROSS(illDir, vhatper, vhatpar) ;
    VECT_NORM(vhatpar, vhatpar) ;

    
    for (int i=0; i<nRays; i++) {
        par = VECT_DOT(tris[i].MP, vhatpar) ;
        per = VECT_DOT(tris[i].MP, vhatper) ;
        VECT_SCMULT(vhatpar, par, parvec) ;
        VECT_SCMULT(vhatper, per, pervec) ;
        VECT_ADD(illOrigin, parvec, rays[i].org) ;
        VECT_ADD(rays[i].org, pervec, rays[i].org) ;
        VECT_CREATE(-illDir.x, -illDir.y, -illDir.z, rays[i].dir) ;
        rays[i].pow = RAYPOW ;
        rays[i].pow = rays[i].pow / (4.0 * SIPC_pi*illRange*illRange);
        rays[i].len = 0.0 ;
        if(fabs(VECT_DOT(rays[i].dir, zhat)) >= 0.99 ){
            VECT_CREATE(1, 0, 0, Hpol);
        }else{
            VECT_CROSS(rays[i].dir, zhat, Hpol) ;
            VECT_NORM(Hpol, Hpol) ;
        }
        VECT_CROSS(Hpol, rays[i].dir, Vpol) ;
        if(TXPOL == "V"){
            rays[i].pol = Vpol ;
        }else{
            rays[i].pol = Hpol ;
        }
    }
    
    */
    
    int iphi, niphis;
    int itheta, nitheta;
    double startPhi,endPhi,startTheta,endTheta;
    niphis     = NPHIS ;
    nitheta    = NTHETAS ;
    startPhi   = DEG2RAD(STARTPHI) ;
    endPhi     = DEG2RAD(ENDPHI) ;
    startTheta = DEG2RAD(STARTTHETA) ;
    endTheta   = DEG2RAD(ENDTHETA) ;
    
    double deltaiphi, deltaitheta;
    deltaiphi = (endPhi-startPhi) / niphis ;
    deltaitheta = (endTheta-startTheta) / (nitheta) ;
    double phi_s, theta_s;
    SPVector RxPnt, obsDir ;
    double obsDist = OBSDIST ;
    double Erx;
    double k = 2 * SIPC_pi / LAMBDA ;

    for(itheta=0; itheta<nitheta; itheta++){
        theta_s = startTheta + itheta * deltaitheta ;
        
        for(iphi=0; iphi < niphis; iphi++){
            phi_s = startPhi + iphi * deltaiphi ;
            
            obsDir.x = sin(theta_s) * cos(phi_s);
            obsDir.y = sin(theta_s) * sin(phi_s);
            obsDir.z = cos(theta_s);
            VECT_SCMULT(obsDir, obsDist, RxPnt) ;
            
            // For debugging - provides a handy way to see the expected phase at the receiver
            //
//            double psum,p,p1;
//            psum=0.0;
//             for (int n=0; n<nRays; n++ ){
//                 SPVector tmp1,tmp2;
//                 VECT_SUB(hits[n], rays[n].org, tmp1);
//                 p = (-(VECT_MAG(tmp1))*2*SIPC_pi/LAMBDA) ;
//                 p1 = (p / (2*SIPC_pi)) - ((int)(p / (2*SIPC_pi)));
//                 p1 = p1 * 2 * SIPC_pi ;
//                 printf("GO Phase TX->Hitpoint : %f\n",RAD2DEG(p1));
//                 
//                 VECT_SUB(RxPnt, hits[n], tmp2);
//                 p = (-VECT_MAG(tmp2)*2*SIPC_pi/LAMBDA) ;
//                 p1 = (p / (2*SIPC_pi)) - ((int)(p / (2*SIPC_pi)));
//                 p1 = p1 * 2 * SIPC_pi ;
//                 printf("GO Phase Hitpoint->RX : %f\n",RAD2DEG(p1));
//                 
//                 p = (-(VECT_MAG(tmp1)+VECT_MAG(tmp2))*2*SIPC_pi/LAMBDA)-SIPC_pi ;
//                 p1 = (p / (2*SIPC_pi)) - ((int)(p / (2*SIPC_pi)));
//                 p1 = p1 * 2 * SIPC_pi ;
//                 printf("GO Phase at RX for ray %d should be %f radians (%f deg)\n",n ,p , RAD2DEG(p1));
//                 psum+=RAD2DEG(p1);
//            }
//            printf("Mean GO Phase at RX is %f deg (+PI = %f)\n",psum/nRays,(psum/nRays)+180);

            
            // Define unit vectors for V & H fields
            // We do this as the definition of H and V from the sensor may not be truly horizontal
            // or vertical from the sensor. By allowing us to define the V & H directions we can
            // accurately model the V & H at the sensor.
            // We also do it here as we dont want POTriangle() calculating this each time it is called
            //
            SPVector RXVdir, RXHdir ;
            VECT_CREATE(0, 0, 1, zhat) ;
            if(VECT_DOT(zhat, obsDir) >= 0.9999999){
                // Ie looking from above - for diagnostic purposes - remove from SAR
                // simulation as this never occurs.
                //
                VECT_CREATE(rays[0].pol.x, rays[0].pol.y, 0, RXVdir) ;
                VECT_CROSS(RXVdir, zhat, RXHdir) ;
            }else{
                VECT_CROSS(zhat, obsDir, RXHdir) ;
                VECT_CROSS(obsDir, RXHdir, RXVdir) ;
            }
            VECT_NORM(RXVdir, RXVdir) ;
            VECT_NORM(RXHdir, RXHdir) ;
            
            SPCmplx EsV1, EsH1, EsV, EsH ;

            EsV.r = EsV.i = EsH.r = EsH.i = 0 ;
            
            if(USEOPENCL){
                oclPOField(tris, ntris, rays, ntris, hits, RxPnt, k, RXVdir, RXHdir, &EsV, &EsH);
            }else if(GOONLY){
                double r,p,phs;
                SPVector v;
                for (int n=0; n<nRays; n++ ){
                    r=0;
                    if(triForHit[n] >= 0){
                        VECT_SUB(RxPnt, hits[n], v);
                        r += VECT_MAG(v);
                        VECT_SUB(hits[n], illOrigin, v);
                        r += VECT_MAG(v);
                        p = sqrt(RAYPOW) / (4 * SIPC_pi * r * r) ;
                        phs = (-2.0 * SIPC_pi * r / LAMBDA) - SIPC_pi ;
                        EsV1.r = p * cos(phs);
                        EsV1.i = p * sin(phs);
                        CMPLX_ADD(EsV, EsV1, EsV);
                        EsH.r = EsH.i = 0.0;
                        // For debugging - provides a way to see individual scattering contributions
                        //
//                        printf("EsV1[%d] : %f deg amp : %e v/m  (%f, %f)\n",n,RAD2DEG(CMPLX_PHASE(EsV1)), CMPLX_MAG(EsV1),EsV1.r,EsV1.i);
                    }
                }
            }else{
                for (int n=0; n<nRays; n++ ){
                    if(triForHit[n] >= 0){
                        POField(tris[triForHit[n]], rays[n], hits[n], RxPnt, k,  RXVdir, RXHdir, &EsV1, &EsH1);
                        CMPLX_SCMULT(1.0/hitsForTri[triForHit[n]], EsV1, EsV1);
                        CMPLX_SCMULT(1.0/hitsForTri[triForHit[n]], EsH1, EsH1);
                        CMPLX_ADD(EsV, EsV1, EsV) ;
                        CMPLX_ADD(EsH, EsH1, EsH) ;
                        // For debugging - provides a way to see individual scattering contributions
                        //
//                        printf("EsV1[%d] : %f deg amp : %e v/m  (%f, %f)\n",n,RAD2DEG(CMPLX_PHASE(EsV1)), CMPLX_MAG(EsV1),EsV1.r,EsV1.i);
                    }
                }
            }
            
            if (RXPOL == "V") {
                double EsV_mag = CMPLX_MAG(EsV);
                if (RCSOUTPUT) {
                    EsV_mag =  4*EsV_mag*EsV_mag * 4.0 * SIPC_pi * OBSDIST * OBSDIST / Ei_pow ;
                }
                if (OUTPUTDBS) {
                    Erx = 10*log10(EsV_mag) ;
                }else{
                    Erx = EsV_mag ;
                }
                if(Erx <  0.0){
                    printf("%f, %f, %f \n",0.0,phi_s,theta_s);
                }else{
                    printf("%f, %f, %f \n",Erx,phi_s,theta_s);
                }
            }else if (RXPOL == "H"){
                double EsH_mag = CMPLX_MAG(EsH) ;
                if (RCSOUTPUT) {
                    EsH_mag =  4*EsH_mag*EsH_mag * 4.0 * SIPC_pi * OBSDIST * OBSDIST / Ei_pow ;
                }
                if (OUTPUTDBS) {
                    Erx = 10*log10(EsH_mag) ;
                }else{
                    Erx = EsH_mag ;
                }
                if(Erx < 0.0){
                    printf("%f, %f, %f \n",0.0,phi_s,theta_s);
                }else{
                    printf("%f, %f, %f \n",Erx,phi_s,theta_s);
                }
            }else{
                SPVector Ev, Eh, E;
                double a, E_mag, p ;
                a = CMPLX_MAG(EsV);
                p = CMPLX_PHASE(EsV) ;
                VECT_SCMULT(RXVdir, a, Ev);
                a = CMPLX_MAG(EsH);
                VECT_SCMULT(RXHdir, a, Eh);
                VECT_ADD(Ev, Eh, E) ;
                E_mag = VECT_MAG(E) ;
                if(RCSOUTPUT){
                    E_mag =  4*E_mag*E_mag * 4.0 * SIPC_pi * OBSDIST * OBSDIST / Ei_pow ;
                }
                if (OUTPUTDBS) {
                    Erx = 10*log10(E_mag) ;
                }else{
                    Erx = E_mag ;
                }
                // For debugging - if you want to see the scattering phase uncomment the next line
                //
                //Erx = RAD2DEG(CMPLX_PHASE(EsV));
                
                if(Erx < 00.0){
                    printf("%f, %f, %f \n",0.0,phi_s,theta_s);
                }else{
                    printf("%f, %f, %f \n",Erx,phi_s,theta_s);
//                    double r; SPVector v,v1;
//                    r=0;
//                    VECT_CREATE(0, 1, 0, v);
//                    VECT_SUB(RxPnt, v, v1);
//                    r += VECT_MAG(v1);
//                    VECT_SUB(v, illOrigin, v1);
//                    r += VECT_MAG(v1);
//                
//                    p = (-(r)*2*SIPC_pi/LAMBDA)-SIPC_pi ;
//                    double p1 = (p / (2*SIPC_pi));
//                    p1 = p1 - ((int)(p / (2*SIPC_pi)));
//                    p1 = p1 * 2 * SIPC_pi ;
//                    printf("GO Phase at RX should be %f radians (%f deg)\n",p , RAD2DEG(p1));
                }
            }
            
        }
    }

    free(tris) ;
    free(rays) ;
    free(hits) ;
    return 0;
}

void readTriFile(triangle **tris, int *ntris, SPVector illuminationDir, const char *fname){
    
    /*
     Output File Structure
     int number_of_triangles
     int sizeof_vertex_information
     int length_of_material_names
     foreach triangle{
     int triangle_number
     <sizeof_vertex_information> vertexA.x
     <sizeof_vertex_information> vertexA.y
     <sizeof_vertex_information> vertexA.z
     <sizeof_vertex_information> vertexB.x
     <sizeof_vertex_information> vertexB.y
     <sizeof_vertex_information> vertexB.z
     <sizeof_vertex_information> vertexC.x
     <sizeof_vertex_information> vertexC.y
     <sizeof_vertex_information> vertexC.z
     <sizeof_vertex_information> normal.x
     <sizeof_vertex_information> normal.y
     <sizeof_vertex_information> normal.z
     <sizeof_vertex_information> midpoint.x
     <sizeof_vertex_information> midpoint.y
     <sizeof_vertex_information> midpoint.z
     <sizeof_vertex_information>  triangleArea
     <sizeof_vertex_information>  globalTolocalMatrix[0]
     <sizeof_vertex_information>  globalTolocalMatrix[1]
     <sizeof_vertex_information>  globalTolocalMatrix[2]
     <sizeof_vertex_information>  globalTolocalMatrix[3]
     <sizeof_vertex_information>  globalTolocalMatrix[4]
     <sizeof_vertex_information>  globalTolocalMatrix[5]
     <sizeof_vertex_information>  globalTolocalMatrix[6]
     <sizeof_vertex_information>  globalTolocalMatrix[7]
     <sizeof_vertex_information>  globalTolocalMatrix[8]
     <sizeof_vertex_information>  localToglobalMatrix[0]
     <sizeof_vertex_information>  localToglobalMatrix[1]
     <sizeof_vertex_information>  localToglobalMatrix[2]
     <sizeof_vertex_information>  localToglobalMatrix[3]
     <sizeof_vertex_information>  localToglobalMatrix[4]
     <sizeof_vertex_information>  localToglobalMatrix[5]
     <sizeof_vertex_information>  localToglobalMatrix[6]
     <sizeof_vertex_information>  localToglobalMatrix[7]
     <sizeof_vertex_information>  localToglobalMatrix[8]
     <sizeof_vertex_information>  localToglobalMatrix[9]
     char[<MATBYTES>] Material_name
     }
     */
    
    int nt;
    int vertexBytes;
    int materialBytes;
    
    FILE *fpin ;
    fpin = fopen(fname, "r");
    if (fpin == NULL) {
        printf("Error : failed to open input triangle file %s\n",fname);
        exit(1);
    }
    fread(&nt, sizeof(int), 1, fpin);
    *ntris = nt ;
    fread(&vertexBytes,sizeof(int), 1, fpin);
    fread(&materialBytes, sizeof(int), 1, fpin);
    
    if (vertexBytes != sizeof(double) ) {
        printf("Error : Triangle vertices of size %d bytes not yet supported\n",vertexBytes);
        exit(1);
    }
    if (materialBytes != MATBYTES ) {
        printf("Error : Material Bytes in file does not match the program\n");
        exit(1);
    }
    
    triangle *tris_tmp = (triangle *)sp_malloc(sizeof(triangle) * nt);
    triangle tri ;
    int trind = 0;
    
    for (int itri=0; itri < nt; itri++ ) {
        fread(&(tri.id), sizeof(int), 1, fpin) ;
        fread(&(tri.AA.x), sizeof(double), 1, fpin) ;
        fread(&(tri.AA.y), sizeof(double), 1, fpin) ;
        fread(&(tri.AA.z), sizeof(double), 1, fpin) ;
        fread(&(tri.BB.x), sizeof(double), 1, fpin) ;
        fread(&(tri.BB.y), sizeof(double), 1, fpin) ;
        fread(&(tri.BB.z), sizeof(double), 1, fpin) ;
        fread(&(tri.CC.x), sizeof(double), 1, fpin) ;
        fread(&(tri.CC.y), sizeof(double), 1, fpin) ;
        fread(&(tri.CC.z), sizeof(double), 1, fpin) ;
        fread(&(tri.NN.x), sizeof(double), 1, fpin) ;
        fread(&(tri.NN.y), sizeof(double), 1, fpin) ;
        fread(&(tri.NN.z), sizeof(double), 1, fpin) ;
        fread(&(tri.area), sizeof(double), 1, fpin) ;
        for (int i=0; i<9; i++){
            fread(&(tri.globalToLocalMat[i]), sizeof(double), 1, fpin) ;
        }
        for (int i=0; i<9; i++){
            fread(&(tri.localToGlobalMat[i]), sizeof(double), 1, fpin) ;
        }
        fread(&(tri.matId) , sizeof(int), 1, fpin) ;
        
        if (VECT_DOT(tri.NN, illuminationDir) > 0.000001) {
            tris_tmp[trind].id = tri.id ;
            tris_tmp[trind].AA.x = tri.AA.x ;
            tris_tmp[trind].AA.y = tri.AA.y ;
            tris_tmp[trind].AA.z = tri.AA.z ;
            tris_tmp[trind].BB.x = tri.BB.x ;
            tris_tmp[trind].BB.y = tri.BB.y ;
            tris_tmp[trind].BB.z = tri.BB.z ;
            tris_tmp[trind].CC.x = tri.CC.x ;
            tris_tmp[trind].CC.y = tri.CC.y ;
            tris_tmp[trind].CC.z = tri.CC.z ;
            tris_tmp[trind].NN.x = tri.NN.x ;
            tris_tmp[trind].NN.y = tri.NN.y ;
            tris_tmp[trind].NN.z = tri.NN.z ;
            tris_tmp[trind].area = tri.area ;
            for (int i=0; i<9; i++){
                tris_tmp[trind].globalToLocalMat[i] = tri.globalToLocalMat[i] ;
            }
            for (int i=0; i<9; i++){
                tris_tmp[trind].localToGlobalMat[i] = tri.localToGlobalMat[i] ;
            }
//            printf("Triangle %d : \n",tri.id) ;
//            printf("    A      : %3.6f,%3.6f,%3.6f\n",tri.AA.x,tri.AA.y,tri.AA.z );
//            printf("    B      : %3.6f,%3.6f,%3.6f\n",tri.BB.x,tri.BB.y,tri.BB.z );
//            printf("    C      : %3.6f,%3.6f,%3.6f\n",tri.CC.x,tri.CC.y,tri.CC.z );
//            printf("    N      : %3.6f,%3.6f,%3.6f\n",tri.NN.x,tri.NN.y,tri.NN.z );
//            printf("    M      : %3.6f,%3.6f,%3.6f\n",tri.MP.x,tri.MP.y,tri.MP.z );
            if(PRINTTRIS){
                printf("%3.6f,%3.6f,%3.6f\n",tri.AA.x,tri.AA.y,tri.AA.z );
                printf("%3.6f,%3.6f,%3.6f\n",tri.BB.x,tri.BB.y,tri.BB.z );
                printf("%3.6f,%3.6f,%3.6f\n",tri.CC.x,tri.CC.y,tri.CC.z );
            }

            trind++;
        }else{
            *ntris = *ntris-1 ;
        }
        
    }
    *tris = tris_tmp ;
    
    fclose(fpin);
}

void buildTriangle(SPVector AA, SPVector BB, SPVector CC, triangle * tri){
    

    SPVector NN, l1,l2, zhat ;
    double alpha, beta;
    double T_dash[9];
    double T_dashdash[9];
    
    tri->id = 0;
    tri->AA = AA ;
    tri->BB = BB ;
    tri->CC = CC ;
    VECT_SUB(BB, AA, l1) ;
    VECT_SUB(CC, AA, l2) ;
    VECT_CROSS(l1, l2, NN) ;
    tri->area = VECT_MAG(NN) * 0.5 ;
    VECT_NORM(NN, tri->NN) ;
    tri->matId = 0;
    
    VECT_CREATE(0, 0, 1, zhat);
    alpha = atan2(tri->NN.y, tri->NN.x);
    beta  = acos(VECT_DOT(zhat, tri->NN));
    T_dash[0] = cos(alpha);
    T_dash[1] = sin(alpha);
    T_dash[2] = 0;
    T_dash[3] = -sin(alpha);
    T_dash[4] = cos(alpha);
    T_dash[5] = 0;
    T_dash[6] = 0;
    T_dash[7] = 0;
    T_dash[8] = 1;
    T_dashdash[0] = cos(beta);
    T_dashdash[1] = 0;
    T_dashdash[2] = -sin(beta);
    T_dashdash[3] = 0;
    T_dashdash[4] = 1;
    T_dashdash[5] = 0;
    T_dashdash[6] = sin(beta);
    T_dashdash[7] = 0;
    T_dashdash[8] = cos(beta);
    
    matmul(T_dashdash, T_dash, tri->globalToLocalMat, 3, 3, 3, 3);
    mat3by3inv(tri->globalToLocalMat, tri->localToGlobalMat);

    return ;
}


int PointInTriangle(SPVector p, SPVector p0, SPVector p1, SPVector p2)
{
    double s = (p0.y * p2.x) - (p0.x * p2.y) + ((p2.y - p0.y) * p.x) + ((p0.x - p2.x) * p.y);
    double t = (p0.x * p1.y) - (p0.y * p1.x) + ((p0.y - p1.y) * p.x) + ((p1.x - p0.x) * p.y);
    
    if ((s < 0) != (t < 0))
        return 0;
    
    double A = (-p1.y * p2.x) + (p0.y * (p2.x - p1.x)) + (p0.x * (p1.y - p2.y)) + (p1.x * p2.y);
    if (A < 0.0){
        s = -s;
        t = -t;
        A = -A;
    }
    if (s > 0 && t > 0 && (s + t) < A){
        return 1 ;
    }else{
        return 0 ;
    }
}
