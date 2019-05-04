/***************************************************************************
 * 
 *           Module :  threadCore.cpp
 *          Program :  bircs
 *       Created by :  Darren Muff on 12/06/2017
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *   Description:
 *     Core thread routine for BIRCS
 * 
 *   (c) Crown Copyright 2018 Defence Science and Technology Laboratory
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software")
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 * 
 ***************************************************************************/
#include <fftw3.h>
extern "C" {
#include "boxMullerRandom.h"
#include "RCS.h"
#include "printProgress.h"

}
#include "ranf.h"
#include "threadCore.hpp"
#include <sys/time.h>
#include "buildKernel.hpp"
#include "rayTrace.hpp"
#include "cpuPOField.hpp"
#include <sarclib/sarclib.h>

#define NOINTERSECTION -1


void * devPulseBlock ( void * threadArg ) {
    
    struct bircsThreadData *td;
    td = (struct bircsThreadData *) threadArg;
    
    int nAzBeam, nElBeam, bounceToShow, nrnpItems, nShadowRays, tid ;
    int nbounce, nxRay, nyRay, nRays, reflectCount, iray, nShadows, interrogate ;
    int maxRaysPerBounce, pol ;
    int reportN = 10 ;
    int rayGenMethod ;

    int nThreads, nPulses, startPulse;
    int hrs,min,sec;
    bool monostatic = true;
    bool calcShadowRays = true ;
    
    double gainRx, PowPerRay, derampRange ;
    double *ranges, *newranges, interogRad, k ;
    double pCentDone, sexToGo ;
    
    struct tm *p;
    time_t current,complete;
    char ct[1000];

    AABB SceneBoundingBox ;
    rangeAndPower * rnp ;
    Ray *rayArray, *newRays, *reflectedRays, *shadowRays, *LRays, *RRays;
    Hit *hitArray, *newHits, *shadowHits;
    
    SPVector aimdir_tx, aimdir_rx, RxPos, TxPos, origin, interogPt;
    SPVector maxRCS = {-9e99,0.0,0.0};
    SPStatus status;
    
    FILE *interogFP ;
    
    Timer threadTimer ;
    ATS *accelTriangles = NULL;
    AABB sceneAABB ;
    
    kdTree::KdData * tree = NULL;
    int treeSize;
    int nTriangles;
    
    tid              = td->tid ;
    nThreads         = td->nThreads ;
    startPulse       = td->startPulse ;
    nPulses          = td->nPulses ;
    gainRx           = td->gainRx ;
    PowPerRay        = td->PowPerRay ;
    bounceToShow     = td->bounceToShow ;
    nAzBeam          = td->nAzBeam ;
    nElBeam          = td->nElBeam ;
    SceneBoundingBox = td->SceneBoundingBox ;
    interrogate      = td->interrogate ;
    interogPt        = td->interogPt ;
    interogRad       = td->interogRad ;
    interogFP        = *(td->interogFP) ;
    k                = 2 * SIPC_pi / (SIPC_c / td->freq_centre) ;
    pol              = td->polarisation ;
    rayGenMethod     = td->rayGenMethod ;

    
    VECT_CREATE(0, 0, 0, origin);
    im_init_status(status, 0) ;
    
    // **** loop  start here
    //
    
    startTimer(&threadTimer, &status) ;
    SPVector *rayAimPoints = NULL ;
    
    tree = *(td->tree) ;
    treeSize = td->treesize ;
    accelTriangles = *(td->accelTriangles) ;
    TriangleMesh newMesh ;
    newMesh = *(td->sceneMesh) ;
    nTriangles   = (int)newMesh.triangles.size() ; // number of triangles in array 'triangles'
    
    for (int pulse=0; pulse<nPulses; pulse++){
        int pulseIndex = pulse + startPulse ;
        if(tid == 0){
            if( pulse % reportN == 0 && nPulses != 1){
                pCentDone = 100.0*pulse/nPulses ;
                printProgress(pCentDone, 60);
                
                if(pulse != 0 ){
                    current  = time(NULL);
                    sexToGo  = estimatedTimeRemaining(&threadTimer, pCentDone, &status);
                    hrs      = (int)floor(sexToGo/60/60) ;
                    min      = (int)floor(sexToGo / 60) - (hrs*60);
                    sec      = (int)sexToGo % 60;
                    complete = current + sexToGo ;
                    p        = localtime(&complete) ;
                    strftime(ct, 1000, "%d/%b/%y %H:%M", p);
                    printf(" ETC %s (in %2.0dh:%2.0dm:%2.0ds) ",ct,hrs,min,sec);
                }else{
                    printf(" Calculating ETC...");
                }
            }
        }
        // Set correct parameters for beam to ray trace
        //
        TxPos = td->TxPositions[pulseIndex] ;
        RxPos = td->RxPositions[pulseIndex] ;
        if( (TxPos.x == RxPos.x) && (TxPos.y == RxPos.y) && (TxPos.z == RxPos.z) )  monostatic = true;
        double phi = atan2(RxPos.y,RxPos.x);
        double theta = atan2(sqrt(RxPos.x*RxPos.x+RxPos.y*RxPos.y),RxPos.z);
        
        
        // Generate a distribution of nAzbeam x nElbeam rays that originate from the TxPosition aiming at the origin. Use beamMax as the std deviation
        // for the distribution
        //
        nbounce = 0;
        nxRay   = nAzBeam ;
        nyRay   = nElBeam ;
        nRays   = nxRay*nyRay;
        buildRays(&rayArray, &nRays, nAzBeam, nElBeam, &newMesh, TxPos, PowPerRay, td->SceneBoundingBox, &rayAimPoints, rayGenMethod, pol);
        maxRaysPerBounce = nRays;  // Use this for memory as nxRay/nyRay may be incorrect if buildRays set to triangle centres
        
        // Set up deramp range for this pulse
        //
        VECT_MINUS(TxPos, aimdir_tx) ;
        VECT_MINUS(RxPos, aimdir_rx) ;
        derampRange = (VECT_MAG(aimdir_tx) + VECT_MAG(aimdir_rx)) / 2.0;
        
        // Use Calloc for rnp as we will be testing for zeroes later on
        //
        rnp = (rangeAndPower *)sp_calloc(maxRaysPerBounce*MAXBOUNCES, sizeof(rangeAndPower));
        
        while ( nbounce < MAXBOUNCES &&  nRays != 0){
            
            // If this pulse is monostatic then the first bounce ray intersections by definition will be
            // visible to the receiver and so no shadowrays are needed to be calculated
            // This is quite a saving
            //
            if(nbounce == 0 && monostatic){
                calcShadowRays = false ;
            }else{
                calcShadowRays = true ;
            }
            
            // Malloc space for hits for this bounce
            //
            hitArray = (Hit *)sp_malloc(nRays * sizeof(Hit));
            
            // Cast the rays in rayArray through the KdTree using a stackless traversal technique
            // return the hit locations in hitArray
            //
            shootRay(tree, accelTriangles, nRays, rayArray, hitArray) ;
            
            // Sort out hits
            //
            reflectCount = 0 ;
            iray         = 0 ;
            int *edgeHit = (int *)sp_malloc(sizeof(int)*nRays);
            
            // How many hits occured on this ray cast
            //
            for(int i=0; i<nRays; i++) if ( hitArray[i].trinum != NOINTERSECTION ) {
                // Check to see if out hit / reflection point was near to a vertex or edge by looking at the hitpoints
                //
                if(hitArray[i].u < 0.001 || hitArray[i].v < 0.001 || hitArray[i].u+hitArray[i].v > 0.999){
                    edgeHit[i] = 1;
                }else{
                    edgeHit[i] = 0;
                }
                if(edgeHit[i] == 0)
                    reflectCount++ ;
            }
            
            if( reflectCount == 0) break ;
            
            // shrink rayArray and hitArray to get rid of misses
            //
            newRays       = (Ray *)sp_malloc(reflectCount * sizeof(Ray)) ;
            newHits       = (Hit *)sp_malloc(reflectCount * sizeof(Hit)) ;
            reflectedRays = (Ray *)sp_malloc(reflectCount * sizeof(Ray)) ;
            
            for (int i=0; i<nRays; i++) {
                if ( hitArray[i].trinum != NOINTERSECTION && edgeHit[i]==0 ){
                    newRays[iray] = rayArray[i] ;
                    newHits[iray] = hitArray[i] ;
                    iray++ ;
                }
            }
            free(rayArray) ;
            free(hitArray) ;
            free(edgeHit) ;
            rayArray = newRays ;
            hitArray = newHits ;
            nRays    = reflectCount ;
            
            // Build forward scattering rays ready for next turn round the loop
            //
            reflect(nRays, rayArray, hitArray, accelTriangles, reflectedRays);

            
            // If debug out is required then capturing here using the origins of the Reflected rays
            // which will save us having to calculate the hit locations again
            //
            if(bounceToShow == nbounce){
                printf("Scene Intersection points for rays on bounce %d\n",nbounce);
                printf("------------------------------------------------\n");
                for (int i=0; i<nRays; i++){
                    printf("%f,%f,%f\n",reflectedRays[i].org.x,reflectedRays[i].org.y,reflectedRays[i].org.z);
                }
            }
            
            // Malloc space for shadowRays (arrays from a hit going back to receiver)
            //
            shadowRays =    (Ray *)sp_malloc(nRays * sizeof(Ray)) ;
            shadowHits =    (Hit *)sp_malloc(nRays * sizeof(Hit)) ;
            ranges     = (double *)sp_malloc(nRays * sizeof(double)) ;
            
            // Build Shadowrays
            //
            buildShadowRays(nRays, RxPos, reflectedRays, shadowRays, ranges) ;
            
            // Work out which rays have a path back to receiver using stackless traverse kernel
            //
            if(calcShadowRays)shootRay(tree, accelTriangles, nRays, shadowRays, shadowHits) ;
            
            // Shrink the shadowRays to only include those that made it back to the sensor
            // in order to calculate power at sensor we also need the Illumination or LRays
            //
            nShadows = 0 ;
            double *facing = (double *)sp_malloc(sizeof(double) * nRays);
            int *hitsOnEachTri  = (int *)sp_calloc(nTriangles, sizeof(int));
            
            SPVector normal;
            for(int i=0; i<nRays; i++){
                // A ray might be on the plane of a triangle having arrived from the side obscured from the receiver
                // (infinitely thin triangle plane). To make sure rays do not pass through a triangle, make sure the normal
                // is pointing in the same direction as the reflected ray
                //
                normal = newMesh.triangles[hitArray[i].trinum].N.asSPVector() ;
                if( VECT_DOT(normal, reflectedRays[i].dir) < 0){
                    VECT_MINUS(normal, normal);
                }
                // shadowRays propagate from RxPoint so when facing[i] is negative then hitpoint is visible to receiver
                //
                facing[i] =VECT_DOT(normal, shadowRays[i].dir);
                // Keep a count of the number of times a triangle is hit as we calculate the facet RCS per intersection
                //
                hitsOnEachTri[ hitArray[i].trinum ]++ ;
                // If the shadow ray intersects a facet further away that the hit point then cull it
                //
                if ( shadowHits[i].dist > ranges[i]) shadowHits[i].trinum = NOINTERSECTION ;
                // Now count the number of shadow rays so that we can resize the arrays
                //
                if ( ((shadowHits[i].trinum == NOINTERSECTION) && (facing[i] < -0.1) && (hitsOnEachTri[hitArray[i].trinum] <= 1)) || !calcShadowRays ){
                    nShadows++ ;
                }
            }
            free(hitsOnEachTri);
            
            if( nShadows != 0){
                
                newRays       = (Ray *)sp_malloc(nShadows * sizeof(Ray)) ;
                newHits       = (Hit *)sp_malloc(nShadows * sizeof(Hit)) ;
                LRays         = (Ray *)sp_malloc(nShadows * sizeof(Ray)) ;
                RRays         = (Ray *)sp_malloc(nShadows * sizeof(Ray)) ;
                newranges     = (double *)sp_malloc(nShadows * sizeof(double)) ;
                hitsOnEachTri = (int *)sp_calloc(nTriangles, sizeof(int));
                
                iray     = 0 ;
                for (int i=0; i<nRays; i++) {
                    hitsOnEachTri[ hitArray[i].trinum ]++ ;
                    // Remove any shadowRays that are from a triangle whose normal is not in the same direction as the Receiver
                    //
                    if ( ((shadowHits[i].trinum == NOINTERSECTION) && (facing[i] < -0.1) && (hitsOnEachTri[hitArray[i].trinum] <= 1)) || !calcShadowRays ){
                        newRays[iray] = shadowRays[i] ;
                        newHits[iray] = hitArray[i] ;
                        LRays[iray]   = rayArray[i] ;
                        RRays[iray]   = reflectedRays[i] ;
                        newranges[iray] = ranges[i] ;
                        iray++ ;
                    }
                }
                free(hitsOnEachTri);
                free(shadowRays) ;
                free(hitArray) ;
                free(ranges);
                shadowRays  = newRays  ;
                hitArray    = newHits  ;
                nShadowRays = nShadows ;
                ranges      = newranges ;
                
                // At this point the size of all arrays is nShadowRays not nRays
                //
                
                // For each ray that isn't occluded back to the receiver, calculate the power and put it into rnp.
                //
                cpuPOField(&newMesh, hitArray, nShadowRays, LRays, shadowRays, RxPos, k, ranges, gainRx, nbounce+1, pol, &(rnp[maxRaysPerBounce*nbounce])) ;
                
                free(LRays);
                free(RRays);
            }
            
            free(facing);
            free(shadowRays);
            free(shadowHits);
            free(ranges);
            
            memcpy(rayArray, reflectedRays, sizeof(Ray)*reflectCount);
            free(reflectedRays);
            
            
            // Reset rayArrays and nRays for next time round loop
            //
            nRays = reflectCount ;
            
            nbounce++;
            
            free(hitArray);
        }
        
        free(rayArray) ;
        
        int cnt=0;
        nrnpItems = 0;
        nRays = maxRaysPerBounce ;
        for (int i=0; i<maxRaysPerBounce*MAXBOUNCES; i++){
            if ((rnp[i].Es.r * rnp[i].Es.r + rnp[i].Es.i * rnp[i].Es.i) != 0 && rnp[i].range !=0) cnt++ ;
            //            printf("[%d] range: %f, %e, %f \n",i,rnp[i].range,rnp[i].Es.r,rnp[i].Es.i);
        }
        if(cnt > 0){
            nrnpItems = cnt;
        }
        
        SPCmplxD targtot = {0.0,0.0};
        double cmplx_mag = 0.0;
        
        if(nrnpItems > 0){
            rnpData_t * rnpData = (rnpData_t *)sp_malloc(nrnpItems * sizeof(rnpData_t));
            cnt = 0;
            for (int i=0; i<maxRaysPerBounce*MAXBOUNCES; i++){
                if (CMPLX_MAG(rnp[i].Es) != 0 && rnp[i].range !=0){
                    rnpData[cnt].Es    = rnp[i].Es ;
                    rnpData[cnt].rdiff = rnp[i].range - derampRange;
                    cnt++ ;
                }
            }
            
            for (int i=0; i<nrnpItems; i++){
                CMPLX_ADD(rnpData[i].Es, targtot, targtot);
            }
            cmplx_mag = 10*log10(RCS(PowPerRay, CMPLX_MAG(targtot), VECT_MAG(aimdir_tx), VECT_MAG(aimdir_rx)));
            free(rnpData);
        }
        
        if(cmplx_mag <0 ){
            td->results[pulseIndex].r     = 0.0 ;
            td->results[pulseIndex].theta = theta ;
            td->results[pulseIndex].phi   = phi ;
        }else{
            td->results[pulseIndex].r     = cmplx_mag ;
            td->results[pulseIndex].theta = theta ;
            td->results[pulseIndex].phi   = phi ;
        }
        
        if (maxRCS.R < cmplx_mag){
            maxRCS.R = cmplx_mag;
            maxRCS.phi = phi;
            maxRCS.theta = theta;
        }

        free(rnp) ;
    }
    
    
    free( rayAimPoints );
    
    // return to parent thread
    //
    endTimer(&threadTimer, &status) ;
    //.printf("Thread %d has completed in %f secs \n",tid, timeElapsedInSeconds(&threadTimer, &status));
    pthread_exit(NULL) ;
}

void packSinc(SPCmplxD point, SPCmplx *outData, double rdiff, double sampleSpacing, long long nxInData, double * ikernel)
{
    
    double diffFromSincCentre = rdiff / sampleSpacing + nxInData/2;
    int nearestSample = (int)(diffFromSincCentre+0.5);
    double nearestSincSample = nearestSample - diffFromSincCentre;
    int offset;
    int iidx = (int)(nearestSincSample*OVERSAMP);
    if (iidx < 0 ){
        iidx += OVERSAMP ;
        offset = 1;
    }else{
        offset = 0;
    }
    
    for(int x = 0; x < NPOINTS; x++)
    {
        outData[nearestSample+x-NPOINTS/2+offset].r += point.r * ikernel[iidx];
        outData[nearestSample+x-NPOINTS/2+offset].i += point.i * ikernel[iidx];
        iidx += OVERSAMP ;
    }
    
    return ;
}

void sinc_kernel(int oversample_factor, int num_of_points, double resolution, double sampleSpacing, double * data) {
    int n = oversample_factor * num_of_points + 1;
    int x;
    double val;
    double k = M_PI / resolution ;
    double xpos ;
    
    if (oversample_factor == 0 || num_of_points == 0)
    {
        fprintf(stderr, "Either oversample_factor (= %d) or number of points (= %d) is zero in generate_sinc_kernel\n",
                oversample_factor, num_of_points);
        exit(803);
    }
    
    for(x = 0; x < n; x++)
    {
        xpos = (x - n/2) * (sampleSpacing/oversample_factor) ;
        val  = k * xpos ;
        data[x] = (fabs(val) < 1.0e-4 ) ? 1.0 : sin(val)/val;
    }
    
    ham1dx(data, n);
    
    return ;
}

void ham1dx(double * data, int nx)
{
    int x;
    
    double ped = 0.08;
    double a = 1.0 - ped;
    double val;
    
    for(x = 0; x < nx; x++)
    {
        val = M_PI * (-0.5 + (double) x / (double) (nx - 1));
        data[x] *= ped + a * cos(val) * cos(val);
    }
}
