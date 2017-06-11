
#include <fftw3.h>
extern "C" {
#include "boxMullerRandom.h"
#include "RCS.h"
}
#include "ranf.h"
#include "threadCore.hpp"
#include <sys/time.h>
#include "rayTrace.hpp"
#include "cpuPOField.hpp"
#include <SILib2/SILib2.h>

#define NOINTERSECTION -1
//#define TOTALRCSINPULSE


void * devPulseBlock ( void * threadArg ) {
    
    struct threadData *td;
    td = (struct threadData *) threadArg;

    int nAzBeam, nElBeam, bounceToShow, nrnpItems, nx=0, nShadowRays ;
    int nbounce, nxRay, nyRay, nRays, reflectCount, iray, nShadows, interrogate ;
    int maxRaysPerBounce ;

    int hrs,min,sec;
    int reportN = 10 ;
    int tid, nThreads, nPulses, startPulse;
    
    double gainRx, PowPerRay, derampRange, derampPhase, pCentDone, sexToGo ;
    double rangeLabel, phasecorr,resolution,sampSpacing=0, *ikernel=NULL, bandwidth ;
    double *ranges, *newranges, interogRad, intMinR, intMaxR, k ;
    
    AABB SceneBoundingBox ;
    rangeAndPower * rnp, irp ;
    Ray *rayArray, *newRays, *reflectedRays, *shadowRays, *LRays, *RRays;
    Hit *hitArray, *newHits, *shadowHits;
    
    SPVector aimdir, RxPos, TxPos, origin, interogPt;
    SPCmplx pcorr, tmp;
    SPCmplxD targ ;
    SPImage pulseLine ;
    SPStatus status;
    
    FILE *interogFP ;
    
    struct tm *p;
    Timer threadTimer ;
    time_t current,complete;
    char ct[1000];
    bool dynamicScene = false;
    kdTree::KdData *node;
    ATS *accelTriangles = NULL;
    int Ropes[6] ;
    AABB sceneAABB ;
    
    kdTree::KdData * tree = NULL;
    int treeSize;
    int nTriangles=0;
    
    tid              = td->tid ;
    nThreads         = td->nThreads ;
    startPulse       = td->startPulse ;
    nPulses          = td->nPulses ;
    CPHDHeader *hdr  = td->cphdhdr ;
    gainRx           = td->gainRx ;
    PowPerRay        = td->PowPerRay ;
    bounceToShow     = td->bounceToShow ;
    nAzBeam          = td->nAzBeam ;
    nElBeam          = td->nElBeam ;
    SceneBoundingBox = td->SceneBoundingBox ;
    bandwidth        = hdr->chirp_gamma * hdr->pulse_length ;
    interrogate      = td->interrogate ;
    interogPt        = td->interogPt ;
    interogRad       = td->interogRad ;
    interogFP        = *(td->interogFP) ;
    k                = 2 * SIPC_pi / (SIPC_c / hdr->freq_centre) ;
    
    VECT_CREATE(0, 0, 0, origin);
    
    im_init_status(status, 0) ;

    
    if(td->phd != NULL){
        // build a sinc kernel
        //
        nx = (int)td->phd->nx ;
        resolution  = SIPC_c / (2.0 * bandwidth ) ;
        sampSpacing = SIPC_c / (2.0 * (nx/(hdr->pulse_length * hdr->clock_speed)) * bandwidth) ;
        ikernel     = (double *)sp_calloc((OVERSAMP * NPOINTS + 1), sizeof(double)) ;
        
        sinc_kernel(OVERSAMP, NPOINTS, resolution, sampSpacing, ikernel);
    }
    
    if(td->moverMesh->triangles.size() != 0){
        dynamicScene = true;
    }
    // **** loop  start here
    //
    
    startTimer(&threadTimer, &status) ;
    SPVector *rayAimPoints = NULL ;
    
    double t0 = hdr->pulses[0].sat_tx_time ;
    
    // set up some timers
    //
    Timer pulseTimer, buildRaysTimer, shootRayTimer, reflectTimer, buildShadowRayTimer, packPulseTimer, POTimer;
    double pulseDur, buildRaysDur, shootRaysDur, reflectDur, buildShadDur, packPulseDur, PODur ;
    pulseDur = buildRaysDur = shootRaysDur = reflectDur = buildShadDur = packPulseDur = PODur = 0.0;
    
    for (int pulse=0; pulse<nPulses; pulse++){
        startTimer(&pulseTimer, &status);
        int pulseIndex = pulse + startPulse ;
#ifndef TOTALRCSINPULSE
        // print out some useful progress information
        //
        
        if(tid == 0){
            if( pulse % reportN == 0 && nPulses != 1){
                pCentDone = 100.0*pulse/nPulses ;
                printf("Processing pulses %6d - %6d out of %6d [%2d%%]",
                       pulse, (pulse+((reportN > (nPulses)) ?  (nPulses*nThreads) : reportN*nThreads)), nPulses*nThreads ,(int)pCentDone);
                if(pulse != 0 ){
                    current  = time(NULL);
                    sexToGo  = estimatedTimeRemaining(&threadTimer, pCentDone, &status);
                    hrs      = (int)floor(sexToGo/60/60) ;
                    min      = (int)floor(sexToGo / 60) - (hrs*60);
                    sec      = (int)sexToGo % 60;
                    complete = current + sexToGo ;
                    p        = localtime(&complete) ;
                    strftime(ct, 1000, "%a %b %d %H:%M", p);
                    printf("  ETC %s (in %2.0dh:%2.0dm:%2.0ds) \n",ct,hrs,min,sec);
                }else{
                    printf("  Calculating ETC...\n");
                }
            }
        }
      
#endif //TOTALRCSINPULSE
        
        // Set correct parameters for beam to ray trace
        //
        TxPos = hdr->pulses[pulseIndex].sat_ps_tx ;
        RxPos = hdr->pulses[pulseIndex].sat_ps_rx ;
        TriangleMesh newMesh ;
        if(dynamicScene){
            // Construct a new kdTree for this pulse taking into account any movers
            //
            double t = hdr->pulses[pulseIndex].sat_tx_time  - t0;
            
            SPVector S0,S1,S2,S;
            VECT_CREATE(0.0, 0.0, 0.0, S0) ;  // Translation
            VECT_CREATE(0.0, 0.0, 0.0, S1) ;  // Velocity
            VECT_CREATE(0.0, 0.0, 0.0, S2) ;  // Acceleration
            
            
            // move the movers to the location for this pulse
            //
            // Linear Motion
            //
            S.x = S0.x + (S1.x * t) + (0.5 * S2.x * t * t) ;
            S.y = S0.y + (S1.y * t) + (0.5 * S2.y * t * t) ;
            S.z = S0.z + (S1.z * t) + (0.5 * S2.z * t * t) ;
            
            // Harmonic Motion
            //
            double A = 0.0125 ;
            double f = 0.33 ;
            S.x = A * sin(2 * SIPC_pi * t * f ) ;
            double Ts = 1/f;
            double frac = (t / Ts) - (int)(t / Ts) ;
            S.x = A * frac ;
            
            // Modification for top hat motion
            //
//            if ( S.x > 0.0) {
//                S.x = A;
//            }else{
//                S.x = -A ;
//            }

            //if ( S.x > A/2) {
            //    S.x = A;
            //}else if (S.x < -A/2){
            //    S.x = -A ;
            //}else{
            //    S.x = 0;
            //}
            
            TriangleMesh mesh_t = *(td->moverMesh) ;
            for(int i=0; i<mesh_t.vertices.size(); ++i){
                mesh_t.vertices[i].x += S.x ;
                mesh_t.vertices[i].y += S.y ;
                mesh_t.vertices[i].z += S.z ;
            }
            
            // Add movers to base scene
            //
            newMesh = td->sceneMesh->add(mesh_t) ;
        }else{
            newMesh = *(td->sceneMesh) ;
        }
        
        // Build kdTree if required
        //
        if (dynamicScene) {
            if (!dynamicScene) {
                endTimer(&threadTimer, &status);
            }
            kdTree::buildTree(newMesh, &tree, &treeSize, (kdTree::TREEOUTPUT)(kdTree::OUTPUTNO)) ;
            accelerateTriangles(&newMesh,&accelTriangles) ;
            
            // Initialise the tree and build ropes and boxes to increase efficiency when traversing
            //
            node = &(tree[0]) ;
            
            // build Ropes and Boxes
            //
            sceneAABB = tree[0].brch.aabb ;
            for(int i=0; i<6; i++) Ropes[i] = NILROPE;
            BuildRopesAndBoxes(node, Ropes, sceneAABB, tree);
            if (!dynamicScene) {
                startTimer(&threadTimer, &status);
            }
            
        }else{
            tree = *(td->tree) ;
            treeSize = td->treesize ;
            accelTriangles = *(td->accelTriangles) ;
        }

        // Allocate memory for items that are needed in all kernels
        //
        nTriangles   = (int)newMesh.triangles.size() ; // number of triangles in array 'triangles'
        
        // Generate a distribution of nAzbeam x nElbeam rays that originate from the TxPosition aiming at the origin. Use beamMax as the std deviation
        // for the distribution
        //
        nbounce = 0;
        nxRay   = nAzBeam ;
        nyRay   = nElBeam ;
        nRays   = nxRay*nyRay; // nRays is the number of rays in each bounce and is rewritten after each reflection
        startTimer(&buildRaysTimer, &status) ;
        buildRays(&rayArray, &nRays, nAzBeam, nElBeam, &newMesh, TxPos, PowPerRay, td->SceneBoundingBox, &rayAimPoints);
        endTimer(&buildRaysTimer, &status);
        buildRaysDur += timeElapsedInMilliseconds(&buildRaysTimer, &status);
        maxRaysPerBounce = nRays;  // Use this for memory as nxRay/nyRay may be incorrect if buildRays set to triangle centres
        
        // Set up deramp range for this pulse
        //
        VECT_MINUS(TxPos, aimdir) ;
        derampRange = VECT_MAG(aimdir);
        derampPhase = -4.0 * SIPC_pi * derampRange * hdr->freq_centre / SIPC_c ; ;
        
        // Use Calloc for rnp as we will be testing for zeroes later on
        //
        rnp = (rangeAndPower *)sp_calloc(maxRaysPerBounce*MAXBOUNCES, sizeof(rangeAndPower));
        
        while ( nbounce < MAXBOUNCES &&  nRays != 0){
            
            // Malloc space for hits for this bounce
            //
            hitArray = (Hit *)sp_malloc(nRays * sizeof(Hit));
            
            // Cast the rays in rayArray through the KdTree using a stackless traversal technique
            // return the hit locations in hitArray
            //
//            oclKdTreeHits(context, commandQ, stackTraverseKL, nRays, stackTraverseLWS, dTriangles, dKdTree, dtriListData, dtriListPtrs, SceneBoundingBox, rayArray, hitArray);
            
            startTimer(&shootRayTimer, &status);
            shootRay(tree, accelTriangles, nRays, rayArray, hitArray) ;
            endTimer(&shootRayTimer, &status);
            shootRaysDur += timeElapsedInMilliseconds(&shootRayTimer, &status);
            
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
            startTimer(&reflectTimer, &status) ;
            reflect(nRays, rayArray, hitArray, accelTriangles, reflectedRays);
            endTimer(&reflectTimer, &status);
            reflectDur += timeElapsedInMilliseconds(&reflectTimer, &status) ;
            
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
            startTimer(&buildShadowRayTimer, &status);
            buildShadowRays(nRays, RxPos, reflectedRays, shadowRays, ranges) ;
            endTimer(&buildShadowRayTimer, &status);
            buildShadDur += timeElapsedInMilliseconds(&buildShadowRayTimer, &status);
            
            // Work out which rays have a path back to receiver using stackless traverse kernel
            //
//            oclKdTreeHits(context, commandQ, stackTraverseKL, nRays, stackTraverseLWS, dTriangles , dKdTree, dtriListData, dtriListPtrs, SceneBoundingBox, shadowRays, shadowHits);
            startTimer(&shootRayTimer, &status);
            shootRay(tree, accelTriangles, nRays, shadowRays, shadowHits) ;
            endTimer(&shootRayTimer, &status);
            shootRaysDur += timeElapsedInMilliseconds(&shootRayTimer, &status);
            
            // Shrink the shadowRays to only include those that made it back to the sensor
            // in order to calculate power at sensor we also need the Illumination or LRays
            //
            nShadows = 0 ;
            double *facing = (double *)sp_malloc(sizeof(double) * nRays);
            int *hitsOnEachTri  = (int *)sp_calloc(nTriangles, sizeof(int));
            
            SPVector normal;
            for(int i=0; i<nRays; i++){
                normal = newMesh.triangles[hitArray[i].trinum].N.asSPVector() ;
                facing[i] =VECT_DOT(normal, shadowRays[i].dir);
            }
            for(int i= 0 ; i<nRays; i++){
                hitsOnEachTri[ hitArray[i].trinum ]++ ;
                if ( (shadowHits[i].trinum == NOINTERSECTION) && (facing[i] > 0.1) && (hitsOnEachTri[hitArray[i].trinum] <= 1) ){
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
                    if ( (shadowHits[i].trinum == NOINTERSECTION) && (facing[i] > 0.1) && (hitsOnEachTri[hitArray[i].trinum] <= 1) ){
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
                startTimer(&POTimer, &status);
                cpuPOField(&newMesh, hitArray, nShadowRays, LRays, shadowRays, RxPos, k, ranges, gainRx, nbounce+1, &(rnp[maxRaysPerBounce*nbounce])) ;
                endTimer(&POTimer, &status);
                PODur += timeElapsedInMilliseconds(&POTimer, &status);
                
                
                // If we are going to interrogate a point then we need to work out the min and max ranges for
                // the point and then check to see if any scatterers are from that range
                //
                if (interrogate) {
                    SPVector intOutRg, intRetRg ;
                    double intRg ;
                    
                    VECT_SUB(interogPt, TxPos, intOutRg);
                    VECT_SUB(RxPos, interogPt, intRetRg);
                    intRg = (VECT_MAG(intOutRg) + VECT_MAG(intRetRg))/2.0;
                    intMinR = intRg - interogRad/2.0 ;
                    intMaxR = intRg + interogRad/2.0 ;
                    
                    for ( int i=0; i<nShadowRays; i++){
                        
                        irp = rnp[(maxRaysPerBounce*nbounce)+i] ;
                        
                        if ( irp.range > intMinR && irp.range < intMaxR ) {
                            fprintf(interogFP, "%f,\t%e,\t%2d,\t%4d,\t%06.3f,%06.3f,%06.3f\n",irp.range,CMPLX_MAG(irp.Es),nbounce,hitArray[i].trinum,shadowRays[i].org.x,shadowRays[i].org.y,shadowRays[i].org.z);
                        }
                    }
                }
                
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
        
        if(td->phd != NULL){
            startTimer(&packPulseTimer, &status);

            int cnt=0;
            nrnpItems = 0;
            nRays = maxRaysPerBounce ;
            for (int i=0; i<maxRaysPerBounce*MAXBOUNCES; i++){
                if ((rnp[i].Es.r * rnp[i].Es.r + rnp[i].Es.i * rnp[i].Es.i) != 0 && rnp[i].range !=0) cnt++ ;
            }
            if(cnt > 0){
                nrnpItems = cnt;
            }else{
                printf("ERROR: No intersections for pulse %d\n",pulseIndex);
                exit(-1);
            }
            
            // DEBUG
            //        nrnpItems = 1;
            // DEBUG
            rnpData_t * rnpData = (rnpData_t *)sp_malloc(nrnpItems * sizeof(rnpData_t));
            cnt = 0;
            for (int i=0; i<maxRaysPerBounce*MAXBOUNCES; i++){
                if (CMPLX_MAG(rnp[i].Es) != 0 && rnp[i].range !=0){
                    rnpData[cnt].Es    = rnp[i].Es ;
                    rnpData[cnt].rdiff = rnp[i].range - derampRange;
                    cnt++ ;
                }
            }
            
            im_init(&pulseLine, &status);
            im_create(&pulseLine, ITYPE_CMPL_FLOAT, nx, 1, 1.0, 1.0, &status);
            double phse;
            SPCmplxD targtot = {0.0,0.0};
            for (int i=0; i<nrnpItems; i++){
                rangeLabel  = (rnpData[i].rdiff/sampSpacing) + (pulseLine.nx / 2) ;
                if (rangeLabel > NPOINTS/2 && rangeLabel < nx - NPOINTS) {
                    
                    phse        = CMPLX_PHASE(rnpData[i].Es) - derampPhase ;
                    targ.r      = CMPLX_MAG(rnpData[i].Es) * cos(phse) ;
                    targ.i      = CMPLX_MAG(rnpData[i].Es) * sin(phse) ;
                    
                    packSinc(targ, pulseLine.data.cmpl_f, rnpData[i].rdiff, sampSpacing, pulseLine.nx, ikernel);
                    
                    CMPLX_ADD(rnpData[i].Es, targtot, targtot);
                }
            }
            
#ifdef TOTALRCSINPULSE
            double cmplx_mag = RCS(PowPerRay, CMPLX_MAG(targtot), derampRange, derampRange);
            double oneOverLambda = td->cphdhdr->freq_centre / SIPC_c ;
            if(pulse%1==0)printf("%d, %e\n", pulseIndex+td->startPulse,cmplx_mag);
            printf("Total RCS for pulse %d is %f m^2 (%f dB m^2)\n",pulse,cmplx_mag,10*log10(cmplx_mag));
            printf("For comparison: \n");
            printf("    1m^2 flat plate : %f (%f dB m^2)\n",1*SIPC_pi*4*oneOverLambda*oneOverLambda,10*log10(4*SIPC_pi*oneOverLambda*oneOverLambda)); // 8620.677
            printf("    1m dihedral     : %f (%f dB m^2)\n",8*SIPC_pi*oneOverLambda*oneOverLambda,10*log10(8*SIPC_pi*oneOverLambda*oneOverLambda));   // 17241.354
            printf("    1m trihedral    : %f (%f dB m^2)\n",12*SIPC_pi*oneOverLambda*oneOverLambda,10*log10(12*SIPC_pi*oneOverLambda*oneOverLambda)); // 25862.031
#endif // TOTALRCSINPULSE
            // perform phase correction to account for deramped jitter in receiver timing
            //
            phasecorr = (((hdr->pulses[pulseIndex].fx0 - hdr->freq_centre) / hdr->pulses[pulseIndex].fx_step_size)) * 2.0 * M_PI / pulseLine.nx;
            
            for(int x = 0; x < pulseLine.nx; x++) {
                pcorr.r =  cos(phasecorr * (x - pulseLine.nx/2)) / (hdr->pulses[pulseIndex].amp_sf0);
                pcorr.i =  sin(phasecorr * (x - pulseLine.nx/2)) / (hdr->pulses[pulseIndex].amp_sf0);
                
                if(CMPLX_MAG(pulseLine.data.cmpl_f[x]) == 0.0 ){
                    tmp.r = tmp.i = 0.0 ;
                }else{
                    CMPLX_MULT(pulseLine.data.cmpl_f[x], pcorr, tmp);
                }
                
                pulseLine.data.cmpl_f[x] = tmp;
            }
            im_circshift(&pulseLine, -(pulseLine.nx/2), 0, &status);
            im_fftw(&pulseLine, (FFTMODE)(FFT_X_ONLY+FWD+NOSCALE), &status);
            im_insert(&pulseLine, 0, pulseIndex, td->phd, &status) ;
            im_destroy(&pulseLine, &status) ;
            
            free(rnpData);
            endTimer(&packPulseTimer, &status);
            packPulseDur+= timeElapsedInMilliseconds(&packPulseTimer, &status);
            
        } // end of pulse loop
        
        free(rnp) ;
        if(dynamicScene){
            free(tree) ;
            delete accelTriangles ;
        }
        endTimer(&pulseTimer, &status);
        pulseDur += timeElapsedInMilliseconds(&pulseTimer, &status) ;
    }
    
    if(td->phd != NULL){
        free(ikernel);
    }
    
    free( rayAimPoints );
    if(!dynamicScene){
        free(tree);
        delete accelTriangles ;
    }
    
    // print timer Summary
    //
    if(tid==0){
        printf("           Timing Summary\n");
        printf("==========================================\n");
        printf("Time spent building rays          : %8.2f ms\n",buildRaysDur * nThreads);
        printf("Time spent shooting rays          : %8.2f ms\n",shootRaysDur * nThreads);
        printf("Time spent reflecting rays        : %8.2f ms\n",reflectDur * nThreads);
        printf("Time spent building shadowRays    : %8.2f ms\n",buildShadDur * nThreads);
        printf("Time spent packing the pulse      : %8.2f ms\n",packPulseDur * nThreads);
        printf("Time spent calculating PO Fields  : %8.2f ms\n",PODur * nThreads);
        printf("Time data processing in loop      : %8.2f ms\n",(pulseDur - (buildRaysDur+shootRaysDur+reflectDur+buildShadDur+packPulseDur+PODur)) * nThreads);
        printf("Number of loops                   : %8d pulses\n",nPulses * nThreads);
        printf("Total time spent on pulses        : %8.2f ms\n",pulseDur * nThreads);
        printf("------------------------------------------\n");
        printf("  Total time per pulse            : %8.2f ms\n", pulseDur/(nPulses * nThreads));
    }

//    return (NULL);
    // return to parent thread
    //
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
