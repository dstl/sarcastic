/***************************************************************************
 *
 *       Module:    threadCore.c
 *      Program:    SARCASTIC
 *   Created by:    Darren on 01/08/2013.
 *                  Copyright (c) 2013 Dstl. All rights reserved.
 *
 *   Description:
 *      This function runs on a single thread and controls a single OpenCL device
 *
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  15/01/2014
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
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
 * THE SOFTWARE IN ITS ENTIRETY OR ANY SUBSTANTIAL PORTION SHOULD NOT BE
 * USED AS A WHOLE OR COMPONENT PART OF DERIVATIVE SOFTWARE THAT IS BEING
 * SOLD TO THE GOVERNMENT OF THE UNITED KINGDOM OF GREAT BRITAIN AND NORTHERN
 * IRELAND.
 *
 ***************************************************************************/
#include "sarcastic.h"
#include <fftw3.h>

#define FASTPULSEBUILD
#define NPOINTS (32)
#define OVERSAMP (512)
#define NOINTERSECTION -1



void packSinc(SPCmplx point, SPCmplx *outData, double rdiff, double sampleSpacing, long long nxInData, double * ikernel);
static void ham1dx(double * data, int nx);
void sinc_kernel(int oversample_factor, int num_of_points, double resolution, double sampleSpacing, double *ikernel);
void generate_gaussian_ray_distribution(double xStdev,double yStdev, SPVector txPos, SPVector GRP, int nx,int ny, Ray *rayArray);
void * safeCalloc(int numItems, int sizeofItems);
void * safeMalloc(int numItems, int sizeofItems);

    
void * devPulseBlock ( void * threadArg ) {
    
    struct threadData *td;
    td = (struct threadData *) threadArg;
    
    cl_device_id devId ;
    cl_context context ;
    cl_command_queue commandQ ;

    int nAzBeam, nElBeam, bounceToShow, nrnpItems, nx, nShadowRays, tid ;
//    int debug, debugX, debugY ;
    int nbounce, nxRay, nyRay, nRays, reflectCount, iray, nShadows, interrogate ;
    int hrs,min,sec;
    int reportN = 500 ;

    double gainRx, PowPerRay, derampRange, derampPhase, pCentDone, sexToGo ;
    double rangeLabel, phasecorr,resolution,sampSpacing, *ikernel, bandwidth ;
    double *ranges, interogRad, intMinR, intMaxR ;

    AABB SceneBoundingBox ;
    rangeAndPower * rnp, irp ;
    Ray *rayArray, *newRays, *reflectedRays, *shadowRays, *LRays, *RRays;
    Hit *hitArray, *newHits, *shadowHits;
    
    SPVector aimdir, RxPos, TxPos, origin, interogPt;
    SPCmplx targ, pcorr, tmp;
    SPImage pulseLine ;
    SPStatus status;
    
    FILE *interogFP ;
   
    struct tm *p;
    Timer threadTimer ;
    time_t current,complete;
    char ct[1000];

    gainRx           = td->gainRx ;
    PowPerRay        = td->PowPerRay ;
    bounceToShow     = td->bounceToShow ;
    nAzBeam          = td->nAzBeam ;
    nElBeam          = td->nElBeam ;
    SceneBoundingBox = td->SceneBoundingBox ;
    bandwidth        = td->chirpRate * td->pulseDuration ;
    tid              = td->devIndex ;
    devId            = td->platform.device_ids[tid] ;
    nx               = (int)td->phd->nx ;
    interrogate      = td->interrogate ;
    interogPt        = td->interogPt ;
    interogRad       = td->interogRad ;
    interogFP        = *(td->interogFP) ;
    
//    debug            = td->debug ;
//    debugX           = td->debugX ;
//    debugY           = td->debugY ;
    
    VECT_CREATE(0, 0, 0, origin);

    fftwf_init_threads();
    
    printf("Compiling OpenCL Kernels...");

    // Create OpenCL context and command queue for this device
    //
    context  = CL_CHECK_ERR(clCreateContext(td->platform.props, 1, &devId, NULL, NULL, &_err));
    commandQ = CL_CHECK_ERR(clCreateCommandQueue(context, devId, 0, &_err));

    // We have the following OpenCL kernels in this thread:
    //  randomRays :  Generates a net of Gaussian distributed rays
    //  stackLessTraverse : Performs a stackless traversal of a KdTree to find the intersection points for rays
    //  reflect : Calculates the reflection ray for a net of rays hitting a surface
    //  buildShadowRays : Calculates the net of rays from a hit point back to the receiver
    //  ( we then call stackLessTraverse again to see which of the shadowRays we can exclude )
    //  reflectPower : Calculates the amount of power, and range at the receiver for a net of rays
    //
    //  Build the kernels now and bail out if any fail to compile
    //
    cl_program randRaysPG,     stackTraversePG,  reflectPG,  buildShadowsPG,  reflectPowerPG ;
    cl_kernel  randRaysKL,     stackTraverseKL,  reflectKL,  buildShadowsKL,  reflectPowerKL ;
    size_t     randRaysLWS[2], stackTraverseLWS, reflectLWS, buildShadowsLWS, reflectPowerLWS ;

    static char *randRaysCode      = "/Users/darren/Development/sarcastic/sarcastic/kernels/randomRays.cl" ;
    static char *stackTraverseCode = "/Users/darren/Development/sarcastic/sarcastic/kernels/stacklessTraverse.cl" ;
    static char *reflectCode       = "/Users/darren/Development/sarcastic/sarcastic/kernels/reflectRays.cl" ;
    static char *buildShadowsCode  = "/Users/darren/Development/sarcastic/sarcastic/kernels/buildShadowRays.cl" ;
    static char *reflectPowCode    = "/Users/darren/Development/sarcastic/sarcastic/kernels/reflectionPowerPO.cl" ;
    
    CL_CHECK(buildKernel(context, randRaysCode,      "randomRays",        devId, 2, &randRaysPG,      &randRaysKL,      randRaysLWS));
    CL_CHECK(buildKernel(context, stackTraverseCode, "stacklessTraverse", devId, 1, &stackTraversePG, &stackTraverseKL, &stackTraverseLWS));
    CL_CHECK(buildKernel(context, reflectCode,       "reflect",           devId, 1, &reflectPG,       &reflectKL,       &reflectLWS));
    CL_CHECK(buildKernel(context, buildShadowsCode,  "buildShadowRays",   devId, 1, &buildShadowsPG,  &buildShadowsKL,  &buildShadowsLWS));
    CL_CHECK(buildKernel(context, reflectPowCode,    "reflectPower",      devId, 1, &reflectPowerPG,  &reflectPowerKL,  &reflectPowerLWS));
  
    // Allocate memory for items that are needed in all kernels
    //
    int                nTriangles       = td->KDT.nTriangles;         // number of triangles in array 'triangles'
    Triangle *         triangles        = td->KDT.triangles;          // Array of triangles of size nTriangles
    int                nTextures        = td->KDT.nTextures;          // number of textures in array 'textures'
    Texture *          textures         = td->KDT.textures;           // Array of textures of size nTextures
    int                nTreeNodes       = td->KDT.nTreeNodes;         // number of nodes in KdTree
    KdData *           KdTree           = td->KDT.KdTree;             // SAH - KdTree to optimise ray traversal through volume
    int                triListDataSize  = td->KDT.triListDataSize;    // size of trianglelist data
    int *              triangleListData = td->KDT.triangleListData;   // array of triangle indices into Triangles
    int                nLeaves          = td->KDT.nLeaves;            // number of leaves (nodes with triangles) in KdTree
    int *              triangleListPtrs = td->KDT.triangleListPtrs;   // array of pointers into triangleListData for each KdTree node

    
    cl_mem dTriangles, dTextures,dKdTree, dtriListData,dtriListPtrs;
    dTriangles   = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(Triangle)*nTriangles, NULL, &_err));
    dTextures    = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(Texture)*nTextures,   NULL, &_err));
    dKdTree      = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(KdData)*nTreeNodes,   NULL, &_err));
    dtriListData = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  triListDataSize,             NULL, &_err));
    dtriListPtrs = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(int)*nLeaves,         NULL, &_err));

    CL_CHECK(clEnqueueWriteBuffer(commandQ, dTriangles,   CL_TRUE, 0, sizeof(Triangle)*nTriangles, triangles,        0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(commandQ, dTextures,    CL_TRUE, 0, sizeof(Texture)*nTextures,   textures,         0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(commandQ, dKdTree,      CL_TRUE, 0, sizeof(KdData)*nTreeNodes,   KdTree,           0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(commandQ, dtriListData, CL_TRUE, 0, triListDataSize,             triangleListData, 0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(commandQ, dtriListPtrs, CL_TRUE, 0, sizeof(int)*nLeaves,         triangleListPtrs, 0, NULL, NULL));
    
    // build a sinc kernel
    //
    resolution  = SIPC_c / (2.0 * bandwidth ) ;
    sampSpacing = SIPC_c / (2.0 * (nx/(td->pulseDuration * td->ADRate)) * bandwidth) ;
    ikernel     = safeCalloc((OVERSAMP * NPOINTS + 1), sizeof(double)) ;

    sinc_kernel(OVERSAMP, NPOINTS, resolution, sampSpacing, ikernel);
    
    printf("...Done\n");
    
    // **** loop  start here
    //
    
    startTimer(&threadTimer, &status) ;

    for (int pulse=0; pulse<td->nPulses; pulse++){
        int pulseIndex = (tid * td->nPulses) + pulse ;
        
        // print out some useful progress information
        //
        if ( td->devIndex == 0 ) {
            if( pulse % reportN == 0 && td->nPulses != 1){
                pCentDone = 100.0*pulse/td->nPulses ;
                printf("Processing pulses %6d - %6d out of %6d [%2d%%]",  pulse*td->nThreads, (pulse+((reportN > td->nPulses) ?  td->nPulses : reportN))*td->nThreads, td->nPulses*td->nThreads,(int)pCentDone);
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
        
        // Set correct parameters for beam to ray trace
        //
        TxPos = td->TxPositions[pulseIndex] ;
        RxPos = td->RxPositions[pulseIndex] ;
        
        // Generate a distribution of nAzbeam x nElbeam rays that originate from the TxPosition aiming at the origin. Use beamMax as the std deviation
        // for the distribution
        //
        rayArray = (Ray *)safeMalloc(nAzBeam*nElBeam, sizeof(Ray));
        
        oclRandomRays(context,commandQ,randRaysKL,nAzBeam,nElBeam,randRaysLWS,td->beamMaxAz,td->beamMaxEl,TxPos, origin, PowPerRay, rayArray);
        
        nbounce = 0;
        nxRay   = nAzBeam ;
        nyRay   = nElBeam ;
        nRays   = nxRay*nyRay;
        
        // Set up deramp range for this pulse
        //
        VECT_MINUS(TxPos, aimdir) ;
        derampRange = VECT_MAG(aimdir);
        derampPhase = -4.0 * SIPC_pi * derampRange * td->oneOverLambda ;
        
        // Use Calloc for rnp as we will be testing for zeroes later on
        //
        rnp = (rangeAndPower *)safeCalloc(nAzBeam*nElBeam*MAXBOUNCES, sizeof(rangeAndPower));
        
        while ( nbounce < MAXBOUNCES &&  nRays != 0){
            
            // Malloc space for hits for this bounce
            //
            hitArray = (Hit *)safeMalloc(nRays,sizeof(Hit));

            // Cast the rays in rayArray through the KdTree using a stackless traversal technique
            // return the hit locations in hitArray
            //
            oclKdTreeHits(context, commandQ, stackTraverseKL, nRays, stackTraverseLWS, dTriangles, dTextures, dKdTree, dtriListData, dtriListPtrs, SceneBoundingBox, rayArray, hitArray);
            
            // Sort out hits
            //
            reflectCount = 0 ;
            iray         = 0 ;
            
            // How many hits occured on this ray cast
            //
            for(int i=0; i<nRays; i++) if ( hitArray[i].trinum != NOINTERSECTION ) reflectCount++ ;
            
            if( reflectCount == 0) break ;
            
            // shrink rayArray and hitArray to get rid of misses
            //
            newRays       = (Ray *)safeMalloc(reflectCount, sizeof(Ray)) ;
            newHits       = (Hit *)safeMalloc(reflectCount, sizeof(Hit)) ;
            reflectedRays = (Ray *)safeMalloc(reflectCount, sizeof(Ray)) ;
            
            for (int i=0; i<nRays; i++) {
                if ( hitArray[i].trinum != NOINTERSECTION ){
                    newRays[iray] = rayArray[i] ;
                    newHits[iray] = hitArray[i] ;
                    iray++ ;
                }
            }
            free(rayArray) ;
            free(hitArray) ;
            rayArray = newRays ;
            hitArray = newHits ;
            nRays    = reflectCount ;
            
            // Build forward scattering rays ready for next turn round the loop
            //
            oclReflect(context, commandQ, reflectKL, dTriangles, dTextures, nTextures, nRays, reflectLWS, rayArray, hitArray, reflectedRays);
            
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
            shadowRays =    (Ray *)safeMalloc(reflectCount, sizeof(Ray)) ;
            shadowHits =    (Hit *)safeMalloc(reflectCount, sizeof(Hit)) ;
            ranges     = (double *)safeMalloc(reflectCount, sizeof(double)) ;
            
            // Build Shadowrays
            //
            oclBuildShadowRays(context, commandQ, buildShadowsKL, buildShadowsLWS, reflectCount, RxPos, reflectedRays, shadowRays, ranges);
            
            
            // Work out which rays have a path back to receiver using stackless traverse kernel
            //
            oclKdTreeHits(context, commandQ, stackTraverseKL, nRays, stackTraverseLWS, dTriangles, dTextures, dKdTree, dtriListData, dtriListPtrs, SceneBoundingBox, shadowRays, shadowHits);
            
            // Shrink the shadowRays to only include those that made it back to the sensor
            // in order to calculate power at sensor we also need the Illumination or LRays
            //
            nShadows = 0 ;
            iray     = 0 ;
            for(int i= 0 ; i<nRays; i++) if ( shadowHits[i].trinum == NOINTERSECTION ) nShadows++ ;
            
            if( nShadows != 0){
                
                newRays       = (Ray *)safeMalloc(nShadows, sizeof(Ray)) ;
                newHits       = (Hit *)safeMalloc(nShadows, sizeof(Hit)) ;
                LRays         = (Ray *)safeMalloc(nShadows, sizeof(Ray)) ;
                RRays         = (Ray *)safeMalloc(nShadows, sizeof(Ray)) ;
                
                for (int i=0; i<nRays; i++) {
                    if ( shadowHits[i].trinum == NOINTERSECTION ){
                        newRays[iray] = shadowRays[i] ;
                        newHits[iray] = hitArray[i] ;
                        LRays[iray]   = rayArray[i] ;
                        RRays[iray]   = reflectedRays[i] ;
                        iray++ ;
                    }
                }
                free(shadowRays) ;
                free(hitArray) ;
                shadowRays  = newRays  ;
                hitArray    = newHits  ;
                nShadowRays = nShadows ;
                
                // For each ray that isn't occluded back to the receiver, calculate the power and put it into rnp.
                //
                oclReflectPower(context, commandQ, reflectPowerKL, reflectPowerLWS, dTriangles, dTextures, hitArray, RxPos, gainRx/(4.0*SIPC_pi), nShadowRays, LRays, RRays, shadowRays, ranges, &(rnp[nAzBeam*nElBeam*nbounce]));
                
                
                // If we are going to interrogate a point then we need to work out the min and max ranges for
                // the point and then check to see if any scatterers are from that range
                //
                if (interrogate) {
                    SPVector intOutRg, intRetRg ;
                    double intRg ;
                    
                    VECT_SUB(interogPt, TxPos, intOutRg);
                    VECT_SUB(RxPos, interogPt, intRetRg);
                    intRg = (VECT_MAG(intOutRg) + VECT_MAG(intRetRg))/2.0;
                    intMinR = intRg - interogRad ;
                    intMaxR = intRg + interogRad ;
                    
                    for ( int i=0; i<nAzBeam*nElBeam; i++){
                        
                        irp = rnp[(nAzBeam*nElBeam*nbounce)+i] ;
                        
                        if ( irp.range > intMinR && irp.range < intMaxR ) {
                            fprintf(interogFP, "%f,\t%e,\t%2d,\t%4d,\t%06.3f,%06.3f,%06.3f\n",irp.range,CMPLX_MAG(irp.Es),nbounce,hitArray[i].trinum,shadowRays[i].org.x,shadowRays[i].org.y,shadowRays[i].org.z);
                        }
                    }
                }
                
                free(LRays);
                free(RRays);
            }
            
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
        for (int i=0; i<nAzBeam*nElBeam*MAXBOUNCES; i++){
            if ((rnp[i].Es.r * rnp[i].Es.r + rnp[i].Es.i * rnp[i].Es.i) != 0 && rnp[i].range !=0) cnt++ ;
        }
        if(cnt > 0){
            nrnpItems = cnt;
        }else{
            printf("ERROR: No intersections on device %d for pulse %d\n",tid, pulseIndex);
            exit(-1);
        }
        
// DEBUG
//        nrnpItems = 1;
// DEBUG
        rnpData_t * rnpData = (rnpData_t *)safeMalloc(nrnpItems, sizeof(rnpData_t));
        
        cnt = 0;
        double totpow = 0.0;
        for (int i=0; i<nAzBeam*nElBeam*MAXBOUNCES; i++){
            if ((rnpData[cnt].Es.r * rnpData[cnt].Es.r + rnpData[cnt].Es.i * rnpData[cnt].Es.i) != 0 && rnp[i].range !=0){
                rnpData[cnt].Es    = rnp[i].Es ;
                rnpData[cnt].rdiff = rnp[i].range - derampRange;
                totpow += rnpData[cnt].Es.r * rnpData[cnt].Es.r + rnpData[cnt].Es.i * rnpData[cnt].Es.i ;
                cnt++ ;
            }
        }
//        printf("Total power : %e, # of intersecting rays : %d, powerperray: %f db (%f)\n",totpow,nrnpItems,10*log(PowPerRay), PowPerRay);
        
        ////////// DEBUG CODE to inject a single point target
        
//        float xpos = -5 ;
//        float ypos = 0 ;
//        SPVector tp ; // Target Point in metres from scene centre in E,N,Alt system
//        SPVector rvect1, rvect2 ;
//        SPVector R,sensor_pos, sc_pos, R_,KDP_,JDP_,IDP_;
//        SPVector rng_dir, azi_dir ;
//        sensor_pos = td->TxPositions[td->nPulses/2] ;
//        VECT_CREATE(0, 0, 0, sc_pos);
//        VECT_SUB(sensor_pos, sc_pos, R);
//        VECT_UNIT(R, R_);
//        VECT_CREATE(0, 0, 1, KDP_);     // == zbar
//        VECT_PROJ(R_, KDP_, JDP_);
//        VECT_CROSS(JDP_, KDP_, IDP_);
//        VECT_SCMULT(JDP_, -1*ypos, rng_dir);
//        VECT_SCMULT(IDP_,  xpos, azi_dir);
//        VECT_ADD(azi_dir, rng_dir, tp) ;
//        VECT_SUB(tp, TxPos, rvect1);
//        VECT_SUB(tp, RxPos, rvect2) ;
//        double r1,r2 ;
//        r1 = VECT_MAG(rvect1);
//        r2 = VECT_MAG(rvect2);
//        int phase_sign = -1.0;
//        rnpData[0].power = 1.0 ;
//        rnpData[0].range = (r1 + r2) / 2.0 ;
//        rnpData[0].rdiff = phase_sign * (derampRange - rnpData[0].range) ;
        
        ////////// DEBUG CODE to inject a single point target
        
        
        im_create(&pulseLine, ITYPE_CMPL_FLOAT, nx, 1, 1.0, 1.0, &status);
        
        for (int i=0; i<nrnpItems; i++){
            
            double phse = CMPLX_PHASE(rnpData[i].Es) - derampPhase ;
            targ.r      = CMPLX_MAG(rnpData[i].Es) * cos(phse) ;
            targ.i      = CMPLX_MAG(rnpData[i].Es) * sin(phse) ;
            
//            double phse  = (-4.0 * SIPC_pi * rnpData[i].rdiff * td->oneOverLambda) ;
//            targ.r = rnpData[i].power*cos(phse) ;
//            targ.i = rnpData[i].power*sin(phse) ;
            
            rangeLabel = (rnpData[i].rdiff/sampSpacing) + (pulseLine.nx / 2) ;
            
            if (rangeLabel > NPOINTS/2 && rangeLabel < nx - NPOINTS) {
                packSinc(targ, pulseLine.data.cmpl_f, rnpData[i].rdiff, sampSpacing, pulseLine.nx, ikernel);
            }
        }
        
        // perform phase correction to account for deramped jitter in receiver timing
        //
        phasecorr = (((td->Fx0s[pulseIndex] - td->freq_centre) / td->FxSteps[pulseIndex])) * 2.0 * M_PI / pulseLine.nx;
        
        for(int x = 0; x < pulseLine.nx; x++) {
            pcorr.r = td->amp_sf0[pulseIndex] * cos(phasecorr * (x - pulseLine.nx/2)) ;
            pcorr.i = td->amp_sf0[pulseIndex] * sin(phasecorr * (x - pulseLine.nx/2)) ;
            
            if(CMPLX_MAG(pulseLine.data.cmpl_f[x]) == 0.0 ){
                tmp.r = tmp.i = 0.0 ;
            }else{
                CMPLX_MULT(pulseLine.data.cmpl_f[x], pcorr, tmp);
            }

            pulseLine.data.cmpl_f[x] = tmp;
        }
        
        im_circshift(&pulseLine, -(pulseLine.nx/2), 0, &status);
        im_fftw(&pulseLine, FFT_X_ONLY+FWD+NOSCALE, &status);
        im_insert(&pulseLine, 0, pulseIndex, td->phd, &status) ;
        im_destroy(&pulseLine, &status) ;
        
        free(rnpData);
        free(rnp);

    } // end of pulse loop
    free(ikernel);
    
    // Clear down OpenCL allocations
    //
    clReleaseKernel(randRaysKL);
    clReleaseKernel(stackTraverseKL);
    clReleaseKernel(reflectKL);
    clReleaseKernel(buildShadowsKL);
    clReleaseKernel(reflectPowerKL);
    clReleaseProgram(randRaysPG) ;
    clReleaseProgram(stackTraversePG) ;
    clReleaseProgram(reflectPG) ;
    clReleaseProgram(buildShadowsPG) ;
    clReleaseProgram(reflectPowerPG) ;
    clReleaseCommandQueue(commandQ);
    clReleaseContext(context);
    
    // return to parent thread
    //
    pthread_exit(NULL) ;
}

void packSinc(SPCmplx point, SPCmplx *outData, double rdiff, double sampleSpacing, long long nxInData, double * ikernel)
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

static void
ham1dx(double * data, int nx)
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

void * safeCalloc(int numItems, int sizeofItems){
    void * ret ;
    ret = calloc(numItems, sizeofItems);
    if (ret == NULL) {
        printf("*** Failed to safeCalloc %d bytes at line %d, file %s\n",sizeofItems*numItems, __LINE__, __FILE__);
        exit(-98);
    }
    return ret ;
}

void * safeMalloc(int numItems, int sizeofItems){
    void * ret ;
    ret = malloc(sizeofItems*numItems);
    if (ret == NULL){
        printf("*** Failed to safeMalloc %d bytes at line %d, file %s\n",sizeofItems*numItems, __LINE__, __FILE__);
        exit(-99);
    }
//    printf(" ++++ Malloced %d bytes (pointer %p)\n",sizeofItems*numItems,ret);

    return ret ;
}