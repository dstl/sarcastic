/***************************************************************************
 *
 *       Module:    threadCore.c
 *      Program:    GOSS
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
#include "GOSS.h"
#include <fftw3.h>

#define FASTPULSEBUILD
#define NPOINTS (32)
#define OVERSAMP (512)
#define NOINTERSECTION -1



void packSinc(SPCmplx point, SPCmplx *outData, double rdiff, double sampleSpacing, long long nxInData, double * ikernel);
static void ham1dx(double * data, int nx);
void sinc_kernel(int oversample_factor, int num_of_points, double resolution, double sampleSpacing, double *ikernel);
void generate_gaussian_ray_distribution(double xStdev,double yStdev, SPVector txPos, SPVector GRP, int nx,int ny, Ray *rayArray);

void * devPulseBlock ( void * threadArg ) {
    
    struct threadData *td;
    td = (struct threadData *) threadArg;
    
    cl_device_id devId ;
    cl_context context ;
    cl_command_queue commandQ ;
    cl_int err ;
    rangeAndPower * rnp ;
    SPVector aimdir ;
    
    int nAzBeam             = td->nAzBeam ;
    int nElBeam             = td->nElBeam ;
    SPVector RxPos;
    SPVector TxPos;
    double gainRx           = td->gainRx ;
    double PowPerRay        = td->PowPerRay ;
    AABB SceneBoundingBox   = td->SceneBoundingBox ;
    int bounceToShow        = td->bounceToShow ;
    int debug               = td->debug;
    int debugX              = td->debugX;
    int debugY              = td->debugY;
    int nrnpItems;
    double derampRange ;
    SPCmplx targ, pcorr, tmp;
    SPImage pulseLine ;
    double rangeLabel, phasecorr,resolution,sampSpacing, *ikernel ;
    double bandwidth = td->chirpRate * td->pulseDuration ;
    
    Timer threadTimer ;
    SPStatus status;
    startTimer(&threadTimer, &status) ;
    int hrs,min,sec;
    time_t current,complete;
    struct tm *p;
    int reportN = 500 ;
    char ct[1000];
    double pCentDone ;
    double sexToGo ;
    int nx;

    int nbounce ;
    int nxRay   = nAzBeam ;
    int nyRay   = nElBeam ;
    int nRays   = nxRay*nyRay;
    int reflectCount ;
    int iray;
    int nShadows ;
    SPVector origin;
    
    VECT_CREATE(0, 0, 0, origin);
    Ray *rayArray, *newRays, *reflectedRays;
    Hit *hitArray, *newHits;
    Ray *shadowRays, *LRays, *RRays;
    Hit *shadowHits;
    double *ranges ;
    double *tempDists ;
    int nShadowRays;
    
    fftwf_init_threads();
    
    printf("Using FAST pulse building routines...\n");

    int tid ;
    tid   = td->devIndex ;
    devId = td->platform.device_ids[tid] ;
    nx    = (int)td->phd->nx ;

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

    static char *randRaysCode      = "/Users/darren/Development/GOSS/GOSS/randomRays.cl" ;
    static char *stackTraverseCode = "/Users/darren/Development/GOSS/GOSS/stacklessTraverse.cl" ;
    static char *reflectCode       = "/Users/darren/Development/GOSS/GOSS/reflectRays.cl" ;
    static char *buildShadowsCode  = "/Users/darren/Development/GOSS/GOSS/buildShadowRays.cl" ;
    static char *reflectPowCode    = "/Users/darren/Development/GOSS/GOSS/reflectionPower.cl" ;
    
    CL_CHECK(buildKernel(context, randRaysCode,      "randomRays",        devId, 2, &randRaysPG,      &randRaysKL,      randRaysLWS));
    CL_CHECK(buildKernel(context, stackTraverseCode, "stacklessTraverse", devId, 1, &stackTraversePG, &stackTraverseKL, &stackTraverseLWS));
    CL_CHECK(buildKernel(context, reflectCode,       "reflect",           devId, 1, &reflectPG,       &reflectKL,       &reflectLWS));
    CL_CHECK(buildKernel(context, buildShadowsCode,  "buildShadowRays",   devId, 1, &buildShadowsPG,  &buildShadowsKL,  &buildShadowsLWS));
    CL_CHECK(buildKernel(context, reflectPowCode,    "reflectPower",      devId, 1, &reflectPowerPG,  &reflectPowerKL,  &reflectPowerLWS));
  

    // build a sinc kernel
    //
    resolution  = SIPC_c / (2.0 * bandwidth ) ;
    sampSpacing = SIPC_c / (2.0 * (nx/(td->pulseDuration * td->ADRate)) * bandwidth) ;
    ikernel = calloc((OVERSAMP * NPOINTS + 1), sizeof(double));
    if (ikernel == NULL){
        fprintf(stderr, "Failed to calloc ikernel %s:%d\n", __FILE__, __LINE__);
        exit(802);
    }
    sinc_kernel(OVERSAMP, NPOINTS, resolution, sampSpacing, ikernel);
    
    // **** loop  start here
    //
    for (int pulse=0; pulse<td->nPulses; pulse++){
        int pulseIndex = (tid * td->nPulses) + pulse ;
        
        // print out some useful progress information
        //
        if ( td->devIndex == 0 ) {
            if( pulse % reportN == 0 && td->nPulses != 1){
                pCentDone = 100.0*pulse/td->nPulses ;
                printf("Processing pulses %6d - %6d out of %6d [%2d%%]",  pulse*td->nThreads, (pulse+reportN)*td->nThreads, td->nPulses*td->nThreads,(int)pCentDone);
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
        
        // Generate a distribution of nAzbeam x nElbeam rays that originate form the TxPosition aiming at the origin. Use beamMax as the std deviation
        // for the distribution
        //
        rayArray = (Ray *)malloc(sizeof(Ray)*nAzBeam*nElBeam);
        
        oclRandomRays(context,commandQ,randRaysKL,nAzBeam,nElBeam,randRaysLWS,td->beamMaxAz,td->beamMaxEl,TxPos, origin, PowPerRay, rayArray);
        
        nbounce = 0;
        nxRay   = nAzBeam ;
        nyRay   = nElBeam ;
        nRays   = nxRay*nyRay;
        
        // Use Calloc for rnp as we will be testing for zeroes later on
        //
        rnp = (rangeAndPower *)calloc(nAzBeam*nElBeam*MAXBOUNCES,sizeof(rangeAndPower)) ;
        
        while ( nbounce < MAXBOUNCES &&  nRays != 0){
            
            // Malloc space for hits for this bounce
            //
            hitArray = (Hit *)malloc(sizeof(Hit)*nRays);
            if(hitArray == NULL){ printf("MALLOC ERROR \n");exit(-42);}

            // Cast the rays in rayArray through the KdTree using a stackless traversal technique
            // return the hit locations in hitArray
            //
            oclKdTreeHits(context, commandQ, stackTraverseKL, nRays, stackTraverseLWS, td->KDT, SceneBoundingBox, rayArray, hitArray) ;
            
            // Sort out hits
            //
            reflectCount = 0 ;
            iray         = 0 ;
            
            // How many hits occured on this ray cast
            //
            for(int i=0; i<nRays; i++) if ( hitArray[i].trinum != NOINTERSECTION ) reflectCount++ ;
            
            if( reflectCount != 0 ){
            
                // shrink rayArray and hitArray to get rid of misses
                //
                newRays       = (Ray *)malloc(sizeof(Ray) * reflectCount) ;
                newHits       = (Hit *)malloc(sizeof(Hit) * reflectCount) ;
                reflectedRays = (Ray *)malloc(sizeof(Ray) * reflectCount) ;

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
                oclReflect(context, commandQ, reflectKL, td->KDT, nRays, reflectLWS, rayArray, hitArray, reflectedRays);
            
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
                shadowRays =    (Ray *)malloc(sizeof(Ray)    * reflectCount) ;
                shadowHits =    (Hit *)malloc(sizeof(Hit)    * reflectCount) ;
                ranges     = (double *)malloc(sizeof(double) * reflectCount) ;
                
                // Build Shadowrays
                //
                oclBuildShadowRays(context, commandQ, buildShadowsKL, buildShadowsLWS, reflectCount, RxPos, reflectedRays, shadowRays, ranges);
            
            
                // Work out which rays have a path back to receiver using stackless traverse kernel
                //
                oclKdTreeHits(context, commandQ, stackTraverseKL, nRays, stackTraverseLWS, td->KDT, SceneBoundingBox, shadowRays, shadowHits) ;
            
                // Shrink the shadowRays to only include those that made it back to the sensor
                // in order to calculate power at sensor we also need the Illumination or LRays
                //
                nShadows = 0 ;
                iray     = 0 ;
                for(int i= 0 ; i<nRays; i++) if ( shadowHits[i].trinum == NOINTERSECTION ) nShadows++ ;
            
                if( nShadows != 0){
                    
                    newRays       = (Ray *)malloc(sizeof(Ray) * nShadows) ;
                    newHits       = (Hit *)malloc(sizeof(Hit) * nShadows) ;
                    LRays         = (Ray *)malloc(sizeof(Ray) * nShadows) ;
                    RRays         = (Ray *)malloc(sizeof(Ray) * nShadows) ;
                    
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
                    oclReflectPower(context, commandQ, reflectPowerKL, reflectPowerLWS, td->KDT, hitArray, RxPos, gainRx/(4.0*SIPC_pi), nShadowRays, LRays, RRays, shadowRays, ranges, &(rnp[nAzBeam*nElBeam*nbounce]));
                    
                    free(LRays);
                    free(RRays);
                }
                
                free(shadowRays);
                free(shadowHits);
                free(ranges);
                
                memcpy(rayArray, reflectedRays, sizeof(Ray)*reflectCount);
                free(reflectedRays);

            }
            // Reset rayArrays and nRays for next time round loop
            //
            nRays = reflectCount ;
            
            nbounce++;
            
            free(hitArray);
        }
        
        int cnt=0;
        for (int i=0; i<nAzBeam*nElBeam*MAXBOUNCES; i++){
            if (rnp[i].power != 0 && rnp[i].range !=0) cnt++ ;
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
        rnpData_t * rnpData = (rnpData_t *)malloc(sizeof(rnpData_t)*nrnpItems) ;
        if(rnpData == NULL){
            printf("Error : malloc fail for rnpData on device %d. request data for %d items\n",tid,nrnpItems) ;
            exit(-15);
        }
        
        VECT_MINUS(TxPos, aimdir) ;
        derampRange = VECT_MAG(aimdir);
        
        cnt = 0;
        double totpow = 0.0;
        for (int i=0; i<nAzBeam*nElBeam*MAXBOUNCES; i++){
            if (rnp[i].power != 0 && rnp[i].range !=0){
                rnpData[cnt].power = rnp[i].power ;
                rnpData[cnt].range = rnp[i].range ; 
                rnpData[cnt].rdiff = rnp[i].range - derampRange;
                totpow += rnpData[cnt].power;
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
            
            double phse  = (-4.0 * SIPC_pi * rnpData[i].rdiff * td->oneOverLambda) ;
            targ.r = rnpData[i].power*cos(phse) ;
            targ.i = rnpData[i].power*sin(phse) ;
            
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
            
            CMPLX_MULT(pulseLine.data.cmpl_f[x], pcorr, tmp);
            pulseLine.data.cmpl_f[x] = tmp;
        }
        
        im_circshift(&pulseLine, -(pulseLine.nx/2), 0, &status);
        im_fftw(&pulseLine, FFT_X_ONLY+FWD+SCALE_N, &status);
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

