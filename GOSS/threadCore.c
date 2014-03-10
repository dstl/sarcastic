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

static char *STkernelCodePath = "/Users/darren/Development/GOSS/GOSS/stacklessTraverse.cl" ;
static char *RRkernelCodePath = "/Users/darren/Development/GOSS/GOSS/randomRays.cl" ;

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
    cl_program STprogram, RRprogram ;
    cl_kernel STkernel, RRkernel ;
    cl_int err ;
    rangeAndPower * rnp ;
    SPVector aimdir ;
    
    int nAzBeam             = td->nAzBeam ;
    int nElBeam             = td->nElBeam ;
    SPVector RxPos;
    SPVector TxPos;
    double gainRx           = td->gainRx ;
    double TxPowPerRay      = td->TxPowPerRay ;
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
    fftwf_init_threads();
    printf("Using FAST pulse building routines...\n");

    int tid ;
    tid   = td->devIndex ;
    devId = td->platform.device_ids[tid] ;
    context  = CL_CHECK_ERR(clCreateContext(td->platform.props, 1, &devId, NULL, NULL, &_err));
    commandQ = CL_CHECK_ERR(clCreateCommandQueue(context, devId, 0, &_err));

    // build and compile the Ray Tracing kernel
    //
    STprogram = CL_CHECK_ERR(clCreateProgramWithSource(context, 1, (const char **) &STkernelCodePath, NULL, &_err));
    
    char STcompilerOptions[255];
    sprintf(STcompilerOptions, "-Werror -I/Users/darren/Development/GOSS/GOSS/") ;
    
    err = clBuildProgram(STprogram, 0, NULL, STcompilerOptions, NULL, NULL);
    if (err != CL_SUCCESS){
        size_t len;
        char buffer[32768];
        printf("[thread:%d], Error: Failed to build program executable!\n",tid);
        clGetProgramBuildInfo(STprogram, devId, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        printf("err: %d. Buffer:\n",err);
        printf("%s\n", buffer);
        exit(-3);
    }
    
    STkernel   = CL_CHECK_ERR(clCreateKernel(STprogram, "stackLessTraverse", &_err));
    
    // build and compile the Random Ray Generation kernel
    //
    RRprogram = CL_CHECK_ERR(clCreateProgramWithSource(context, 1, (const char **) &RRkernelCodePath, NULL, &_err));
    
    char RRcompilerOptions[255];
    sprintf(RRcompilerOptions, "-Werror -I/Users/darren/Development/GOSS/GOSS/") ;
    err = clBuildProgram(RRprogram, 0, NULL, RRcompilerOptions, NULL, NULL);
    if (err != CL_SUCCESS){
        size_t len;
        char buffer[32768];
        printf("[thread:%d], Error: Failed to build program executable!\n",tid);
        clGetProgramBuildInfo(RRprogram, devId, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        printf("err: %d. Buffer:\n",err);
        printf("%s\n", buffer);
        exit(-3);
    }
    
    RRkernel   = CL_CHECK_ERR(clCreateKernel(RRprogram, "randomRays", &_err));
    
    
    // Calculate global and local work sizes
    //
    int tx,ty ;
    size_t STlocalWorkSize [2] ;
    size_t STglobalWorkSize[2] ;
    STglobalWorkSize[0] = nAzBeam ;
    STglobalWorkSize[1] = nElBeam ;
    best2DWorkSize(STkernel, devId, STglobalWorkSize[0], STglobalWorkSize[1], &tx, &ty, &td->status) ;
    STlocalWorkSize[0]  = tx ;
    STlocalWorkSize[1]  = ty ;
    
    size_t RRlocalWorkSize [2] ;
    size_t RRglobalWorkSize[2] ;
    RRglobalWorkSize[0] = nAzBeam ;
    RRglobalWorkSize[1] = nElBeam ;
    best2DWorkSize(RRkernel, devId, RRglobalWorkSize[0], RRglobalWorkSize[1], &tx, &ty, &td->status) ;
    RRlocalWorkSize[0]  = tx ;
    RRlocalWorkSize[1]  = ty ;
    
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
    nx = (int)td->phd->nx ;

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
    
    
    // generate a random gaussian distribution of rays and pass these to the openCL kernel
    //
    SPVector origin;
    VECT_CREATE(0, 0, 0, origin);
    Ray *rayArray;
    Hit *hitArray;
    Ray *shadowRays;
    double *distnces ;
    double *tempDists ;
    
    rnp      = (rangeAndPower *)malloc(sizeof(rangeAndPower)*nAzBeam*nElBeam*MAXBOUNCES);
    
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
        hitArray = (Hit *)malloc(sizeof(Hit)*nAzBeam*nElBeam);
        distnces = (double *)malloc(sizeof(double)*nAzBeam*nElBeam);
        
        oclRandomRays(context,commandQ,RRkernel,RRglobalWorkSize,RRlocalWorkSize,td->beamMaxAz,td->beamMaxEl,TxPos,origin, rayArray);
        
        int nbounce = 0;
        int nxRay   = nAzBeam ;
        int nyRay   = nElBeam ;
        int nRays   = nxRay*nyRay;
        int hitcount ;
        int iray;
        SPVector hp,dir;

        while ( nbounce < MAXBOUNCES &&  nRays != 0){
            
            
            // Calculate work size of OpenCL
            //
            
            
            oclKdTreeHits(context, commandQ, STkernel, STglobalWorkSize, STlocalWorkSize, td->KDT, SceneBoundingBox, rayArray, hitArray) ;
            
            // Sort out hits
            //
            hitcount = 0;
            iray = 0;
            // How many hits occured on this ray cast
            //
            for(int i=0; i<nxRay*nyRay; i++) if ( hitArray[i].trinum != NOINTERSECTION ) hitcount++ ;
            
            // Malloc space for shadowRays (arrays form a hit going back to receiver
            // Malloc space for tempDists (array to temporarily store distances for forward casted rays
            //
            shadowRays = (Ray *)malloc(sizeof(Ray)*hitcount) ;
            tempDists  = (double *)malloc(sizeof(double)*hitcount) ;
            
            // Build Shadowrays and update distances
            //
            for (int i=0; i<nxRay*nyRay; i++) {
                if ( hitArray[i].trinum != NOINTERSECTION ){
                    VECT_SCMULT(rayArray[i].dir, hitArray[i].dist,hp);
                    VECT_ADD(rayArray[i].org, hp, hp);
                    shadowRays[iray].org   = hp ;
                    VECT_SUB(RxPos, hp, dir);
                    VECT_NORM(dir, shadowRays[iray].dir);
                    tempDists[iray] = distnces[i] + hitArray[i].dist ;
                    iray++ ;
                }
            }
            
            // Build forward scattering rays ready for next turn round the loop
            //
            
            oclReflect(rayArray, hitArray, newRayArray, &newRayArraySize);
            free(rayArray) ;
            rayArray = newRayarray;

            
            // Work out which rays have a path back to receiver and calculate power vs range for the origin of each ray
            //
            
            
            
            
            
            nbounce++;
            free(shadowRays);
            free(shadDists);
        }
        
        free(rayArray);
        free(hitArray);
        free(distnces);
        
        // Now use OpenCl to cast the rays through the KdTree and return an array of intersections with values for teh range
        // of each intersection and the power at that intersection (stored in array 'rnp'
        //
        /*oclRayTrace(context,
                    commandQ,
                    RTkernel,
                    RTglobalWorkSize,
                    RTlocalWorkSize,
                    td->nTriangles,
                    td->Triangles,
                    td->nTextures,
                    td->Textures,
                    td->nTreeNodes,
                    td->KdTree,
                    td->triListDataSize,
                    td->triListData,
                    td->nLeaves,
                    td->triPtrs,
                    RxPos,
                    gainRx,
                    TxPowPerRay,
                    SceneBoundingBox,
                    bounceToShow,
                    pulseIndex,
                    rayArray,
                    rnp);*/
        
        // Some thoughts on how it should be
        //
        // generate random rays
        //
        // While bounces < max bounces and nRays > 0
        // order random rays in some way
        // stackless traversal for all rays in parallel - generates hit list
        // each hit contains inbound ray and triangle
        // In thread 1:
        // Generate array of rays back to receiver
        // sort rays in some way
        // Ray cast arrays back to receiver to check for occlusion. If a ray is occluded
        // then remove it. If it makes it back then calculate power from power at
        // reflection point and range back to receiver
        // (power is calculated using a nominal area for each ray, the phong properties
        // of the scattering triangle and the phong equation)
        // Add range and power to a linked list
        //
        // In thread 2:
        // calculate all forward ray bounces
        // store total range to ray origin and power of ray
        // (power calculated from reflection)
        // end while
        // Build output pulse from RVP's
        
        int cnt=0;
        for (int i=0; i<nAzBeam*nElBeam*MAXBOUNCES; i++){
            if (rnp[i].power != 0 && rnp[i].range !=0) cnt++ ;
        }
        if(cnt > 0){
            nrnpItems = cnt;
        }else{
            printf("ERROR: No intersections on device %d\n",tid);
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
        
        VECT_MINUS(td->RxPositions[pulseIndex], aimdir) ;
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
//        printf("Total power : %e, # of intersecting rays : %d, powerperray: %f db (%f)\n",totpow,nrnpItems,10*log(TxPowPerRay), TxPowPerRay);
        
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
        
    } // end of pulse loop
    free(rayArray);
    free(rnp);
    free(ikernel);
    
    // Clear down OpenCL allocations
    //
    clReleaseCommandQueue(commandQ);
    clReleaseContext(context);
    clReleaseKernel(RTkernel);
    clReleaseProgram(RTprogram);
    
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

