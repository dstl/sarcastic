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
#include "boxMullerRandom.h"
#include "ranf.h"

#define FASTPULSEBUILD
#define NPOINTS (32)
#define OVERSAMP (512)
#define NOINTERSECTION -1
//#define TOTALRCSINPULSE

void packSinc(SPCmplxD point, SPCmplx *outData, double rdiff, double sampleSpacing, long long nxInData, double * ikernel);
static void ham1dx(double * data, int nx);
void sinc_kernel(int oversample_factor, int num_of_points, double resolution, double sampleSpacing, double *ikernel);
void generate_gaussian_ray_distribution(double xStdev,double yStdev, SPVector txPos, SPVector GRP, int nx,int ny, Ray *rayArray);
void buildRays(Ray **rayArray, int *nRays, int nAzRays, int nElRays, int nTriangles, Triangle *triangles, SPVector TxPos, double PowPerRay, AABB SceneBoundingBox,
               cl_context context, cl_command_queue commandQ, cl_kernel  randRaysKL, size_t     randRaysLWS[2], SPVector **rayAimPoints );

typedef struct HitPoint {
    SPVector hit;       // Location of hitpoint in x,y,z
    int tri;            // index of triangle that this hit is on
} HitPoint ;

    
void * devPulseBlock ( void * threadArg ) {
    
    struct threadData *td;
    td = (struct threadData *) threadArg;
    
    cl_device_id devId ;
    cl_context context ;
    cl_command_queue commandQ ;

    int nAzBeam, nElBeam, bounceToShow, nrnpItems, nx, nShadowRays, tid ;
    int nbounce, nxRay, nyRay, nRays, reflectCount, iray, nShadows, interrogate ;
    int hrs,min,sec;
    int reportN = 500 ;

    double gainRx, PowPerRay, derampRange, derampPhase, pCentDone, sexToGo ;
    double rangeLabel, phasecorr,resolution,sampSpacing, *ikernel, bandwidth ;
    double *ranges, interogRad, intMinR, intMaxR, k ;

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

    gainRx           = td->gainRx ;
    PowPerRay        = td->PowPerRay ;
    bounceToShow     = td->bounceToShow ;
    nAzBeam          = td->nAzBeam ;
    nElBeam          = td->nElBeam ;
    SceneBoundingBox = td->SceneBoundingBox ;
    bandwidth        = td->chirpRate * td->pulseDuration ;
    tid              = td->devIndex ;
    devId            = td->platform.device_ids[tid] ;
    interrogate      = td->interrogate ;
    interogPt        = td->interogPt ;
    interogRad       = td->interogRad ;
    interogFP        = *(td->interogFP) ;
    k                = 2 * SIPC_pi / (SIPC_c / td->freq_centre) ;
    
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
    cl_program randRaysPG,     stackTraversePG,  reflectPG,  buildShadowsPG,  POFieldPG ;
    cl_kernel  randRaysKL,     stackTraverseKL,  reflectKL,  buildShadowsKL,  POFieldKL ;
    size_t     randRaysLWS[2], stackTraverseLWS, reflectLWS, buildShadowsLWS, POFieldLWS ;

    static char *randRaysCode      = OCLKERNELSPATH"/randomRays.cl" ;
    static char *stackTraverseCode = OCLKERNELSPATH"/stacklessTraverse.cl" ;
    static char *reflectCode       = OCLKERNELSPATH"/reflectRays.cl" ;
    static char *buildShadowsCode  = OCLKERNELSPATH"/buildShadowRays.cl" ;
    static char *reflectPowCode    = OCLKERNELSPATH"/POField.cl" ;
    
    CL_CHECK(buildKernel(context, randRaysCode,      "randomRays",        devId, 2, &randRaysPG,      &randRaysKL,      randRaysLWS));
    CL_CHECK(buildKernel(context, stackTraverseCode, "stacklessTraverse", devId, 1, &stackTraversePG, &stackTraverseKL, &stackTraverseLWS));
    CL_CHECK(buildKernel(context, reflectCode,       "reflect",           devId, 1, &reflectPG,       &reflectKL,       &reflectLWS));
    CL_CHECK(buildKernel(context, buildShadowsCode,  "buildShadowRays",   devId, 1, &buildShadowsPG,  &buildShadowsKL,  &buildShadowsLWS));
    CL_CHECK(buildKernel(context, reflectPowCode,    "POField",           devId, 1, &POFieldPG,       &POFieldKL,       &POFieldLWS));
  
    // Allocate memory for items that are needed in all kernels
    //
    int                nTriangles       = td->KDT.nTriangles;         // number of triangles in array 'triangles'
    ATS *              accelTriangles   = td->KDT.accelTriangles;     // Array of triangles of size nTriangles
    int                nTreeNodes       = td->KDT.nTreeNodes;         // number of nodes in KdTree
    KdData *           KdTree           = td->KDT.KdTree;             // SAH - KdTree to optimise ray traversal through volume
    int                triListDataSize  = td->KDT.triListDataSize;    // size of trianglelist data
    int *              triangleListData = td->KDT.triangleListData;   // array of triangle indices into Triangles
    int                nLeaves          = td->KDT.nLeaves;            // number of leaves (nodes with triangles) in KdTree
    int *              triangleListPtrs = td->KDT.triangleListPtrs;   // array of pointers into triangleListData for each KdTree node

    cl_mem dTriangles,dKdTree, dtriListData,dtriListPtrs;
    dTriangles   = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(ATS)*nTriangles, NULL, &_err));
    dKdTree      = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(KdData)*nTreeNodes,   NULL, &_err));
    dtriListData = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  triListDataSize,             NULL, &_err));
    dtriListPtrs = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(int)*nLeaves,         NULL, &_err));

    CL_CHECK(clEnqueueWriteBuffer(commandQ, dTriangles,   CL_TRUE, 0, sizeof(ATS)*nTriangles, accelTriangles,        0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(commandQ, dKdTree,      CL_TRUE, 0, sizeof(KdData)*nTreeNodes,   KdTree,           0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(commandQ, dtriListData, CL_TRUE, 0, triListDataSize,             triangleListData, 0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(commandQ, dtriListPtrs, CL_TRUE, 0, sizeof(int)*nLeaves,         triangleListPtrs, 0, NULL, NULL));
    
    if(td->phd != NULL){
        // build a sinc kernel
        //
        nx = (int)td->phd->nx ;
        resolution  = SIPC_c / (2.0 * bandwidth ) ;
        sampSpacing = SIPC_c / (2.0 * (nx/(td->pulseDuration * td->ADRate)) * bandwidth) ;
        ikernel     = sp_calloc((OVERSAMP * NPOINTS + 1), sizeof(double)) ;
        
        sinc_kernel(OVERSAMP, NPOINTS, resolution, sampSpacing, ikernel);
        
        printf("...Done\n");
    }
    
    // **** loop  start here
    //
    
    startTimer(&threadTimer, &status) ;
    SPVector *rayAimPoints = NULL ;

    for (int pulse=0; pulse<td->nPulses; pulse++){
        int pulseIndex = (tid * td->nPulses) + pulse ;
#ifndef TOTALRCSINPULSE
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
#endif //TOTALRCSINPULSE
        
        // Set correct parameters for beam to ray trace
        //
        TxPos = td->TxPositions[pulseIndex] ;
        RxPos = td->RxPositions[pulseIndex] ;
        
        // Generate a distribution of nAzbeam x nElbeam rays that originate from the TxPosition aiming at the origin. Use beamMax as the std deviation
        // for the distribution
        //
        
        nbounce = 0;
        nxRay   = nAzBeam ;
        nyRay   = nElBeam ;
        nRays   = nxRay*nyRay;
        buildRays(&rayArray, &nRays, nAzBeam, nElBeam, nTriangles, td->triangles, TxPos, PowPerRay, td->SceneBoundingBox, context, commandQ, randRaysKL, randRaysLWS, &rayAimPoints);
        
        // Set up deramp range for this pulse
        //
        VECT_MINUS(TxPos, aimdir) ;
        derampRange = VECT_MAG(aimdir);
        derampPhase = -4.0 * SIPC_pi * derampRange * td->oneOverLambda ;
        
        // Use Calloc for rnp as we will be testing for zeroes later on
        //
        rnp = (rangeAndPower *)sp_calloc(nRays*MAXBOUNCES, sizeof(rangeAndPower));
        
        while ( nbounce < MAXBOUNCES &&  nRays != 0){
            
            // Malloc space for hits for this bounce
            //
            hitArray = (Hit *)sp_malloc(nRays * sizeof(Hit));

            // Cast the rays in rayArray through the KdTree using a stackless traversal technique
            // return the hit locations in hitArray
            //
            oclKdTreeHits(context, commandQ, stackTraverseKL, nRays, stackTraverseLWS, dTriangles, dKdTree, dtriListData, dtriListPtrs, SceneBoundingBox, rayArray, hitArray);
            
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
            newRays       = (Ray *)sp_malloc(reflectCount * sizeof(Ray)) ;
            newHits       = (Hit *)sp_malloc(reflectCount * sizeof(Hit)) ;
            reflectedRays = (Ray *)sp_malloc(reflectCount * sizeof(Ray)) ;
            
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
            oclReflect(context, commandQ, reflectKL, dTriangles, nRays, reflectLWS, rayArray, hitArray, reflectedRays);
            
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
            shadowRays =    (Ray *)sp_malloc(reflectCount * sizeof(Ray)) ;
            shadowHits =    (Hit *)sp_malloc(reflectCount * sizeof(Hit)) ;
            ranges     = (double *)sp_malloc(reflectCount * sizeof(double)) ;
            
            // Build Shadowrays
            //
            oclBuildShadowRays(context, commandQ, buildShadowsKL, buildShadowsLWS, reflectCount, RxPos, reflectedRays, shadowRays, ranges);
            
            // Work out which rays have a path back to receiver using stackless traverse kernel
            //
            oclKdTreeHits(context, commandQ, stackTraverseKL, nRays, stackTraverseLWS, dTriangles , dKdTree, dtriListData, dtriListPtrs, SceneBoundingBox, shadowRays, shadowHits);
            
            // Shrink the shadowRays to only include those that made it back to the sensor
            // in order to calculate power at sensor we also need the Illumination or LRays
            //
            nShadows = 0 ;
            iray     = 0 ;
            for(int i= 0 ; i<nRays; i++) if ( shadowHits[i].trinum == NOINTERSECTION ) nShadows++ ;
            
            if( nShadows != 0){
                
                newRays       = (Ray *)sp_malloc(nShadows * sizeof(Ray)) ;
                newHits       = (Hit *)sp_malloc(nShadows * sizeof(Hit)) ;
                LRays         = (Ray *)sp_malloc(nShadows * sizeof(Ray)) ;
                RRays         = (Ray *)sp_malloc(nShadows * sizeof(Ray)) ;
                
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
                oclPOField(context, commandQ, POFieldKL, POFieldLWS, td->triangles, nTriangles, hitArray, nShadowRays, LRays, shadowRays, RxPos, k, ranges, gainRx, nbounce+1, &(rnp[nRays*nbounce])) ;
                
                
                
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
                    
                    for ( int i=0; i<nRays; i++){
                        
                        irp = rnp[(nRays*nbounce)+i] ;
                        
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
        
        if(td->phd != NULL){

            int cnt=0;
            for (int i=0; i<nRays*MAXBOUNCES; i++){
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
            rnpData_t * rnpData = (rnpData_t *)sp_malloc(nrnpItems * sizeof(rnpData_t));
            cnt = 0;
            for (int i=0; i<nRays*MAXBOUNCES; i++){
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
            printf("%d, %e\n",pulse,cmplx_mag);
//            printf("Total RCS for pulse %d is %f m^2 (%f dB m^2)\n",pulse,cmplx_mag,10*log(cmplx_mag));
//            printf("For comparison: \n");
//            printf("    1m^2 flat plate : %f (%f dB m^2)\n",4*SIPC_pi*td->oneOverLambda*td->oneOverLambda,10*log(4*SIPC_pi*td->oneOverLambda*td->oneOverLambda));
//            printf("    1m dihedral     : %f (%f dB m^2)\n",8*SIPC_pi*td->oneOverLambda*td->oneOverLambda,10*log(8*SIPC_pi*td->oneOverLambda*td->oneOverLambda));
//            printf("    1m trihedral    : %f (%f dB m^2)\n",12*SIPC_pi*td->oneOverLambda*td->oneOverLambda,10*log(12*SIPC_pi*td->oneOverLambda*td->oneOverLambda));
#endif // TOTALRCSINPULSE

            // perform phase correction to account for deramped jitter in receiver timing
            //
            phasecorr = (((td->Fx0s[pulseIndex] - td->freq_centre) / td->FxSteps[pulseIndex])) * 2.0 * M_PI / pulseLine.nx;
            
            for(int x = 0; x < pulseLine.nx; x++) {
                pcorr.r =  cos(phasecorr * (x - pulseLine.nx/2)) / (td->amp_sf0[pulseIndex]);
                pcorr.i =  sin(phasecorr * (x - pulseLine.nx/2)) / (td->amp_sf0[pulseIndex]);
                
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
            
        } // end of pulse loop
        
        free(rnp) ;
    }
    
    if(td->phd != NULL){
        free(ikernel);
    }
    
    free( rayAimPoints );
    // Clear down OpenCL allocations
    //
    clReleaseKernel(randRaysKL);
    clReleaseKernel(stackTraverseKL);
    clReleaseKernel(reflectKL);
    clReleaseKernel(buildShadowsKL);
    clReleaseKernel(POFieldKL);
    clReleaseProgram(randRaysPG) ;
    clReleaseProgram(stackTraversePG) ;
    clReleaseProgram(reflectPG) ;
    clReleaseProgram(buildShadowsPG) ;
    clReleaseProgram(POFieldPG) ;
    clReleaseCommandQueue(commandQ);
    clReleaseContext(context);
    
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

void buildRays(Ray **rayArray, int *nRays, int nAzRays, int nElRays, int nTriangles, Triangle *triangles, SPVector TxPos,
               double PowPerRay, AABB SceneBoundingBox,
               cl_context context, cl_command_queue commandQ, cl_kernel  randRaysKL, size_t randRaysLWS[2],
               SPVector **rayAimPoints
               ){
    
    int METHOD  = 4;
    
    // 1 - each ray aimed at triangle centre
    // 2 - random rays on each call across scene
    // 3 - random rays created first time but the same hitpoints used for each subsequent call
    // 4 - like 2 (random rays on each call across the scene) but rays are parallel from Tx
    
    if(METHOD == 1){
        *nRays = nTriangles;
        *rayArray = (Ray *)sp_malloc(*nRays * sizeof(Ray));

        for(int i=0;i<nTriangles; i++){
            Triangle t = triangles[i];
            SPVector mean;

            for(int j=0; j<3; j++)mean.cell[j] = (t.AA.cell[j]+t.BB.cell[j]+t.CC.cell[j]) / 3.0 ;
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
        double centreRange = VECT_MAG(TxPos);
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
        
        double maxBeamUsedAz = (maxAz - minAz) / centreRange ;
        double maxBeamUsedEl = (maxEl - minEl) / centreRange ;
        SPVector aimpoint;
        VECT_CREATE(SceneBoundingBox.AA.x+((SceneBoundingBox.BB.x - SceneBoundingBox.AA.x)/2),
                    SceneBoundingBox.AA.y+((SceneBoundingBox.BB.y - SceneBoundingBox.AA.y)/2),
                    SceneBoundingBox.AA.z+((SceneBoundingBox.BB.z - SceneBoundingBox.AA.z)/2), aimpoint);
        
        *rayArray = (Ray *)sp_malloc(*nRays * sizeof(Ray));
        
        oclRandomRays(context,commandQ,randRaysKL,nAzBeam,nElBeam,randRaysLWS,maxBeamUsedAz/2,maxBeamUsedEl/2,TxPos, aimpoint, PowPerRay, *rayArray);
        
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