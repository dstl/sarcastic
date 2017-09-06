/***************************************************************************
 *
 *       Module:    threadCore.hpp
 *      Program:    bircs
 *   Created by:    Darren Muff on 12/06/2017.
 *                  Copyright (c) 2017 [Dstl]. All rights reserved.
 *
 *   Description:
 *     Core thread routine for BIRCS
 *
 *
 *   CLASSIFICATION        :  OFFICIAL
 *   Date of CLASSN        :  12/06/2017
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
#include <fftw3.h>
extern "C" {
#include "boxMullerRandom.h"
#include "RCS.h"
}
#include "ranf.h"
#include "threadCore.hpp"
#include <sys/time.h>
#if defined (__APPLE__) || defined(MACOSX)
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif
#include "buildKernel.hpp"
#include "rayTrace.hpp"
#include "cpuPOField.hpp"
#include <SILib2/SILib2.h>

#define NOINTERSECTION -1


void * devPulseBlock ( void * threadArg ) {
    
    struct bircsThreadData *td;
    td = (struct bircsThreadData *) threadArg;
    
    cl_device_id devId ;
    cl_context context ;
    cl_command_queue commandQ ;
    cl_mem dTriangles;
    
    int nAzBeam, nElBeam, bounceToShow, nrnpItems, nShadowRays, tid ;
    int nbounce, nxRay, nyRay, nRays, reflectCount, iray, nShadows, interrogate ;
    int maxRaysPerBounce, pol ;

    
    double gainRx, PowPerRay, derampRange ;
    double *ranges, *newranges, interogRad, k ;
    
    AABB SceneBoundingBox ;
    rangeAndPower * rnp ;
    Ray *rayArray, *newRays, *reflectedRays, *shadowRays, *LRays, *RRays;
    Hit *hitArray, *newHits, *shadowHits;
    
    SPVector aimdir, RxPos, TxPos, origin, interogPt;
    SPVector maxRCS = {-9e99,0.0,0.0};
    SPStatus status;
    
    FILE *interogFP ;
    
    Timer threadTimer ;
    kdTree::KdData *node;
    ATS *accelTriangles = NULL;
    int Ropes[6] ;
    AABB sceneAABB ;
    
    kdTree::KdData * tree = NULL;
    int treeSize;
    int nTriangles;
    
    gainRx           = td->gainRx ;
    PowPerRay        = td->PowPerRay ;
    bounceToShow     = td->bounceToShow ;
    nAzBeam          = td->nAzBeam ;
    nElBeam          = td->nElBeam ;
    SceneBoundingBox = td->SceneBoundingBox ;
    tid              = td->devIndex ;
    devId            = td->platform.device_ids[tid] ;
    interrogate      = td->interrogate ;
    interogPt        = td->interogPt ;
    interogRad       = td->interogRad ;
    interogFP        = *(td->interogFP) ;
    k                = 2 * SIPC_pi / (SIPC_c / td->freq_centre) ;
    pol              = td->polarisation ;
    
    VECT_CREATE(0, 0, 0, origin);
    
    fftwf_init_threads();
    im_init_status(status, 0) ;
    
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
    cl_kernel  randRaysKL,     reflectKL,  buildShadowsKL,  POFieldKL ;
    size_t     randRaysLWS[2], /*stackTraverseLWS,*/ reflectLWS, buildShadowsLWS, POFieldLWS ;
    
    static const char *randRaysCode      = OCLKERNELSPATH"/randomRays.cl" ;
    static const char *reflectCode       = OCLKERNELSPATH"/reflectRays.cl" ;
    static const char *buildShadowsCode  = OCLKERNELSPATH"/buildShadowRays.cl" ;
    static const char *reflectPowCode    = OCLKERNELSPATH"/POField2.cl" ;
    
    CL_CHECK(buildKernel(context, randRaysCode,      "randomRays",        devId, 2, &randRaysPG,      &randRaysKL,      randRaysLWS));
    CL_CHECK(buildKernel(context, reflectCode,       "reflect",           devId, 1, &reflectPG,       &reflectKL,       &reflectLWS));
    CL_CHECK(buildKernel(context, buildShadowsCode,  "buildShadowRays",   devId, 1, &buildShadowsPG,  &buildShadowsKL,  &buildShadowsLWS));
    CL_CHECK(buildKernel(context, reflectPowCode,    "POField",           devId, 1, &POFieldPG,       &POFieldKL,       &POFieldLWS));
    
    // Allocate memory for items that are needed in all kernels
    //
    //    int nTriangles       = (int)td->sceneMesh->triangles.size() ; // number of triangles in array 'triangles'
    //    ATS *accelTriangles  = td->accelTriangles;     // Array of triangles of size nTriangles
    //    int nTreeNodes = td->nNodes;         // number of nodes in KdTree
    //    kdTree::KdData *KdTree = td->tree ;             // SAH - KdTree to optimise ray traversal through volume
    
    //    cl_mem dTriangles,dKdTree, dtriListData,dtriListPtrs;
    //    dTriangles   = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(ATS)*nTriangles, NULL, &_err));
    //    dKdTree      = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(kdTree::KdData)*nTreeNodes,   NULL, &_err));
    //
    //    CL_CHECK(clEnqueueWriteBuffer(commandQ, dTriangles,   CL_TRUE, 0, sizeof(ATS)*nTriangles, accelTriangles, 0, NULL, NULL));
    //    CL_CHECK(clEnqueueWriteBuffer(commandQ, dKdTree,      CL_TRUE, 0, sizeof(kdTree::KdData)*nTreeNodes,KdTree,  0, NULL, NULL));
    

    // **** loop  start here
    //
    
    startTimer(&threadTimer, &status) ;
    SPVector *rayAimPoints = NULL ;
    
    TriangleMesh newMesh ;
    
    newMesh = *(td->sceneMesh) ;
    
    // Build kdTree if required
    //
    
    kdTree::buildTree(newMesh, &tree, &treeSize, (kdTree::TREEOUTPUT)(kdTree::OUTPUTSUMM)) ;
    accelerateTriangles(&newMesh,&accelTriangles) ;
    // Initialise the tree and build ropes and boxes to increase efficiency when traversing
    //
    node = &(tree[0]) ;
    
    // build Ropes and Boxes
    //
    sceneAABB = tree[0].brch.aabb ;
    for(int i=0; i<6; i++) Ropes[i] = NILROPE;
    BuildRopesAndBoxes(node, Ropes, sceneAABB, tree);
    
    // Allocate memory for items that are needed in all kernels
    //
    nTriangles   = (int)newMesh.triangles.size() ; // number of triangles in array 'triangles'
    dTriangles   = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(ATS)*nTriangles, NULL, &_err));
    CL_CHECK(clEnqueueWriteBuffer(commandQ, dTriangles,   CL_TRUE, 0, sizeof(ATS)*nTriangles, accelTriangles, 0, NULL, NULL));

    
    for (int pulse=0; pulse<td->nPulses; pulse++){
        int pulseIndex = (tid * td->nPulses) + pulse ;
        
        // Set correct parameters for beam to ray trace
        //
        TxPos = td->TxPositions[pulseIndex] ;
        RxPos = td->RxPositions[pulseIndex] ;
        double phi = atan2(RxPos.y,RxPos.x);
        double theta = atan2(sqrt(RxPos.x*RxPos.x+RxPos.y*RxPos.y),RxPos.z);
        
        
        // Generate a distribution of nAzbeam x nElbeam rays that originate from the TxPosition aiming at the origin. Use beamMax as the std deviation
        // for the distribution
        //
        nbounce = 0;
        nxRay   = nAzBeam ;
        nyRay   = nElBeam ;
        nRays   = nxRay*nyRay;
        buildRays(&rayArray, &nRays, nAzBeam, nElBeam, &newMesh, TxPos, PowPerRay, td->SceneBoundingBox, &rayAimPoints, TRIANGLECENTRE, pol);
        maxRaysPerBounce = nRays;  // Use this for memory as nxRay/nyRay may be incorrect if buildRays set to triangle centres
        
        // Set up deramp range for this pulse
        //
        VECT_MINUS(TxPos, aimdir) ;
        derampRange = VECT_MAG(aimdir);
        
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
            //            oclReflect(context, commandQ, reflectKL, dTriangles, nRays, reflectLWS, rayArray, hitArray, reflectedRays);
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
            //            oclBuildShadowRays(context, commandQ, buildShadowsKL, buildShadowsLWS, nRays, RxPos, reflectedRays, shadowRays, ranges);
            
            // Work out which rays have a path back to receiver using stackless traverse kernel
            //
            //            oclKdTreeHits(context, commandQ, stackTraverseKL, nRays, stackTraverseLWS, dTriangles , dKdTree, dtriListData, dtriListPtrs, SceneBoundingBox, shadowRays, shadowHits);
            shootRay(tree, accelTriangles, nRays, shadowRays, shadowHits) ;
            
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
                //  oclPOField(context, commandQ, POFieldKL, POFieldLWS, td->triangles, nTriangles, hitArray, nShadowRays, LRays, shadowRays, RxPos, k, ranges, gainRx, nbounce+1, &(rnp[nRays*nbounce])) ;
                //                cpuPOField(td->triangles, nTriangles, hitArray, nShadowRays, LRays, shadowRays, RxPos, k, ranges, gainRx, nbounce+1, &(rnp[nxRay*nyRay*nbounce])) ;
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
            
            cmplx_mag = 10*log10(RCS(PowPerRay, CMPLX_MAG(targtot), derampRange, derampRange));
            free(rnpData);
        }
        
        if(cmplx_mag <0 ){
            printf("%f, %f, %f\n",0.0,phi,theta);
        }else{
            printf("%f, %f, %f\n",cmplx_mag,phi,theta);
        }
        
        if (maxRCS.R < cmplx_mag){
            maxRCS.R = cmplx_mag;
            maxRCS.phi = phi;
            maxRCS.theta = theta;
        }

        free(rnp) ;
    }
    
    double lambda = (SIPC_c / td->freq_centre);
    double lambda_squared = lambda * lambda ;
    printf("Maximum RCS was %f dB m^2 (%f m^2) at az:%fdeg, incidence: %fdeg\n",maxRCS.R,pow(10.0,maxRCS.R/10.0),RAD2DEG(maxRCS.phi), RAD2DEG(maxRCS.theta) );
    printf("    1m^2 flat plate               : %f (%f dB m^2)\n",4*SIPC_pi / lambda_squared,10*log10(4*SIPC_pi / lambda_squared)); // 8620.677
    printf("    1m dihedral                   : %f (%f dB m^2)\n",8*SIPC_pi / lambda_squared,10*log10(8*SIPC_pi / lambda_squared)); // 17241.354
    printf("    1m trihedral                  : %f (%f dB m^2)\n",12*SIPC_pi/ lambda_squared,10*log10(12*SIPC_pi/ lambda_squared)); // 25862.031
    printf("    1m dia Sphere                 : \n");
    printf("        Rayleigh Region r<<lambda : %f ( %f dB m^2)\n", 9*SIPC_pi*0.5*0.5*pow(k*0.5,4), 10*log10(9*SIPC_pi*0.5*0.5*pow(k*0.5,4)));
    printf("        Optical Region  r>>lambda : %f ( %f dB m^2)\n", SIPC_pi*0.5*0.5, 10*log10(SIPC_pi*0.5*0.5));
    printf("         ( r = 0.5 m, lambda = %5.3f m )\n",lambda);

    
    free( rayAimPoints );
    free(tree);
    delete accelTriangles ;
    clReleaseMemObject(dTriangles);
    // Clear down OpenCL allocations
    //
    clReleaseKernel(randRaysKL);
    //    clReleaseKernel(stackTraverseKL);
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

/*void buildRays(Ray **rayArray, int *nRays, int nAzRays, int nElRays, TriangleMesh *mesh, SPVector TxPos,
               double PowPerRay, AABB SceneBoundingBox,
               cl_context context, cl_command_queue commandQ, cl_kernel  randRaysKL, size_t randRaysLWS[2],
               SPVector **rayAimPoints)
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
    
}*/
