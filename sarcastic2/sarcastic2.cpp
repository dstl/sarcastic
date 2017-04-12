/***************************************************************************
 *
 *       Module:    sarcastic2.cpp
 *      Program:    SARCASTIC
 *   Created by:    Darren on 15/03/2017.
 *                  Copyright (c) 2017 Dstl. All rights reserved.
 *
 *   Description:
 *
 *
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  15/03/2017
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

#include <iostream>
#include "sarcastic2.hpp"
#include <SIlib2/SIlib2.h>
#include "tryReadFile.hpp"
#include "TriangleMesh.hpp"
#include "buildTree.hpp"
#include "buildRopesAndBoxes.hpp"
#include "accelerateTriangles.hpp"
extern "C" {
#include "TxPowerPerRay.h"
#include "ecef2SceneCoords.h"
#include "OpenCLUtils.h"
}
#include "getUserInput.hpp"
#include "threadCore.hpp"

int main(int argc, const char * argv[]) {
    
    int rc;

    
    SPStatus status ;
    im_init_status(status, 0) ;
    im_init_lib(&status, (char *)"sarcastic", argc, (char **)argv);
    CHECK_STATUS_NON_PTR(status);
    
    CPHDHeader hdr ;
    TriangleMesh baseMesh ;
    TriangleMesh moverMesh ;
    char *outCPHDFile ;
    int startPulse, nPulses, bounceToShow, nAzBeam, nElBeam, interrogate, pulseUndersampleFactor ;
    SPVector interogPt ;
    double interogRad ;
    FILE *interrogateFP ;
    
    getUserInput(&hdr, &baseMesh, &moverMesh, &outCPHDFile,
                     &startPulse, &nPulses, &bounceToShow, &nAzBeam, &nElBeam, &interrogate, &interogPt, &interogRad,
                 &interrogateFP, &pulseUndersampleFactor, &status) ;
    
    // Start timing after user input
    //
    Timer runTimer ;
    startTimer(&runTimer, &status) ;
    
    // Reduce the Cphd data we have to deal with to make things quicker.
    //
    CPHDHeader newhdr = hdr ;
    newhdr.num_azi = nPulses / pulseUndersampleFactor ;
    newhdr.pulses  = (CPHDPulse *)sp_calloc(newhdr.num_azi, sizeof(CPHDPulse));
    for(int p = 0; p<newhdr.num_azi; ++p){
        newhdr.pulses[p] = hdr.pulses[(p*pulseUndersampleFactor)+startPulse] ;
    }
    nPulses = newhdr.num_azi ;


    // print out info about the collection
    //
    printCPHDCollectionInfo(&hdr, &status) ;
    printf("Wavelength                  : %f m\n",SIPC_c/newhdr.freq_centre);
    printf("Slant range resolution      : %f m\n",SIPC_c/(2*newhdr.pulse_length*newhdr.chirp_gamma));
    
    // Rotate Rx and Tx Coords to be relative to scene centre
    //
    SPVector  *RxPos, *TxPos, SRP, zHat ;
    double    *Fx0s, *FxSteps, *amp_sf0, lambda ;
    
    RxPos = (SPVector *)sp_malloc(sizeof(SPVector)*nPulses) ;
    TxPos = (SPVector *)sp_malloc(sizeof(SPVector)*nPulses) ;
//    Fx0s    = (double *)sp_malloc(sizeof(double)*nPulses)   ;
//    FxSteps = (double *)sp_malloc(sizeof(double)*nPulses)   ;
//    amp_sf0 = (double *)sp_malloc(sizeof(double)*nPulses)   ;
   
    for (int p = 0; p < nPulses; p++){
        RxPos[p] = newhdr.pulses[p].sat_ps_rx ;
        TxPos[p] = newhdr.pulses[p].sat_ps_tx ;
        
//        // While we are looping through pulses also load up the Fx0 and Fx®StepSize values
//        //
//        Fx0s[p]    = newhdr.pulses[p].fx0 ;
//        FxSteps[p] = newhdr.pulses[p].fx_step_size ;
//        amp_sf0[p] = newhdr.pulses[p].amp_sf0 ;
    }
    
    SRP = newhdr.pulses[(nPulses/2)].srp ;
    
    ecef2SceneCoords(nPulses, RxPos, SRP);
    ecef2SceneCoords(nPulses, TxPos, SRP);
    ecef2SceneCoords(1,  &interogPt, SRP);
    
    SPVector tmp;
    for (int p = 0; p < nPulses; p++){
        tmp = newhdr.pulses[p].sat_ps_rx ;
        newhdr.pulses[p].sat_ps_rx = RxPos[p] ;
        RxPos[p]  = tmp;
        tmp = newhdr.pulses[p].sat_ps_tx ;
        newhdr.pulses[p].sat_ps_tx = TxPos[p] ;
        TxPos[p] = tmp;
    }
    
    
    double centreRange, maxEl, maxAz, minEl, minAz, dAz, dEl, maxBeamUsedAz, maxBeamUsedEl ;
    SPVector rVect, unitBeamAz, unitBeamEl  ;
    AABB SceneBoundingBox ;
    centreRange = VECT_MAG(TxPos[nPulses/2]);
    VECT_MINUS( TxPos[nPulses/2], rVect ) ;
    VECT_CREATE(0, 0, 1., zHat) ;
    VECT_CROSS(rVect, zHat, unitBeamAz);
    VECT_NORM(unitBeamAz, unitBeamAz) ;
    VECT_CROSS(unitBeamAz, rVect, unitBeamEl) ;
    VECT_NORM(unitBeamEl, unitBeamEl) ;
    
    maxEl = maxAz = minEl = minAz = 0.0 ;
    SPVector min,max; VECT_CREATE(9e10, 9e10, 9e10, min); VECT_CREATE(-9e10, -9e10, -9e10, max);
    for(int i=0; i<baseMesh.triangles.size(); ++i){
        for(int j=0; j<3; ++j){
            min.cell[j] = (baseMesh.AABBs[i].AA.cell[j] < min.cell[j]) ? baseMesh.AABBs[i].AA.cell[j] : min.cell[j] ;
            max.cell[j] = (baseMesh.AABBs[i].BB.cell[j] > max.cell[j]) ? baseMesh.AABBs[i].BB.cell[j] : max.cell[j] ;
        }
    }
    SceneBoundingBox.AA = min ; SceneBoundingBox.BB = max ;
    SPVector boxPts[8] ;
    
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
    maxBeamUsedAz = (maxAz - minAz) / centreRange ;
    maxBeamUsedEl = (maxEl - minEl) / centreRange ;

    collectionGeom cGeom;
    collectionGeometry(&hdr, nPulses/2, hdr.grp, &cGeom, &status);
    printf("RayTracing beam (Az,El)     : %3.2e deg x %3.2e deg (%3.2f x %3.2f metres (ground plane))\n",
           GeoConsts_RADTODEG*maxBeamUsedAz,GeoConsts_RADTODEG*maxBeamUsedEl,centreRange*maxBeamUsedAz,centreRange*maxBeamUsedEl/sin(cGeom.grazingRad));
    dAz = maxBeamUsedAz / nAzBeam;
    dEl = maxBeamUsedEl / nElBeam;
    lambda = SIPC_c / newhdr.freq_centre ;
    printf("Ray density                 : %3.1f x %3.1f [ %3.1f x %3.1f ground plane] rays per metre\n",
           1.0/(dAz * centreRange), 1.0/(dEl * centreRange),
           1.0/(dAz * centreRange), 1.0/(dEl * centreRange / sin(cGeom.grazingRad)));
    printf("Ray seperation              : %4.3f x %4.3f [ %4.3f x %4.3f ground plane] metres between adjacent rays\n",
           (dAz * centreRange), (dEl * centreRange),
           (dAz * centreRange), (dEl * centreRange / sin(cGeom.grazingRad)));
    
    double TxPowPerRay, gainRx ;
    TxPowPerRay = TxPowerPerRay(dAz, dEl, &gainRx);
    printf("EIRP                        : %e Watts (%f dBW)\n",TxPowPerRay,10*log(TxPowPerRay));
    double TB = TxPowPerRay * hdr.pulse_length * hdr.chirp_gamma * hdr.pulse_length ;
    TxPowPerRay = TxPowPerRay * TB ;
    
    // Initialise OpenCL and load the relevent information into the
    // platform structure. OCL tasks will use this later
    //
    cl_int      err;
    cl_uint     ndevs;
    cl_ulong    memSizeTmp,devMemSize;
#define GPUCAPABILITY_MAJOR 2
#define GPUCAPABILITY_MINOR 0
    OCLPlatform platform;

    platform.clSelectedPlatformID = NULL;
    err = oclGetPlatformID (&platform, &status);
    if(err != CL_SUCCESS){
        printf("Error: Failed to find a suitable OpenCL Launch platform\n");
        exit(-1);
    }
    oclPrintPlatform(platform);
    cl_context_properties props[] = {
        CL_CONTEXT_PLATFORM, (cl_context_properties)platform.clSelectedPlatformID,
        0
    };
    platform.props = props ;
    
    // Find number of devices on this platform
    //
    int useGPU = 0;
    if (useGPU){
        err = oclGetNamedGPUDevices(&platform, "NVIDIA", "",&platform.device_ids , &ndevs, &status);
        if(err == CL_SUCCESS){
            char cbuf[1024];
            cl_uint max_compute_units;
            printf("GPU DEVICES                 : %d\n",ndevs);
            for(int d=0; d<ndevs; d++){
                printf("DEVICE                      : %d\n",d);
                clGetDeviceInfo(platform.device_ids[d], CL_DEVICE_VENDOR, sizeof(cbuf), &cbuf, NULL);
                printf("  DEVICE VENDOR             : %s\n",cbuf);
                clGetDeviceInfo(platform.device_ids[d], CL_DEVICE_NAME, sizeof(cbuf), &cbuf, NULL);
                printf("  DEVICE NAME               : %s\n",cbuf);
                clGetDeviceInfo(platform.device_ids[d], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &max_compute_units, NULL);
                printf("  DEVICE MAX COMPUTE UNITS  : %d\n",max_compute_units);
            }
        }else{
            printf("No GPU devices available with compute capability: %d.%d\n", GPUCAPABILITY_MAJOR, GPUCAPABILITY_MINOR);
            exit(1);
        }
    }else{
        err = oclGetCPUDevices(&platform, &platform.device_ids, &ndevs, &status);
        if(err == CL_SUCCESS){
            printf("CPU DEVICES                 : %d\n",ndevs);
        }
    }
    
    // Find out the maximum amount of memory that all OpenCL devices have
    // For this task
    //
    devMemSize = 100e9;    // big number 100 GB
    for (int dev=0; dev<ndevs; dev++){
        err = clGetDeviceInfo(platform.device_ids[dev], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &memSizeTmp, NULL);
        if(err != CL_SUCCESS){
            printf("Error [%d] : couldn't obtain device info\n",err);
        }
        devMemSize = (memSizeTmp < devMemSize) ? memSizeTmp : devMemSize;
    }
    if (bounceToShow != 0) ndevs = 1;
    platform.nDevs = ndevs ;
    
    printf("DEVICES USED                : %d\n",ndevs);
    
    pthread_t *threads;
    threads = (pthread_t *)malloc(sizeof(pthread_t)*ndevs) ;
    if (threads == NULL) {
        printf("Error : Failed to malloc %d threads\n",ndevs);
        exit(-1);
    }
    
    threadData *threadDataArray ;
    threadDataArray = (threadData *)malloc(sizeof(threadData)*ndevs);
    if (threadDataArray == NULL) {
        printf("Error : Failed to malloc %d threadDataArrays\n",ndevs);
        exit(-1);
    }
    
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    SPImage cphd ;

    if(outCPHDFile != NULL){
        im_init(&cphd, &status) ;
        newhdr.data_type = ITYPE_CMPL_FLOAT ;
        im_create(&cphd, ITYPE_CMPL_FLOAT, newhdr.nsamp, newhdr.num_azi, 1.0, 1.0, &status) ;
    }
    // Calculate how many pulses can be processed on a device at a time
    //
    int pulsesPerDevice;
    
    while (nPulses % ndevs != 0) nPulses-- ;
    pulsesPerDevice = nPulses / ndevs ;
    
    for (int dev=0; dev<ndevs; dev++) {
        
        threadDataArray[dev].devIndex               = dev ;
        threadDataArray[dev].nThreads               = ndevs ;
        threadDataArray[dev].platform               = platform ;
        threadDataArray[dev].SceneBoundingBox       = SceneBoundingBox;
        threadDataArray[dev].startPulse             = startPulse + dev*pulsesPerDevice ;
        threadDataArray[dev].nPulses                = pulsesPerDevice ;
        threadDataArray[dev].nAzBeam                = nAzBeam ;
        threadDataArray[dev].nElBeam                = nElBeam ;
        threadDataArray[dev].cphdhdr                = &newhdr ;
//        threadDataArray[dev].TxPositions            = TxPos ;           // Pointer to beginning of TxPos data
//        threadDataArray[dev].RxPositions            = RxPos ;           // Pointer to beginning of RxPos data
//        threadDataArray[dev].Fx0s                   = Fx0s ;            // Pointer to beginning of Fx0s data
//        threadDataArray[dev].FxSteps                = FxSteps;          // Pointer to beginning of FxSteps data
//        threadDataArray[dev].amp_sf0                = amp_sf0 ;         // Pointer to beginning of amp_sf0 data
        threadDataArray[dev].gainRx                 = gainRx;
        threadDataArray[dev].PowPerRay              = TxPowPerRay ;
        if(outCPHDFile == NULL){
            threadDataArray[dev].phd                    = NULL ;
        }else{
            threadDataArray[dev].phd                    = &cphd;            // Pointer to beginning of cphd data
        }
//        threadDataArray[dev].chirpRate              = hdr.chirp_gamma;
//        threadDataArray[dev].ADRate                 = hdr.clock_speed ;
//        threadDataArray[dev].pulseDuration          = hdr.pulse_length ;
//        threadDataArray[dev].oneOverLambda          = hdr.freq_centre / SIPC_c ;
//        threadDataArray[dev].freq_centre            = hdr.freq_centre ;
//        threadDataArray[dev].StartFrequency         = hdr.freq_centre - (hdr.pulse_length * hdr.chirp_gamma / 2) ;
        threadDataArray[dev].bounceToShow           = bounceToShow-1;
        threadDataArray[dev].status                 = status ;
        threadDataArray[dev].interrogate            = interrogate ;
        threadDataArray[dev].interogPt              = interogPt ;
        threadDataArray[dev].interogRad             = interogRad ;
        threadDataArray[dev].interogFP              = &interrogateFP ;
        threadDataArray[dev].sceneMesh              = &baseMesh ;
        threadDataArray[dev].moverMesh              = &moverMesh ;
        
        if (bounceToShow)printf("\n+++++++++++++++++++++++++++++++++++++++\n");
        if (interrogate){
            time_t rawtime;
            struct tm * timeinfo;
            time ( &rawtime );
            timeinfo = localtime ( &rawtime );
            SPVector intOutRg, intRetRg ;
            double intRg, intMinR, intMaxR ;
            
            VECT_SUB(interogPt, TxPos[0], intOutRg);
            VECT_SUB(RxPos[0], interogPt, intRetRg);
            intRg = (VECT_MAG(intOutRg) + VECT_MAG(intRetRg))/2.0;
            intMinR = intRg - interogRad/2.0 ;
            intMaxR = intRg + interogRad/2.0 ;
            
            fprintf(interrogateFP, "\tInterrogate Output (Sarcastic %s)\n",FULL_VERSION);
            fprintf(interrogateFP, "Time of run             : %s\n",asctime (timeinfo));
            fprintf(interrogateFP, "Interrogation point     : %06.3f,%06.3f,%06.3f\n",interogPt.x,interogPt.y,interogPt.z);
            fprintf(interrogateFP, "Slant range to point (m): %06.3f --  %06.3f -- %06.3f\n", intMinR, intRg, intMaxR);
            fprintf(interrogateFP, "Interrogation Pt Radius : %06.3f\n",interogRad);
            fprintf(interrogateFP, "Interrogation Pulse(s)  : %d - %d\n",startPulse/pulseUndersampleFactor,(startPulse/pulseUndersampleFactor)+nPulses);
            fprintf(interrogateFP, "Range\t\tPower\t\tbounce\tTriangle\tHitPoint\n");
            fprintf(interrogateFP, "--------------------------------------------------------------------\n");
        }
        // Create thread data for each device
        //
        rc = pthread_create(&threads[dev], NULL, devPulseBlock, (void *) &threadDataArray[dev]) ;
        if (rc){
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }
    
    pthread_attr_destroy(&attr);
    void * threadStatus ;
    for (int dev=0; dev<ndevs; dev++) {
        rc = pthread_join(threads[dev], &threadStatus);
        if (rc) {
            printf("ERROR; return code from pthread_join() is %d\n", rc);
            exit(-1);
        }
    }
    if (bounceToShow)printf("\n+++++++++++++++++++++++++++++++++++++++\n");
    if (interrogate){
        fprintf(interrogateFP, "--------------------------------------------------------------------\n");
        fclose(interrogateFP);
    }
    
    endTimer(&runTimer, &status) ;
    printf("Done in " BOLD BLINK GREEN " %f " RESETCOLOR "seconds \n",timeElapsedInSeconds(&runTimer, &status)) ;
 
    if (outCPHDFile != NULL ) {
        FILE *fp;
        fp = fopen(outCPHDFile, "w") ;
        if (fp == NULL) {
            fprintf(stderr, "Failed to open file %s\n", outCPHDFile);
            perror("Error opening file");
            exit(888);
        }
        for (int p = 0; p < nPulses; p++){
            newhdr.pulses[p].sat_ps_rx = RxPos[p] ;
            newhdr.pulses[p].sat_ps_tx = TxPos[p] ;
        }
        printf("Writing CPHD File \"%s\"....",outCPHDFile);
        writeCPHD3Header( &newhdr, fp, &status ) ;
        write_cphd3_nb_vectors(&newhdr, 0, fp, &status) ;
        write_cphd3_wb_vectors(&newhdr, &cphd, 0, fp, &status) ;
        printf("...Done\n");
    }
    
    if(outCPHDFile != NULL){
        im_destroy(&cphd, &status);
    }
    im_close_lib(&status);
    free ( threadDataArray );
    free ( threads ) ;
    free(RxPos) ;
    free(TxPos) ;
    
    return 0;
    
}

    /*
    
    // Define these arrays here. The Raytrace routibe will expand then as we
    // proceed through pulses. They will then be handed to the PO Code
    // to perform the field calculations all at once.
    //
    HitPoint *hitPoints ;
    Ray *incidentRays, *observationRays ;
    int nHits ;
    
    // for each pulse in file
    //
    double t0 = newhdr.pulses[0].sat_tx_time ;
    for (int p=0; p<nPulses; ++p) {
        
        double t = newhdr.pulses[p].sat_tx_time  - t0;
        
        SPVector S0,S1,S2,S;
        VECT_CREATE(0.0, 0.0, 0.0, S0) ;  // Translation
        VECT_CREATE(0.0, 1.0, 0.0, S1) ;  // Velocity
        VECT_CREATE(0.0, 0.0, 0.0, S2) ;  // Acceleration
        
        
        // move the movers to the location for this pulse
        //
        S.x = S0.x + (S1.x * t) + (0.5 * S2.x * t * t) ;
        S.y = S0.y + (S1.y * t) + (0.5 * S2.y * t * t) ;
        S.z = S0.z + (S1.z * t) + (0.5 * S2.z * t * t) ;
        
        TriangleMesh mesh_t = moversMesh ;
        for(int i=0; i<mesh_t.vertices.size(); ++i){
            mesh_t.vertices[i].x += S.x ;
            mesh_t.vertices[i].y += S.y ;
            mesh_t.vertices[i].z += S.z ;
        }
    
        // Add movers to base scene
        //
        TriangleMesh newMesh = baseMesh.add(&mesh_t) ;
    
        // Build kdTree
        //
        kdTree::KdData * tree;
        int treeSize;

        kdTree::buildTree(&newMesh, &tree, &treeSize, (kdTree::TREEOUTPUT)(kdTree::OUTPUTDATA | kdTree::OUTPUTSUMM)) ;
    
        // Initialise the tree and build ropes and boxes to increase efficiency when traversing
        //
        kdTree::KdData *node;
        
        node = &(kdTree[0]) ;
        
        // build Ropes and Boxes
        //
        int Ropes[6] ;
        AABB sceneAABB = kdTree[0].brch.aabb ;
        for(int i=0; i<6; i++) Ropes[i] = NILROPE;
        BuildRopesAndBoxes(node, Ropes, sceneAABB, kdTree);
        
        ATS *accelTriangles;
        accelerateTriangles(mesh,&accelTriangles) ;
        
        
        
        threadCore 
        
        
        // Ray trace it
        //
       
        
          // end for
    //
    }
    
// Write pulse to file
    //

    

    return 0;
}
*/



