/***************************************************************************
 *
 *       Module:    main.cpp
 *      Program:    BIRCS
 *   Created by:    Darren on 30/04/2017.
 *                  Copyright (c) 2013 Dstl. All rights reserved.
 *
 *   Description:
 *      Programm to calculate the bistatic RCS of a target model
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

#include <iostream>
#include "sarcastic2.hpp"
#include <SIlib2/SIlib2.h>
#include "tryReadFile.hpp"
#include "TriangleMesh.hpp"
#include "buildTree.hpp"
#include "buildRopesAndBoxes.hpp"
#include "accelerateTriangles.hpp"
#include "threadCore.hpp"

extern "C" {
#include "TxPowerPerRay.h"
#include "ecef2SceneCoords.h"
#include "OpenCLUtils.h"
#include "bircsBanner.h"
}

#define ROOTPATH "/tmp"

double TxPowerPerRay(double rayWidthRadians, double rayHeightRadians, double *receiverGain);

int main (int argc, char **argv){
    
    SPStatus status ;
    im_init_status(status, 0) ;
    im_init_lib(&status, (char *)"bircs", argc, (char **)argv);
    CHECK_STATUS_NON_PTR(status);
  
    TriangleMesh baseMesh ;
    AABB        SceneBoundingBox;
    OCLPlatform platform;
#define GPUCAPABILITY_MAJOR 2
#define GPUCAPABILITY_MINOR 0
    
    SPVector unitBeamAz, unitBeamEl, rVect, zHat, interogPt ;
    double centreRange, maxEl, maxAz, minEl, minAz, dAz, dEl, lambda, maxBeamUsedAz, maxBeamUsedEl, interogRad ;
    char *baseScene ;
    int startPulse, nPulses, bounceToShow, interrogate ;
    int nAzBeam, nElBeam, nVec, dev, rc ;
    FILE *interrogateFP = NULL ;
    
    cl_int      err;
    cl_uint     ndevs;
    cl_ulong    memSizeTmp,devMemSize;
    
    SPVector  * RxPos         = NULL ;
    SPVector  * TxPos         = NULL ;
    
    bircsBanner() ;
    
    bounceToShow = 0;
    nAzBeam = 100;
    nElBeam = 100;
    interrogate = 0;
    
    double phistart =  0.0;
    double phiend = 360.0;
    double thetastart = 0.0;
    double thetaend = 180.0;
    double obsdist = 20000.0;
    double ILLAZ = 180.0;
    double ILLINC = 45.0;
    double ILLRANGE = 20000.0;
    int nphis = 180;
    int nthetas = 90;
    double freq_centre = 1e9 ;
    
    
    ILLAZ = input_int("IllAz", (char *)"ILLAZ", (char *)"Illumination azimuth", ILLAZ);
    ILLINC = input_int("ILLINC", (char *)"ILLINC", (char *)"Illumination inclination", ILLINC);
    freq_centre = input_dbl("Centre Frequency", (char *)"freq", (char *)"Transmit centre frequency", freq_centre);
    obsdist = input_dbl("Slant range to scene centre", (char *)"SRdist", (char *)"Slant Range to scene centre", obsdist);
    phistart = input_dbl("phistart", (char *)"PHISTART", (char *)"Starting azimuth angle", phistart);
    phiend = input_dbl("phiend", (char *)"PHIEND", (char *)"end azimuth angle", phiend);
    nphis = input_int("nphis", (char *)"NPHIS", (char *)"Number of azimuth samples", nphis);
    thetastart = input_dbl("thetastart", (char *)"THETASTART", (char *)"starting angle of incidence", thetastart);
    thetaend = input_dbl("thetaend", (char *)"THETAEND", (char *)"Ending angle of incidence", thetaend);
    nthetas = input_int("nthetas", (char *)"NTHETAS", (char *)"Number of elevation samples", nthetas);
    nAzBeam = input_int("nAzBeam", (char *)"NAZBEAM", (char *)"Number of azimuth rays to cast for each observation position", nAzBeam);
    nElBeam = input_int("nElBeam", (char *)"NELBEAM", (char *)"Number of elevation rays to cast for each observation position", nElBeam);
    bounceToShow = input_int("bounceToShow", (char *)"BOUNCETOSHOW", (char *)"Output just the ray intersections corresponding to this bounce. (0=no bounce output)", bounceToShow);
    char *polstr ;
    bool validpol ;
    int pol = 0 ;
    do{
        validpol = false ;
        polstr = input_string("Enter polarisation to simulate", "Polarisation",
                              "Options are \'VV\',\'VH\',\'HV\',\'HH\',\'V_\', and \'H_\'. If one of the last two are used then the received H and V fields will be combined",polstr ) ;
        if (!strcasecmp(polstr, "vv")) {
            pol = VV ;
            validpol = true ;
        }else if (!strcasecmp(polstr, "vh")){
            pol = VH ;
            validpol = true ;
        }else if (!strcasecmp(polstr, "hv")){
            pol = HV ;
            validpol = true ;
        }else if (!strcasecmp(polstr, "hh")){
            pol = HH ;
            validpol = true ;
        }else if (!strcasecmp(polstr, "v_")){
            pol = V_ ;
            validpol = true ;
        }else if (!strcasecmp(polstr, "h_")){
            pol = H_ ;
            validpol = true ;
        }else{
            printf("Invalid polarisation. Options are \'VV\',\'VH\',\'HV\',\'HH\',\'V_\' or \'H_\'\n");
        }
    }while(!validpol);
    baseScene = tryReadFile("Name of Base scene", "baseScene",
                            "Enter the name of a file that will be the base scene to be raytraced. The file must be in a .PLY file format"
                            ,ROOTPATH"/delaunay.ply") ;
    
    // Read in the triangle mesh from the input plyfile and check it's
    // integrity
    //
    printf("Reading in file %s...\n",baseScene);
    baseMesh.readPLYFile(baseScene);
    printf("Done. Checking file Integrity...\n");
    baseMesh.checkIntegrityAndRepair();
    baseMesh.buildTriangleAABBs();
    baseMesh.buildTrianglelCentres();
    printf("Done \n");
    
    ILLINC = DEG2RAD(ILLINC);
    ILLAZ = DEG2RAD(ILLAZ);
    phistart = DEG2RAD(phistart);
    phiend = DEG2RAD(phiend);
    thetastart = DEG2RAD(thetastart);
    thetaend = DEG2RAD(thetaend) ;
    
    nPulses = nphis * nthetas ;
    double deltaiphi, deltaitheta;
    deltaiphi = (phiend-phistart) / nphis ;
    deltaitheta = (thetaend-thetastart) / (nthetas) ;
    SPVector illOrigin, illDir;
    double illRange, illAz, illInc ;
    illRange = ILLRANGE ;
    illAz    = ILLAZ ;
    illInc   = ILLINC ;
    VECT_CREATE(illRange*sin(illInc)*cos(illAz), illRange*sin(illInc)*sin(illAz), illRange*cos(illInc), illOrigin) ;
    VECT_NORM(illOrigin, illDir) ;
    
    /*
    kdTree::buildTree(baseMesh, &tree, &treeSize, (kdTree::TREEOUTPUT)(kdTree::OUTPUTSUMM)) ;
    accelerateTriangles(&baseMesh,&accelTriangles) ;
    // Initialise the tree and build ropes and boxes to increase efficiency when traversing
    //
    kdTree::KdData *node;
    node = &(tree[0]) ;
    // build Ropes and Boxes
    //
    AABB sceneAABB ;
    sceneAABB = tree[0].brch.aabb ;
    int Ropes[6] ;
    for(int i=0; i<6; i++) Ropes[i] = NILROPE;
    BuildRopesAndBoxes(node, Ropes, sceneAABB, tree);
    */
    
    // Start timing after user input
    //
    Timer runTimer ;
    startTimer(&runTimer, &status) ;
    
    nVec = nphis*nthetas ;
    
    // Rotate Rx and Tx Coords to be relative to scene centre
    //
    RxPos = (SPVector *)sp_malloc(sizeof(SPVector)*nPulses);
    TxPos = (SPVector *)sp_malloc(sizeof(SPVector)*nPulses);
    
    double theta_s, phi_s;
    SPVector obsDir;
    for (int t=0; t<nthetas; t++){
        theta_s = thetastart + t * deltaitheta ;
        for(int p=0; p<nphis; p++){
            phi_s = phistart + p * deltaiphi ;
            obsDir.x = sin(theta_s) * cos(phi_s);
            obsDir.y = sin(theta_s) * sin(phi_s);
            obsDir.z = cos(theta_s);
            VECT_SCMULT(obsDir, obsdist, RxPos[t*nphis+p]) ;
            TxPos[t*nphis+p] = illOrigin ;
        }
    }
    
    // Calculate azbeamwidth from scene
    // Calculate elbeamwidth from scene
    // use largest so that beam sampling is uniform
    // multiply by sceneBeamMargin
    // set to be beamwidth
    
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
    
    dAz = maxBeamUsedAz / nAzBeam;
    dEl = maxBeamUsedEl / nElBeam;
    lambda = SIPC_c / freq_centre ;
    
    double TxPowPerRay, gainRx ;
    TxPowPerRay = TxPowerPerRay(dAz, dEl, &gainRx);
    //    printf("EIRP                        : %e Watts (%f dBW)\n",TxPowPerRay,10*log(TxPowPerRay));
    
    // Initialise OpenCL and load the relevent information into the
    // platform structure. OCL tasks will use this later
    //
    
    platform.clSelectedPlatformID = NULL;
    err = oclGetPlatformID (&platform, &status);
    if(err != CL_SUCCESS){
        printf("Error: Failed to find a suitable OpenCL Launch platform\n");
        exit(-1);
    }
    //    oclPrintPlatform(platform);
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
            //            printf("GPU DEVICES                 : %d\n",ndevs);
            for(int d=0; d<ndevs; d++){
                //                printf("DEVICE                      : %d\n",d);
                clGetDeviceInfo(platform.device_ids[d], CL_DEVICE_VENDOR, sizeof(cbuf), &cbuf, NULL);
                //                printf("  DEVICE VENDOR             : %s\n",cbuf);
                clGetDeviceInfo(platform.device_ids[d], CL_DEVICE_NAME, sizeof(cbuf), &cbuf, NULL);
                //                printf("  DEVICE NAME               : %s\n",cbuf);
                clGetDeviceInfo(platform.device_ids[d], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &max_compute_units, NULL);
                //                printf("  DEVICE MAX COMPUTE UNITS  : %d\n",max_compute_units);
            }
        }else{
            printf("No GPU devices available with compute capability: %d.%d\n", GPUCAPABILITY_MAJOR, GPUCAPABILITY_MINOR);
            exit(1);
        }
    }else{
        err = oclGetCPUDevices(&platform, &platform.device_ids, &ndevs, &status);
        if(err == CL_SUCCESS){
            //            printf("CPU DEVICES                 : %d\n",ndevs);
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
    
    //    printf("DEVICES USED                : %d\n",ndevs);
    
    pthread_t *threads;
    threads = (pthread_t *)malloc(sizeof(pthread_t)*ndevs) ;
    if (threads == NULL) {
        printf("Error : Failed to malloc %d threads\n",ndevs);
        exit(-1);
    }
    
    bircsThreadData *threadDataArray ;
    threadDataArray = (bircsThreadData *)sp_malloc(sizeof(bircsThreadData)*ndevs);
    if (threadDataArray == NULL) {
        printf("Error : Failed to malloc %d threadDataArrays\n",ndevs);
        exit(-1);
    }
    
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    
    // Calculate how many pulses can be processed on a device at a time
    //
    int pulsesPerDevice;
    
    while (nPulses % ndevs != 0) nPulses-- ;
    pulsesPerDevice = nPulses / ndevs ;
    
    for (dev=0; dev<ndevs; dev++) {
        
        
        threadDataArray[dev].devIndex               = dev ;
        threadDataArray[dev].nThreads               = ndevs ;
        threadDataArray[dev].platform               = platform ;
        threadDataArray[dev].SceneBoundingBox       = SceneBoundingBox;
        threadDataArray[dev].startPulse             = startPulse + dev*pulsesPerDevice ;
        threadDataArray[dev].nPulses                = pulsesPerDevice ;
        threadDataArray[dev].nAzBeam                = nAzBeam ;
        threadDataArray[dev].nElBeam                = nElBeam ;
        threadDataArray[dev].TxPositions            = TxPos ;           // Pointer to beginning of TxPos data
        threadDataArray[dev].RxPositions            = RxPos ;           // Pointer to beginning of RxPos data
        threadDataArray[dev].gainRx                 = gainRx;
        threadDataArray[dev].PowPerRay              = TxPowPerRay ;
        threadDataArray[dev].freq_centre            = freq_centre ;
        threadDataArray[dev].bounceToShow           = bounceToShow-1;
        threadDataArray[dev].status                 = status ;
        threadDataArray[dev].interrogate            = interrogate ;
        threadDataArray[dev].interogPt              = interogPt ;
        threadDataArray[dev].interogRad             = interogRad ;
        threadDataArray[dev].interogFP              = &interrogateFP ;
        threadDataArray[dev].sceneMesh              = &baseMesh ;
        threadDataArray[dev].polarisation           = pol ;
        
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
            fprintf(interrogateFP, "Interrogation Pulse(s)  : %d - %d\n",startPulse,startPulse+nPulses);
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
    for (dev=0; dev<ndevs; dev++) {
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
    //    printf("Done in " BOLD BLINK GREEN " %f " RESETCOLOR "seconds \n",timeElapsedInSeconds(&runTimer, &status)) ;
    
    
    im_close_lib(&status);
    free ( threadDataArray );
    free ( threads ) ;
    free(RxPos) ;
    free(TxPos) ;
    
    return 0;
    
}

