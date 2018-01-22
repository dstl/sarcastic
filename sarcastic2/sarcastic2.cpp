/***************************************************************************
 *
 *       Module:    sarcastic2.cpp
 *      Program:    SARCASTIC
 *   Created by:    Darren on 15/03/2017.
 *                  Copyright (c) 2017 Dstl. All rights reserved.
 *
 *   Description:
 *      Version of sarcastic that uses threads rather than openCL calls
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
}
#include "getUserInput.hpp"
#include "threadCore.hpp"
#include <thread>
#include <fftw3.h>


int main(int argc, const char * argv[]) {
    
    int rc;
    
    SPStatus status ;
    im_init_status(status, 0) ;
    im_init_lib(&status, (char *)"sarcastic2", argc, (char **)argv);
    CHECK_STATUS_NON_PTR(status);
    
    banner() ;

    CPHDHeader hdr ;
    TriangleMesh baseMesh ;
    TriangleMesh moverMesh ;
    char *outCPHDFile ;
    int startPulse, nPulses, bounceToShow, nAzBeam, nElBeam, interrogate, pulseUndersampleFactor ;
    SPVector interogPt ;
    double interogRad ;
    FILE *interrogateFP ;
    kdTree::KdData * tree = NULL;
    ATS *accelTriangles = NULL;;
    int treeSize = 0;
    int polarisation ;
    int rayGenMethod ;
    int interogX,interogY ;
    
    getUserInput(&hdr, &baseMesh, &moverMesh, &outCPHDFile,
                     &startPulse, &nPulses, &bounceToShow, &nAzBeam, &nElBeam, &interrogate, &interogPt, &interogRad, &interogX, &interogY,
                 &interrogateFP, &pulseUndersampleFactor, &polarisation, &rayGenMethod, &status) ;
    
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
    
    // Set the polarisation in the output CPHDFile
    //
    if (polarisation == VV) {
        strcpy(*(hdr.polarisation), "VV");
    }else if (polarisation == VH){
        strcpy(*(hdr.polarisation), "VH");
    }else if (polarisation == HV){
        strcpy(*(hdr.polarisation), "HV");
    }else if (polarisation == HH){
        strcpy(*(hdr.polarisation), "HH");
    }else if (polarisation == V_){
        strcpy(*(hdr.polarisation), "V_");
    }else if (polarisation == H_){
        strcpy(*(hdr.polarisation), "H_");
    }else{
        printf("ERROR: Unknown polarisation type : %d \n", polarisation);
        exit(1);
    }

    // print out info about the collection
    //
    printCPHDCollectionInfo(&hdr, &status) ;
    printf("Wavelength                  : %f m\n",SIPC_c/newhdr.freq_centre);
    printf("Slant range resolution      : %f m\n",SIPC_c/(2*newhdr.pulse_length*newhdr.chirp_gamma));
    
    if(moverMesh.triangles.size() == 0){
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
    }
    
    // Rotate Rx and Tx Coords to be relative to scene centre
    //
    SPVector  *RxPos, *TxPos, SRP, zHat ;
    double    lambda ;
    
    RxPos = (SPVector *)sp_malloc(sizeof(SPVector)*nPulses) ;
    TxPos = (SPVector *)sp_malloc(sizeof(SPVector)*nPulses) ;

    for (int p = 0; p < nPulses; p++){
        RxPos[p] = newhdr.pulses[p].sat_ps_rx ;
        TxPos[p] = newhdr.pulses[p].sat_ps_tx ;
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
    
    SPImage cphd ;
    
    unsigned nThreads;
    if (bounceToShow != 0){
        nThreads = 1;
    }else{
        nThreads = std::thread::hardware_concurrency() ;
    }
    int pulsesPerThread;
    int nVec = newhdr.num_azi ;
    if (nVec < nThreads) {
        nThreads = nVec;
    }
        
    while (nVec % nThreads != 0) nVec-- ;
    newhdr.num_azi = nVec ;
    pulsesPerThread = nVec / nThreads ;
    
    if(outCPHDFile != NULL){
        im_init(&cphd, &status) ;
        newhdr.data_type = ITYPE_CMPL_FLOAT ;
        im_create(&cphd, ITYPE_CMPL_FLOAT, newhdr.nsamp, newhdr.num_azi, 1.0, 1.0, &status) ;
    }
    
    pthread_t *threads;
    threads = (pthread_t *)sp_malloc(sizeof(pthread_t)*nThreads) ;
    threadData *coreData;
    coreData = (threadData *)sp_malloc(sizeof(threadData) * nThreads) ;
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    fftwf_init_threads();

    
    for(int t=0; t<nThreads; ++t){
        
        coreData[t].tid                    = t ;
        coreData[t].nThreads               = nThreads ;
        coreData[t].SceneBoundingBox       = SceneBoundingBox;
        coreData[t].startPulse             = t*pulsesPerThread;
        coreData[t].nPulses                = pulsesPerThread ;
        coreData[t].nAzBeam                = nAzBeam ;
        coreData[t].nElBeam                = nElBeam ;
        coreData[t].cphdhdr                = &newhdr ;
        coreData[t].gainRx                 = gainRx;
        coreData[t].PowPerRay              = TxPowPerRay ;
        if(outCPHDFile == NULL){
            coreData[t].phd                    = NULL ;
        }else{
            coreData[t].phd                    = &cphd;
        }
        coreData[t].bounceToShow           = bounceToShow-1;
        coreData[t].status                 = status ;
        coreData[t].interrogate            = interrogate ;
        coreData[t].interogPt              = interogPt ;
        coreData[t].interogRad             = interogRad ;
        coreData[t].interogFP              = &interrogateFP ;
        coreData[t].sceneMesh              = &baseMesh ;
        coreData[t].moverMesh              = &moverMesh ;
        coreData[t].tree                   = &tree ;
        coreData[t].accelTriangles         = &accelTriangles ;
        coreData[t].treesize               = treeSize ;
        coreData[t].polarisation           = polarisation ;
        coreData[t].rayGenMethod           = rayGenMethod ;
        
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
            fprintf(interrogateFP, "Time of run                                 : %s",asctime (timeinfo));
            fprintf(interrogateFP, "Interrogation point                         : %06.3f,%06.3f,%06.3f\n",interogPt.x,interogPt.y,interogPt.z);
            fprintf(interrogateFP, "Interrogation point in image (x,y) pixels   : %03d,%03d\n",interogX,interogY);
            fprintf(interrogateFP, "Slant range to point (near - mid - far) (m) : %06.3f --  %06.3f -- %06.3f\n", intMinR, intRg, intMaxR);
            fprintf(interrogateFP, "Transmitter location                        : %06.3f,%06.3f,%06.3f\n",
                    newhdr.pulses[0].sat_ps_tx.x,newhdr.pulses[0].sat_ps_tx.y,newhdr.pulses[0].sat_ps_tx.z);
            fprintf(interrogateFP, "Interrogation Pt Radius                     : %06.3f\n",interogRad);
            fprintf(interrogateFP, "Interrogation Pulse                         : %d\n",startPulse/pulseUndersampleFactor);
            fprintf(interrogateFP, "--------------------------------------------------------------------\n");
        }
        
        rc = pthread_create(&threads[t], NULL, devPulseBlock, (void *) &coreData[t]) ;
        if (rc){
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
//        devPulseBlock(&coreData) ;
    }
    
    pthread_attr_destroy(&attr);
    void * threadStatus ;
    for (int t=0; t<nThreads; t++) {
        rc = pthread_join(threads[t], &threadStatus);
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
        
        if (newhdr.version == '3') {
            writeCPHD3Header( &newhdr, fp, &status ) ;
            write_cphd3_nb_vectors(&newhdr, 0, fp, &status) ;
            write_cphd3_wb_vectors(&newhdr, &cphd, 0, fp, &status) ;
        }else if (newhdr.version == 'x'){
            writeCPHDXHeader(&newhdr, fp, &status);
            writeCPHDXNarrowband(&newhdr, 0, fp, &status);
            writeCPHDXWideband(&newhdr, &cphd, fp, &status) ;
        }else{
            printf("Unknown cphd version requested : %c\n",newhdr.version) ;
            printf("Not writing output file...\n");
        }
        printf("...Done\n");
        
        fclose(fp);

    }
    
    if(outCPHDFile != NULL){
        im_destroy(&cphd, &status);
    }
    if(moverMesh.triangles.size() == 0){
        free(tree);
        delete accelTriangles ;
    }
    
    im_close_lib(&status);
    free(RxPos) ;
    free(TxPos) ;
    free(outCPHDFile);
    
    return 0;
    
}



