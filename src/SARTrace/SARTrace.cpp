/***************************************************************************
 *
 *       Module:    SARTrace.c
 *      Program:    SARTRACE
 *   Created by:    Darren on 27/07/2013.
 *                  Copyright (c) 2013 Dstl. All rights reserved.
 *
 *   Description:
 *      Uses the engine in Sarcastic just to ray trace intersection points 
 *      in a KdTree volume
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
#include "SARTrace.hpp"
extern "C" {
#include "ecef2SceneCoords.h"
}
#include "colourCodes.h"
#include "SartraceVersion.h"
#include "TriangleMesh.hpp"
#include "tryReadFile.hpp"

#define ROOTPATH "/tmp"

int main (int argc, char **argv){
    
    CPHDHeader  hdr ;
    OCLPlatform platform;
    
    double maxEl, maxAz, minEl, minAz ;
    int nPulses ;
    int useGPU ;
    
    cl_int      err;
    cl_uint     ndevs;
    cl_ulong    memSizeTmp,devMemSize;

    nPulses = 1;
    useGPU  = 0;
    
    SPStatus status ;
    im_init_status(status, 0) ;
    im_init_lib(&status, (char *)"sartrace", argc, (char **)argv);
	CHECK_STATUS_NON_PTR(status);
    
    SARTracebanner();
    
    char *outDir;
    char *meshFile, *inCPHDFile ;
    int pulseToTrace;
    int nRaysX, nRaysY;
    getSARTraceUserInput(&inCPHDFile, &meshFile, &outDir, &pulseToTrace,&nRaysX,&nRaysY, &status) ;
    
    // Start timing after user input
    //
    Timer runTimer ;
    startTimer(&runTimer, &status) ;
    
    // Read in Triangle Mesh
    //
    TriangleMesh    mesh ;
    AABB            SceneBoundingBox;
    printf("Reading in file %s...\n",meshFile);
    mesh.readPLYFile(meshFile);
    printf("Done. Checking file Integrity...\n");
    mesh.checkIntegrityAndRepair();
    mesh.buildTriangleAABBs();
    mesh.buildTrianglelCentres();
    printf("Done \n");

    // Find bounding box and Rx Gain
    //
    maxEl = maxAz = minEl = minAz = 0.0 ;
    SPVector min,max; VECT_CREATE(9e10, 9e10, 9e10, min); VECT_CREATE(-9e10, -9e10, -9e10, max);
    for(int i=0; i<mesh.triangles.size(); ++i){
        for(int j=0; j<3; ++j){
            min.cell[j] = (mesh.AABBs[i].AA.cell[j] < min.cell[j]) ? mesh.AABBs[i].AA.cell[j] : min.cell[j] ;
            max.cell[j] = (mesh.AABBs[i].BB.cell[j] > max.cell[j]) ? mesh.AABBs[i].BB.cell[j] : max.cell[j] ;
        }
    }
    SceneBoundingBox.AA = min ; SceneBoundingBox.BB = max ;
    
    // Convert CPHDFile ecef coords to scene centred coordinates
    //
    readCPHDHeader(inCPHDFile, &hdr, &status);
    SPVector TxPos[1] ;
    SPVector RxPos[1] ;
    SPVector SRP;
    TxPos[0] = hdr.pulses[pulseToTrace].sat_ps_tx ;
    RxPos[0] = hdr.pulses[pulseToTrace].sat_ps_rx ;
    SRP      = hdr.pulses[pulseToTrace].srp ;
    
    ecef2SceneCoords(1, TxPos, SRP) ;
    ecef2SceneCoords(1, RxPos, SRP) ;
    

    // Initialise OpenCL and load the relevent information into the
    // platform structure. OCL tasks will use this later
    //
    platform.clSelectedPlatformID = NULL;
    err = oclGetPlatformID (&platform, &status);
    if(err != CL_SUCCESS){
        printf("Error: Failed to find a suitable OpenCL Launch platform\n");
        exit(-1);
    }
    cl_context_properties props[] = {
        CL_CONTEXT_PLATFORM, (cl_context_properties)platform.clSelectedPlatformID,
        0
    };
    platform.props = props ;
        
    // Find number of devices on this platform
    //
    err = oclGetCPUDevices(&platform, &platform.device_ids, &ndevs, &status);

    
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
    ndevs = 1;
    platform.nDevs = ndevs ;
    
    traceThreadData threadDataArray ;

    // Calculate how many pulses can be processed on a device at a time
    //
    threadDataArray.devIndex               = 0 ;
    threadDataArray.nThreads               = 1 ;
    threadDataArray.platform               = platform ;
    threadDataArray.SceneBoundingBox       = SceneBoundingBox;
    threadDataArray.startPulse             = pulseToTrace ;
    threadDataArray.nPulses                = 1 ;
    threadDataArray.nAzBeam                = nRaysX ;
    threadDataArray.nElBeam                = nRaysY ;
    threadDataArray.TxPositions            = TxPos ;           // Pointer to beginning of TxPos data
    threadDataArray.RxPositions            = RxPos ;           // Pointer to beginning of RxPos data
    threadDataArray.gainRx                 = 0 ;
    threadDataArray.PowPerRay              = 0 ;
    threadDataArray.freq_centre            = 0 ;
    threadDataArray.bounceToShow           = 0 ;
    threadDataArray.status                 = status ;
    threadDataArray.interrogate            = false ;
    threadDataArray.sceneMesh              = &mesh ;
    
    printf("\n+++++++++++++++++++++++++++++++++++++++\n");
    
    sartraceCore(threadDataArray, outDir);
    
    printf("\n+++++++++++++++++++++++++++++++++++++++\n");

    endTimer(&runTimer, &status) ;
    printf("Raycasting Completed in %f seconds \n",timeElapsedInSeconds(&runTimer, &status)) ;
    
    im_close_lib(&status);

    return 0;

}

