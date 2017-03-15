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

int main(int argc, const char * argv[]) {
    
    char *prompt  = (char *)malloc(sizeof(char)*256);
    char * file ;

    
    SPStatus status ;
    im_init_status(status, 0) ;
    im_init_lib(&status, (char *)"sarcastic", argc, (char **)argv);
    CHECK_STATUS_NON_PTR(status);
    
    // Read in base scene
    //
    char *baseScene ;
    baseScene = tryReadFile("Name of Base scene", "baseScene",
                            "Enter the name of a file that will be the base scene to be raytraced. The file must be in a .PLY file format"
                            , baseScene) ;
    
    // Read in the triangle mesh from the input plyfile and check it's
    // integrity
    //
    TriangleMesh baseMesh;
    baseMesh.readPLYFile(baseScene);
    baseMesh.checkIntegrityAndRepair();
    baseMesh.buildTriangleAABBs();
    
    // Read in the movers scene
    //
    char *moversScene ;
    moversScene = tryReadFile("Name of movers scene", "moversScene",
                            "Enter the name of a file that will contain the things that move. The file must be in a .PLY file format"
                            , moversScene) ;
    
    TriangleMesh moversMesh;
    moversMesh.readPLYFile(moversScene);
    moversMesh.checkIntegrityAndRepair();
    moversMesh.buildTriangleAABBs();
    
    // read in the base cphdFile
    //
    char *inCPHDFile ;
    inCPHDFile = tryReadFile("Name of CPHD File", "CPHDFile", "Enter the name of teh CPHD file to use as input", inCPHDFile) ;
    
    CPHDHeader hdr;
    int startPulse, nPulses, pulseUnderSampleFactor, nAzBeam, nElBeam ;

    readCPHDHeader(inCPHDFile, &hdr, &status) ;
    
    // Get the start Pulse
    //
    do{
        startPulse = 0;
        sprintf(prompt, "Start Pulse (0-%d)",hdr.num_azi-2);
        startPulse = input_int(prompt, (char *)"startPulse", (char *)"Start pulse in CPHD file to reconstruct", startPulse);
    } while (!(startPulse >=0 && startPulse <= hdr.num_azi-2)) ;

    // Get the number of pulses
    //
    sprintf(prompt, "Number of pulses (1-%d)",hdr.num_azi - startPulse);
    nPulses = hdr.num_azi ;
    nPulses = input_int(prompt, (char *)"nPulses", (char *)"Number of pulses to reconstruct in cphdFile", nPulses);
    
    // Get the pulse undersample factor
    //
    do {
        pulseUnderSampleFactor = input_int("Pulse undersampling factor", (char *)"pulseUndersampFact",
                                           (char *)"Reduces the number of azimuth pulses that are processed. This effectively reduces the collection PRF. If the simulated scene size is small then the azimuth ambiguities will not fold in far enough to affect the simulated scene.",
                                           pulseUnderSampleFactor);
    } while (pulseUnderSampleFactor <= 0) ;
    
    // Reduce the Cphd data we have to deal with to make things quicker.
    //
    CPHDHeader newhdr = hdr ;
    newhdr.num_azi = nPulses / pulseUnderSampleFactor ;
    newhdr.pulses  = (CPHDPulse *)calloc(newhdr.num_azi, sizeof(CPHDPulse));
    for(int p=startPulse; p<startPulse+nPulses; p+=pulseUnderSampleFactor){
        newhdr.pulses[(p-startPulse)/pulseUnderSampleFactor] = hdr.pulses[p] ;
    }
    nPulses = newhdr.num_azi ;
    
    // Get the number of rays to trace
    //
    nAzBeam = nElBeam = 100 ;
    nAzBeam = input_int((char *)"Azimuth rays in radar beam?", (char *)"nAzBeam",
                         (char *)"Number of azimuth rays to use to construct radar beam. More is better but slower",nAzBeam);
    nElBeam = input_int((char *)"Elevation rays in radar beam?", (char *)"nElBeam",
                         (char *)"Number of elevation rays to use to construct radar beam. More is better but slower",nElBeam);

    // change extension
    //
    sprintf(prompt, "%s",inCPHDFile);
    char *pExt = strrchr(file, '.');
    if (pExt != NULL)
        strcpy(pExt, ".sarc.cph");
    else
        strcat(prompt, ".sarc.cph");
    
    SPStatus fileStat ;
    char *outCPHDFile;
    FILE *fp;
    do {
        im_init_status(fileStat, 0) ;
        outCPHDFile = input_string((char *)"Output CPHD Filename", (char *)"CPHDOut",
                                (char *)"The name of a CPHD file to create.",
                                    prompt);
        if ( (fp = fopen(outCPHDFile, "w")) == NULL){
            printf(RED "Cannot access file %s\n" RESETCOLOR,outCPHDFile);
            fileStat.status = BAD_FILE ;
        }else fclose(fp) ;
        
    } while(fileStat.status != NO_ERROR);
    
    
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
    
        // Ray trace it
        //
        HitPoint *hitPoints ;
        Ray *incidentRays, *observationRays ;
        int nHits ;
        rangeAndPower *rnp;
        
        rayTrace(tree, treeSize, &newMesh, &hitPoints, &incidentRays, &observationRays, &nHits) ;

        // calculate Fields
        //
        rnp = new rangeAndPower [nHits] ;
        
        POFields(&newMesh, hitPoints, incidentRays, observationRays, nHits, &rnp) ;
        
        // Write pulse to file
        //
    
    // end for
    //
    }
    

    return 0;
}

void rayTrace(kdTree::KdData *tree, int treeSize, TriangleMesh *newMesh, SPVector **hitPoints, Ray **incidentRays, Ray **observationRays, int *nHits) {
    
}

void POFields(const TriangleMesh *mesh, const HitPoint *hitPoints, const Ray *incidentRays, const Ray *observationRays, const int nHits, rangeAndPower **rnp) {
    
    
    
}



