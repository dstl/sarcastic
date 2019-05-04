/***************************************************************************
 * 
 *           Module :  main.cpp
 *          Program :  bircs
 *       Created by :  Darren Muff on  30/04/2017
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *   Description:
 *      Programm to calculate the bistatic RCS of a target model
 * 
 *   (c) Crown Copyright 2018 Defence Science and Technology Laboratory
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software")
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
 ***************************************************************************/

#include <iostream>
#include "sarcastic2.hpp"
#include <sarclib/sarclib.h>
#include "tryReadFile.hpp"
#include "TriangleMesh.hpp"
#include "buildTree.hpp"
#include "buildRopesAndBoxes.hpp"
#include "accelerateTriangles.hpp"
#include "threadCore.hpp"
#include "readMaterialFile.hpp"
#include <thread>

extern "C" {
// #include "TxPowerPerRay.h"
#include "ecef2SceneCoords.h"
//#include "OpenCLUtils.h"
#include "bircsBanner.h"
}

double TxPowerPerRay(double rayWidthRadians, double rayHeightRadians, double *receiverGain);

int main (int argc, char **argv){
    
    SPStatus status ;
    im_init_status(status, 0) ;
    im_init_lib(&status, (char *)"bircs", argc, (char **)argv);
    CHECK_STATUS_NON_PTR(status);
  
    TriangleMesh baseMesh ;
    AABB        SceneBoundingBox;

    
    SPVector unitBeamAz, unitBeamEl, rVect, zHat, interogPt ;
    double centreRange, maxEl, maxAz, minEl, minAz, dAz, dEl, lambda, maxBeamUsedAz, maxBeamUsedEl, interogRad ;
    char *baseScene ;
    int startPulse, bounceToShow, interrogate ;
    int nAzBeam, nElBeam, nVec, rc ;
    FILE *interrogateFP = NULL ;
    kdTree::KdData * tree = NULL;
    ATS *accelTriangles = NULL;;
    int treeSize = 0;
    
    SPVector  * RxPos         = NULL ;
    SPVector  * TxPos         = NULL ;
    SPVector  * results       = NULL ;
    
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
    int rayGenMethod = 1 ;
    
    // Get user input
    //
    
    ILLAZ = input_int("IllAz", (char *)"ILLAZ", (char *)"Illumination azimuth", ILLAZ);
    ILLINC = input_int("ILLINC", (char *)"ILLINC", (char *)"Illumination inclination", ILLINC);
    ILLRANGE = input_dbl("Illumination Range", (char *)"ILLRANGE", (char *)"Slant Range to scene centre for illumination", ILLRANGE);
    freq_centre = input_dbl("Centre Frequency", (char *)"freq", (char *)"Transmit centre frequency", freq_centre);
    obsdist = input_dbl("Slant range to scene centre", (char *)"SRdist", (char *)"Slant Range to scene centre", obsdist);
    phistart = input_dbl("phistart", (char *)"PHISTART", (char *)"Starting azimuth angle", phistart);
    phiend = input_dbl("phiend", (char *)"PHIEND", (char *)"end azimuth angle", phiend);
    nphis = input_int("nphis", (char *)"NPHIS", (char *)"Number of azimuth samples", nphis);
    thetastart = input_dbl("thetastart", (char *)"THETASTART", (char *)"starting angle of incidence", thetastart);
    thetaend = input_dbl("thetaend", (char *)"THETAEND", (char *)"Ending angle of incidence", thetaend);
    nthetas = input_int("nthetas", (char *)"NTHETAS", (char *)"Number of elevation samples", nthetas);
    if(nthetas * nphis == 1){
    bounceToShow = input_int("bounceToShow", (char *)"BOUNCETOSHOW", (char *)"Output just the ray intersections corresponding to this bounce. (0=no bounce output)", bounceToShow);
    }
    do{
        printf("There are four ways to cast rays in SARCASTIC :\n");
        printf("\t1 : TRIANGLECENTRE\n");
        printf("\t2 : RANDOMRAYS\n");
        printf("\t3 : FIRSTTIMERANDOM\n");
        printf("\t4 : PARALLELRANDOM\n");
        rayGenMethod = input_int("Enter number for ray generation method or '?' for help", "RayGenMethod",
          "SARCASTIC can generate rays in one of 4 ways. Here is an explanation of each method and why you should use it: \n\t1 : TRIANGLECENTRE \n\t\tThis method uses the centre of each triangle in the input mesh to determine the direction that\n\t\teach ray should be cast in. The origin is obviously the transmitter location for a given pulse.\n\t\tIt doesnt not try to perform any Z-buffering of triangles and so triangles that are behind another \n\t\tone will still generate a ray and will result in the nearer triangle being hit many times. (this gets \n\t\taccounted for in the RCS calculation and so doesnt affect the output.) Use this method to guarantee that\n\t\tevery triangle in the mesh (that can be illuminated by transmitted ray) is illuminated by the transmitted\n\t\tray. Use this method if: \n\t\t    You have a large scene with a smaller number of large triangles. \n\t\t    You want to make sure that every part of your model is illuminated \n\t\t    You have small triangles and so secondary and higher bounces will be incident on all faces of the mesh. \n\t2 : RANDOMRAYS \n\t\tThis method generates a Gaussian distribution of random rays. The size of the scene is measured first\n\t\tso that the entire scene is illuminated. If this method is selected then the number of rays in \n\t\tazimuth/cross-range and elevation will be asked for. Use this method if: \n\t\t    You want to perform a quick run and are not that bothered about illuminating every part of teh input mesh. \n\t\t    You have a small scene with large triangles. \n\t\t    You have large triangles and want to make sure there are many reflections from the surface of each triangle. \n\t3 : FIRSTIMERANDOM \n\t\tThis method is similar to method 2 in that it generates a Gaussian distribution of random rays. The difference\n\t\thowever is that after generating the aim point for the initial rays it then remembers the aim point for future\n\t\tpulses. This is useful if you want to make sure that pulse scattering centres are coherent from pulse to pulse.\n\t\tUse this method if: \n\t\t    You want to guarantee that a triangle correlates pulse to pulse over the entire SAR aperture \n\t4 : PARALLELRANDOM \n\t\tThis method is the same as method 2 in that it generates a Gaussian distribution of random rays. The difference\n\t\there is that the origin of each pulse is adjusted so that all the rays are parallel. Use this method if: \n\t\t    You are simulating a scene in the near field but would like it to be imaged in the far field. \n\n ",
            rayGenMethod);
    }while(rayGenMethod < 1 || rayGenMethod > 4);
    
    if (rayGenMethod == RANDOMRAYS || rayGenMethod == FIRSTTIMERANDOM || rayGenMethod == PARALLELRANDOM) {
        nAzBeam = input_int((char *)"Azimuth rays in radar beam?", (char *)"nAzBeam",
                             (char *)"Number of azimuth rays to use to construct radar beam. More is better but slower",nAzBeam);
        nElBeam = input_int((char *)"Elevation rays in radar beam?", (char *)"nElBeam",
                             (char *)"Number of elevation rays to use to construct radar beam. More is better but slower",nElBeam);
    }

    char *polstr ;
    bool validpol ;
    int pol = 0 ;
    do{
        validpol = false ;
        polstr = input_string("Enter polarisation to simulate", "Polarisation",
                              "Options are \'VV\',\'VH\',\'HV\',\'HH\',\'V_\', and \'H_\'. If one of the last two are used then the received H and V fields will be combined","VV" ) ;
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
                            ,"delaunay.ply") ;
    
    // Read in the material properties file if required
    //
    char *matfile = input_string((char *)"Input materialfile filename", (char *)"materialfilename",
                                 (char *)"The name of a 'materialfile' or 'none' (defaults used)",
                                 (char *) MATERIALPROPS);
    initialiseMaterials(matfile, true);

    // Get the output name from the user, open the file for write and check that it all worked.
    //
    char * outfname = input_string("Output filename", "outputFileName",
				   "The name of the file to dump the results into - will be ASCII in form of RCS, Phi, Theta",
				   "bircs_output.txt");
    FILE * fp = fopen(outfname, "w");
    if (fp == NULL) {
      fprintf(stderr, "Failed to open file %s for write.\n", outfname);
      perror("Open failed: ");
      exit(1);
    }
    
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
    

    // Start timing after user input
    //
    Timer runTimer ;
    startTimer(&runTimer, &status) ;
    
    nVec = nphis*nthetas ;
    
    // Rotate Rx and Tx Coords to be relative to scene centre
    //
    RxPos   = (SPVector *)sp_malloc(sizeof(SPVector)*nVec);
    TxPos   = (SPVector *)sp_malloc(sizeof(SPVector)*nVec);
    results = (SPVector *)sp_malloc(sizeof(SPVector)*nVec);
    
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
    
    centreRange = VECT_MAG(TxPos[nVec/2]);
    VECT_MINUS( TxPos[nVec/2], rVect ) ;
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
    
    // double TxPowPerRay, gainRx ;
    // TxPowPerRay = TxPowerPerRay(dAz, dEl, &gainRx);
    double TxPowPerRay = 1.0, gainRx = 1.0 ;
    
    unsigned nThreads;
    if (bounceToShow != 0){
        nThreads = 1;
    }else{
        nThreads = std::thread::hardware_concurrency() ;
    }
    int pulsesPerThread;
    if (nVec < nThreads) {
        nThreads = nVec;
    }
    
    while (nVec % nThreads != 0) nVec-- ;
    pulsesPerThread = nVec / nThreads ;
    
    pthread_t *threads;
    threads = (pthread_t *)sp_malloc(sizeof(pthread_t)*nThreads) ;
    
    bircsThreadData *threadDataArray ;
    threadDataArray = (bircsThreadData *)sp_malloc(sizeof(bircsThreadData)*nThreads);
    
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    
    for (int t=0; t<nThreads; ++t) {
        
        threadDataArray[t].tid                    = t ;
        threadDataArray[t].nThreads               = nThreads ;
        threadDataArray[t].SceneBoundingBox       = SceneBoundingBox;
        threadDataArray[t].startPulse             = t*pulsesPerThread ;
        threadDataArray[t].nPulses                = pulsesPerThread ;
        threadDataArray[t].nAzBeam                = nAzBeam ;
        threadDataArray[t].nElBeam                = nElBeam ;
        threadDataArray[t].TxPositions            = TxPos ;           // Pointer to beginning of TxPos data
        threadDataArray[t].RxPositions            = RxPos ;           // Pointer to beginning of RxPos data
        threadDataArray[t].gainRx                 = gainRx;
        threadDataArray[t].PowPerRay              = TxPowPerRay ;
        threadDataArray[t].freq_centre            = freq_centre ;
        threadDataArray[t].bounceToShow           = bounceToShow-1;
        threadDataArray[t].status                 = status ;
        threadDataArray[t].interrogate            = interrogate ;
        threadDataArray[t].interogPt              = interogPt ;
        threadDataArray[t].interogRad             = interogRad ;
        threadDataArray[t].interogFP              = &interrogateFP ;
        threadDataArray[t].sceneMesh              = &baseMesh ;
        threadDataArray[t].tree                   = &tree ;
        threadDataArray[t].accelTriangles         = &accelTriangles ;
        threadDataArray[t].treesize               = treeSize ;
        threadDataArray[t].polarisation           = pol ;
        threadDataArray[t].results                = results ;
        threadDataArray[t].rayGenMethod           = rayGenMethod ;
        
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
            fprintf(interrogateFP, "Interrogation Pulse(s)  : %d - %d\n",startPulse,startPulse+nVec);
            fprintf(interrogateFP, "Range\t\tPower\t\tbounce\tTriangle\tHitPoint\n");
            fprintf(interrogateFP, "--------------------------------------------------------------------\n");
        }
        // Create thread data for each device
        //
        rc = pthread_create(&threads[t], NULL, devPulseBlock, (void *) &threadDataArray[t]) ;
        if (rc){
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
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
    
    float maxR = 0.0 ;
    float maxphi = 0.0 ;
    float maxtheta = 0.0;
    printf("\n--------------------------------------------------------------------\n");
    for (int i=0; i<nVec; ++i){
        if (results[i].r > maxR ) {
            maxR     = results[i].r ;
            maxphi   = results[i].phi ;
            maxtheta = results[i].theta ;
        }
        fprintf(fp, "%f, %f, %f\n",results[i].r,results[i].phi,results[i].theta);
    }

    fclose(fp);
    
    double lambda_squared = lambda * lambda ;
    printf("Maximum RCS was %f dB m^2 (%f m^2) at az:%fdeg, incidence: %fdeg\n",maxR,pow(10.0,maxR/10.0),RAD2DEG(maxphi), RAD2DEG(maxtheta) );
    printf("    1m^2 flat plate               : %f (%f dB m^2)\n",4*SIPC_pi / lambda_squared,10*log10(4*SIPC_pi / lambda_squared)); // 8620.677
    printf("    1m dihedral                   : %f (%f dB m^2)\n",8*SIPC_pi / lambda_squared,10*log10(8*SIPC_pi / lambda_squared)); // 17241.354
    printf("    1m trihedral                  : %f (%f dB m^2)\n",12*SIPC_pi/ lambda_squared,10*log10(12*SIPC_pi/ lambda_squared)); // 25862.031
    printf("    1m dia Sphere                 : \n");
    // printf("        Rayleigh Region r<<lambda : %f ( %f dB m^2)\n", 9*SIPC_pi*0.5*0.5*pow(k*0.5,4), 10*log10(9*SIPC_pi*0.5*0.5*pow(k*0.5,4)));
    // printf("        Optical Region  r>>lambda : %f ( %f dB m^2)\n", SIPC_pi*0.5*0.5, 10*log10(SIPC_pi*0.5*0.5));
    // printf("         ( r = 0.5 m, lambda = %5.3f m )\n",lambda);
    printf("\n--------------------------------------------------------------------\n");


    double runtime = timeElapsedInSeconds(&runTimer, &status);
    printf("BIRCS completed in      : %4.4f secs \n",runtime);
    printf("Measurements per second : %f \n", runtime / (nphis * nthetas)) ;
    if(rayGenMethod == 1){
        printf("Total Rays cast         : %3.2f million \n", (double)((int)baseMesh.triangles.size() * nphis * nthetas) / 1000000.0);
        printf("Rays per second         : %4.2e \n", ((int)baseMesh.triangles.size() * nphis * nthetas) / runtime );
    }else{
        printf("Total Rays cast         : %3.2f million \n", (double)(nAzBeam * nElBeam) * nphis * nthetas / 1000000.0);
        printf("Rays per second         : %4.2e \n", ((nAzBeam * nElBeam) * nphis * nthetas) / runtime );
    }
    printf("\n--------------------------------------------------------------------\n");

    im_close_lib(&status);
    free ( threadDataArray );
    free ( threads ) ;
    free(RxPos) ;
    free(TxPos) ;
    free(results) ;
    free(tree);
    delete accelTriangles ;
    
    return 0;
    
}

