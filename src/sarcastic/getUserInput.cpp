/***************************************************************************
 * 
 *           Module :  getUserInput.cpp
 *          Program :  sarcastic
 *       Created by :  Darren Muff on 15/03/2017
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      Simple routone to interrogate the user and see if teh programme can get
 *      some sense out of them regarding what they would like to achieve
 *
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

#include "getUserInput.hpp"
#include "colourCodes.h"
#include "tryReadFile.hpp"
#include "readMaterialFile.hpp"
#include "sceneExtent.hpp"

int getUserInput(CPHDHeader *hdr, TriangleMesh *baseMesh, TriangleMesh *moverMesh, char **outCPHDFile,
                 int *startPulse, int *nPulses,
                 int *bounceToShow, int *nAzBeam, int *nElBeam,
                 int *interrogate, SPVector *interogPt, double *interogRad, int *interogX, int *interogY,
                 FILE **interrogateFP, int *pulseUndersampleFactor, int *polarisation, int *rayGenMethod, SPStatus *status){
    
    char * prompt;
    char * file ;
    char * interrogFname ;
    SPStatus fileStat ;
    FILE *fp;
    int createNewCPHD ;
    
    *bounceToShow           =  0 ;
    createNewCPHD           =  1 ;
    *pulseUndersampleFactor =  1 ;
    
    prompt  = (char *)sp_malloc(sizeof(char)*256);
    
    char *baseScene ;
    baseScene = tryReadFile("Name of Base scene", "baseScene",
                            "Enter the name of a file that will be the base scene to be raytraced. The file must be in a .PLY file format"
                            ,"delaunay.ply") ;
    
    // Read in the triangle mesh from the input plyfile and check it's
    // integrity
    //
    printf("Reading in file %s...\n",baseScene);
    baseMesh->readPLYFile(baseScene);
    printf("Done. Checking file Integrity...\n");
    baseMesh->checkIntegrityAndRepair();
    baseMesh->buildTriangleAABBs();
    baseMesh->buildTrianglelCentres();
    printf("Done \n");
    free(baseScene);

#ifdef DEVBUILD
    char *moversScene ;
    moversScene = input_string("Name of movers scene", "moversScene",
                              "Enter the name of a file that will contain the things that move. The file must be in a .PLY file format"
                              ,"none") ;

    if (strcmp( moversScene, "none" )) {
        moverMesh->readPLYFile(moversScene);
        moverMesh->checkIntegrityAndRepair();
        moverMesh->buildTriangleAABBs();
    }else{
        moverMesh = NULL ;
    }
    free(moversScene) ;
#else
    moverMesh = NULL ;
#endif
    
    sprintf(prompt, "cphdFile.cph");
    char *inCPHDFile ;
    inCPHDFile = tryReadFile((char *)"CPHD Filename", (char *)"CPHDFile",
                               (char *)"The name of a CPHD file to use.",
                               prompt);

    
    // Change the file extension to let user know this is a modified cphd
    //
    // find file name
    //
    strcpy(prompt, inCPHDFile) ;
    file   = strrchr(prompt, '/');
    file = (file == NULL) ? prompt : file+1;
    
    // Read in CPHD Header data to give the user a hand in choosing correct parameters
    //
    readCPHDHeader(inCPHDFile, hdr, status) ;
    
    do{
        *startPulse = 0;
        sprintf(prompt, "Start Pulse (0-%d)",hdr->num_azi-2);
        *startPulse = input_int(prompt, (char *)"startPulse", (char *)"Start pulse in CPHD file to reconstruct", *startPulse);
    } while (!(*startPulse >=0 && *startPulse <= hdr->num_azi-2)) ;
    
    sprintf(prompt, "Number of pulses (1-%d)",hdr->num_azi - *startPulse);
    *nPulses = hdr->num_azi ;
    *nPulses = input_int(prompt, (char *)"nPulses", (char *)"Number of pulses to reconstruct in cphdFile", *nPulses);
    
    if(*nPulses == 1){
        createNewCPHD = 0;
        *bounceToShow = input_int((char *)"Which bounce number to show (1-8)", (char *)"bounceToShow",
                                  (char *)"Which radar bounce should be displayed? (<1 is do not show bounce info)", *bounceToShow);
        if (*bounceToShow > MAXBOUNCES || *bounceToShow < 1) *bounceToShow = 0 ;
    }else{
        collectionGeom cGeom;
        double maxEL, maxAz, minEl, minAz, centreRange, maxBeamUsedAz, PRF, BDop, sceneAz, SADistance, azResolution, maxAzUndersamp;
        SPVector diff, Pos;
        AABB sbb;

        collectionGeometry(hdr, *startPulse + (*nPulses/2), hdr->grp, &cGeom, status);
        Pos = hdr->pulses[*startPulse + (*nPulses/2)].sat_ps_tx ;
        sceneExtent(Pos, *baseMesh, maxEL, maxAz, minEl, minAz, sbb);
        centreRange    = VECT_MAG(Pos);
        maxBeamUsedAz  = (maxAz - minAz) / centreRange ;
        PRF            = hdr->num_azi / (hdr->pulses[hdr->num_azi-1].sat_tx_time - hdr->pulses[0].sat_tx_time) ;
        BDop           = 2. * VECT_MAG(cGeom.vel) * maxBeamUsedAz * cos(cGeom.squintRad) * hdr->freq_centre / SIPC_c ;
        sceneAz        = centreRange * maxBeamUsedAz ;
        VECT_SUB(hdr->pulses[hdr->num_azi-1].sat_ps_tx, hdr->pulses[0].sat_ps_tx, diff);
        SADistance     = VECT_MAG(diff) ;
        azResolution   = fabs(cGeom.range * (SIPC_c / hdr->freq_centre) / (2.0 * SADistance * cos(cGeom.squintRad)));
        maxAzUndersamp = hdr->num_azi / (sceneAz/ azResolution) ;
        printf("Scene Doppler bandwidth     : %7.2f Hz\n",BDop);
        printf("Mean PRF is                 : %7.2f Hz\n",PRF);
        printf("Scene extent in azimuth     : %7.2f m\n",sceneAz);
        printf("Finest azimuth resolution   : %7.2f m\n",azResolution) ;
        printf("Min Az image pixels         : %5d pix \n",(int)(sceneAz / azResolution));
        printf("Max undersample for image   : %5dx\n",(int)maxAzUndersamp);
        *pulseUndersampleFactor = 1 ;
        do {
            *pulseUndersampleFactor = input_int("Pulse undersampling factor", (char *)"pulseUndersampFact",
                                                (char *)"Reduces the number of azimuth pulses that are processed. This effectively reduces the collection PRF. If the simulated scene size is small then the azimuth ambiguities will not fold in far enough to affect the simulated scene.",
                                                *pulseUndersampleFactor);
        } while (*pulseUndersampleFactor <= 0) ;
        printf("Effective PRF is            : %f Hz\n", PRF / *pulseUndersampleFactor) ;
        if (*pulseUndersampleFactor > maxAzUndersamp) {
            printf(" *** Warning : undersampling factor is larger than required to unambiguously sample the scene at the finest resolution. Ambiguities may result ! ***\n");
        }
    }
    
    *rayGenMethod = 1 ;
    do{
    printf("There are four ways to cast rays in SARCASTIC :\n");
    printf("\t1 : TRIANGLECENTRE\n");
    printf("\t2 : RANDOMRAYS\n");
    printf("\t3 : FIRSTTIMERANDOM\n");
    printf("\t4 : PARALLELRANDOM\n");
    *rayGenMethod = input_int("Enter number for ray generation method or '?' for help", "RayGenMethod",
              "SARCASTIC can generate rays in one of 4 ways. Here is an explanation of each method and why you should use it: \n\t1 : TRIANGLECENTRE \n\t\tThis method uses the centre of each triangle in the input mesh to determine the direction that\n\t\teach ray should be cast in. The origin is obviously the transmitter location for a given pulse.\n\t\tIt doesnt not try to perform any Z-buffering of triangles and so triangles that are behind another \n\t\tone will still generate a ray and will result in the nearer triangle being hit many times. (this gets \n\t\taccounted for in the RCS calculation and so doesnt affect the output.) Use this method to guarantee that\n\t\tevery triangle in the mesh (that can be illuminated by transmitted ray) is illuminated by the transmitted\n\t\tray. Use this method if: \n\t\t    You have a large scene with a smaller number of large triangles. \n\t\t    You want to make sure that every part of your model is illuminated \n\t\t    You have small triangles and so secondary and higher bounces will be incident on all faces of the mesh. \n\t2 : RANDOMRAYS \n\t\tThis method generates a Gaussian distribution of random rays. The size of the scene is measured first\n\t\tso that the entire scene is illuminated. If this method is selected then the number of rays in \n\t\tazimuth/cross-range and elevation will be asked for. Use this method if: \n\t\t    You want to perform a quick run and are not that bothered about illuminating every part of teh input mesh. \n\t\t    You have a small scene with large triangles. \n\t\t    You have large triangles and want to make sure there are many reflections from the surface of each triangle. \n\t3 : FIRSTIMERANDOM \n\t\tThis method is similar to method 2 in that it generates a Gaussian distribution of random rays. The difference\n\t\thowever is that after generating the aim point for the initial rays it then remembers the aim point for future\n\t\tpulses. This is useful if you want to make sure that pulse scattering centres are coherent from pulse to pulse.\n\t\tUse this method if: \n\t\t    You want to guarantee that a triangle correlates pulse to pulse over the entire SAR aperture \n\t4 : PARALLELRANDOM \n\t\tThis method is the same as method 2 in that it generates a Gaussian distribution of random rays. The difference\n\t\there is that the origin of each pulse is adjusted so that all the rays are parallel. Use this method if: \n\t\t    You are simulating a scene in the near field but would like it to be imaged in the far field. \n\n ",
                              *rayGenMethod);
    }while(*rayGenMethod < 1 || *rayGenMethod > 4);
    
    *nAzBeam = *nElBeam = 100 ;
    if (*rayGenMethod == RANDOMRAYS || *rayGenMethod == FIRSTTIMERANDOM || *rayGenMethod == PARALLELRANDOM) {
        *nAzBeam = input_int((char *)"Azimuth rays in radar beam?", (char *)"nAzBeam",
                             (char *)"Number of azimuth rays to use to construct radar beam. More is better but slower",*nAzBeam);
        *nElBeam = input_int((char *)"Elevation rays in radar beam?", (char *)"nElBeam",
                             (char *)"Number of elevation rays to use to construct radar beam. More is better but slower",*nElBeam);
        
    }
    
    // Find out what polarisation the sarcastic CPHD file will have
    //
    printf("Input file has polarisation set to \'%s\'\n",*(hdr->polarisation));
    char *polstr ;
    bool validpol ;
    do{
        validpol = false ;
        polstr = input_string("Enter polarisation to simulate", "Polarisation",
                              "Options are \'VV\',\'VH\',\'HV\',\'HH\',\'V_\', and \'H_\'. If one of the last two are used then the received H and V fields will be combined",*(hdr->polarisation) ) ;
        if (!strcasecmp(polstr, "vv")) {
            *polarisation = VV ;
            validpol = true ;
        }else if (!strcasecmp(polstr, "vh")){
            *polarisation = VH ;
            validpol = true ;
        }else if (!strcasecmp(polstr, "hv")){
            *polarisation = HV ;
            validpol = true ;
        }else if (!strcasecmp(polstr, "hh")){
            *polarisation = HH ;
            validpol = true ;
        }else if (!strcasecmp(polstr, "v_")){
            *polarisation = V_ ;
            validpol = true ;
        }else if (!strcasecmp(polstr, "h_")){
            *polarisation = H_ ;
            validpol = true ;
        }else{
            printf("Invalid polarisation. Options are \'VV\',\'VH\',\'HV\',\'HH\',\'V_\' or \'H_\'\n");
        }
    }while(!validpol);
    free(polstr);
    
    // To save processing time only read in the CPHDFile PHD if we are going to create a new
    // CPHD file
    //
    if( createNewCPHD ){
        // change extension
        //
        sprintf(prompt, "%s",inCPHDFile);
        char *pExt = strrchr(file, '.');
        if (pExt != NULL)
            strcpy(pExt, ".sarc.cph");
        else
            strcat(prompt, ".sarc.cph");
        
        do {
            im_init_status(fileStat, 0) ;
            *outCPHDFile = input_string((char *)"Output CPHD Filename", (char *)"CPHDOut",
                                        (char *)"The name of a CPHD file to create.",
                                        prompt);
            if ( (fp = fopen(*outCPHDFile, "w")) == NULL){
                printf(RED "Cannot access file %s\n" RESETCOLOR,*outCPHDFile);
                fileStat.status = BAD_FILE ;
            }else fclose(fp) ;
            
        } while(fileStat.status != NO_ERROR);
    }else{
        *outCPHDFile = NULL ;
    }
    
    // Ask the user if he would like to interrogate a scattering point in a previously ray traced
    // image. If so, ask for surface file and image coords of point in image so that we can convert the
    // image point to a scene coord vector. Only do this if nPulses is '1'
    //
    char *surfaceFile ;
//    int  interogX, interogY ;
    
    *interrogate = 0;
    if (*nPulses == 1) {
        
        *interrogate = input_yesno((char *)"Interrogate a point in scene?", (char *)"interogPt",
                                   (char *)"Interrogate a point in the scene to find out which scattering primitives made it.",*interrogate);
        
        if(*interrogate){
            sprintf(prompt, "surface.dat");
            do {
                im_init_status(fileStat, 0) ;
                surfaceFile = input_string((char *)"Name of surface file", (char *)"surfaceFile",
                                           (char *)"The name of a surface file to read from",
                                           prompt);
                if ( (fp = fopen(surfaceFile, "r")) == NULL){
                    printf(RED "Cannot access file %s\n" RESETCOLOR,surfaceFile);
                    fileStat.status = BAD_FILE ;
                }else fclose(fp) ;
                
            } while(fileStat.status != NO_ERROR);
            
            SPImage surface ;
            im_load(&surface, surfaceFile, &fileStat) ;
            
            *interogX  = (int)surface.nx / 2 ;
            *interogY  = (int)surface.ny / 2;
            
            *interogX  = input_int((char *)"X coordinate of image pixel to interrogate", (char *)"interogX",
                                  (char *)"location in x direction of the point in the SAR image to interrogate", *interogX) ;
            *interogY  = input_int((char *)"Y coordinate of image pixel to interrogate", (char *)"interogY",
                                  (char *)"location in y direction of the point in the SAR image to interrogate", *interogY) ;
            
            *interogRad = 1.0;
            *interogRad = input_dbl((char *)"Radius (in m)  of region around point to interrogate", (char *)"interogRad",
                                    (char *)"Radius (in m) of region around point to interrogate", *interogRad);
            
            
             do {
                sprintf(prompt, "interrogate.txt");
                
                im_init_status(fileStat, 0) ;
                interrogFname = input_string((char *)"Name of file for interrogation output", (char *)"interrogFname",
                                             (char *)"Name of file for interrogation output",
                                             prompt);
                if ( (*interrogateFP = fopen(interrogFname, "w")) == NULL){
                    printf(RED "Cannot access file %s\n" RESETCOLOR,interrogFname);
                    fileStat.status = BAD_FILE ;
                }
                
            } while(fileStat.status != NO_ERROR);
             
            *interogPt = surface.data.vect[*interogY*surface.nx+ *interogX] ;
            
            im_destroy(&surface, &fileStat);
        }
    }
    free(prompt);
    free(inCPHDFile) ;

    // Read in the material properties file if required
    //
    char *matfile = input_string((char *)"Input materialfile filename", (char *)"materialfilename",
                                 (char *)"The name of a 'materialfile' or 'none' (defaults used)",
                                 (char *) MATERIALPROPS);
    initialiseMaterials(matfile, true);
    
    return (status->status) ;
    
}
