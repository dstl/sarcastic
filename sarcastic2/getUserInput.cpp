//
//  getUserInput.cpp
//  sarcastic
//
//  Created by Darren Muff on 17/03/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#include "getUserInput.hpp"
#include "colourCodes.h"
#include "tryReadFile.hpp"

#define ROOTPATH "/tmp"

int getUserInput(CPHDHeader *hdr, TriangleMesh *baseMesh, TriangleMesh *moverMesh, char **outCPHDFile,
                 int *startPulse, int *nPulses,
                 int *bounceToShow, int *nAzBeam, int *nElBeam,
                 int *interrogate, SPVector *interogPt, double *interogRad,
                 FILE **interrogateFP, int * pulseUndersampleFactor, SPStatus *status){
    
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
                            ,ROOTPATH"/delaunay.ply") ;
    
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
    
    
    sprintf(prompt, "%s/cphdFile.cph",ROOTPATH);
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
        
        do {
            *pulseUndersampleFactor = input_int("Pulse undersampling factor", (char *)"pulseUndersampFact",
                                                (char *)"Reduces the number of azimuth pulses that are processed. This effectively reduces teh collection PRF. If the simulated scene size is small then the azimuth ambiguities will not fold in far enough to affect the simulated scene.",
                                                *pulseUndersampleFactor);
        } while (*pulseUndersampleFactor <= 0) ;
        
    }
    
    *nAzBeam = *nElBeam = 100 ;
    *nAzBeam = input_int((char *)"Azimuth rays in radar beam?", (char *)"nAzBeam",
                         (char *)"Number of azimuth rays to use to construct radar beam. More is better but slower",*nAzBeam);
    *nElBeam = input_int((char *)"Elevation rays in radar beam?", (char *)"nElBeam",
                         (char *)"Number of elevation rays to use to construct radar beam. More is better but slower",*nElBeam);
    
    
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
    // image point to a scene coord vector.
    //
    char *surfaceFile ;
    int  interogX, interogY ;
    
    *interrogate = 0;
    *interrogate = input_yesno((char *)"Interrogate a point in scene?", (char *)"interogPt",
                               (char *)"Interrogate a point in the scene to find out which scattering primitives made it.",*interrogate);
    
    if(*interrogate){
        sprintf(prompt, "/local_storage/DGM/surface.dat");
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
        
        interogX  = (int)surface.nx / 2 ;
        interogY  = (int)surface.ny / 2;
        
        interogX  = input_int((char *)"X coordinate of image pixel to interrogate", (char *)"interogX",
                              (char *)"location in x direction of the point in the SAR image to interrogate", interogX) ;
        interogY  = input_int((char *)"Y coordinate of image pixel to interrogate", (char *)"interogY",
                              (char *)"location in y direction of the point in the SAR image to interrogate", interogY) ;
        
        *interogRad = 1.0;
        *interogRad = input_dbl((char *)"Radius (in m)  of region around point to interrogate", (char *)"interogRad",
                                (char *)"Radius (in m) of region around point to interrogate", *interogRad);
        
        do {
            sprintf(prompt, "/local_storage/DGM/interrogate.txt");
            
            im_init_status(fileStat, 0) ;
            interrogFname = input_string((char *)"Name of file for interrogation output", (char *)"interrogFname",
                                         (char *)"Name of file for interrogation output",
                                         prompt);
            if ( (*interrogateFP = fopen(interrogFname, "w")) == NULL){
                printf(RED "Cannot access file %s\n" RESETCOLOR,interrogFname);
                fileStat.status = BAD_FILE ;
            }
            
        } while(fileStat.status != NO_ERROR);
        
        *interogPt = surface.data.vect[interogY*surface.nx+interogX] ;
        
        im_destroy(&surface, &fileStat);
    }
    
    free(prompt);

    return (status->status) ;
    
}
