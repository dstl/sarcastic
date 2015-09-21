/***************************************************************************
 *
 *       Module:    getUserInput.c
 *      Program:    SARCASTIC
 *   Created by:    Darren on 30/07/2013.
 *                  Copyright (c) 2013 Dstl. All rights reserved.
 *
 *   Description:
 *   <ENTER DESCRIPTION HERE>
 *
 *
 *   CLASSIFICATION        :  <PENDING>
 *   Date of CLASSN        :  14/03/2013
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
#include "colourCodes.h"

static char *rootpath = "/local_storage/DGM" ;

int getUserInput(char **inCPHDFile, char **KdTreeFile, char **outCPHDFile,
                 int *startPulse, int *nPulses,
                 int *bounceToShow, int *nAzBeam, int *nElBeam, int *useGPU,
                 int *interrogate, SPVector *interogPt, double *interogRad,
                 FILE **interrogateFP, SPStatus *status){
    
    char * prompt;
    char * file ;
    char * interrogFname ;
    SPStatus fileStat ;
    FILE *fp;
    int createNewCPHD ;
    
    *useGPU             =  0 ;
    *bounceToShow       =  0 ;
    createNewCPHD       =  1 ;
    
    prompt  = (char *)malloc(sizeof(char)*256);

    do {
        im_init_status(fileStat, 0) ;
        sprintf(prompt, "%s/KdTree.kdt",rootpath);
        *KdTreeFile = input_string((char *)"KdTree Filename", (char *)"KdTreeFile",
                                  (char *)"The name of a KdTree file containing the scene.",
                                  prompt);
        if( access( *KdTreeFile, F_OK ) == -1 ) {
            printf(RED "Cannot access file %s\n" RESETCOLOR,*KdTreeFile);
            fileStat.status = BAD_FILE ;
        }
    } while(fileStat.status != NO_ERROR);
    
    do {
        im_init_status(fileStat, 0) ;
        sprintf(prompt, "%s/cphdFile.cph",rootpath);
        *inCPHDFile = input_string((char *)"CPHD Filename", (char *)"CPHDFile",
                                  (char *)"The name of a CPHD file to use.",
                                  prompt);
        if( access( *inCPHDFile, R_OK ) == -1 ){
            printf(RED "Cannot access file %s\n" RESETCOLOR,*inCPHDFile);
            fileStat.status = BAD_FILE ;
        }
    } while(fileStat.status != NO_ERROR);
     
    // Change the file extension to let user know this is a modified cphd
    //
    // find file name
    //
    strcpy(prompt, *inCPHDFile) ;
    file   = strrchr(prompt, '/');
    file = (file == NULL) ? prompt : file+1;
    
    // Read in CPHD Header data to give the user a hand in choosing correct parameters
    //
    CPHDHeader hdr;
    readCPHDHeader(*inCPHDFile, &hdr, status) ;
    
    do{
        *startPulse = 0;
        sprintf(prompt, "Start Pulse (0-%d)",hdr.num_azi-2);
        *startPulse = input_int(prompt, (char *)"startPulse", (char *)"Start pulse in CPHD file to reconstruct", *startPulse);
    } while (!(*startPulse >=0 && *startPulse <= hdr.num_azi-2)) ;
    
    sprintf(prompt, "Number of pulses (1-%d)",hdr.num_azi - *startPulse);
    *nPulses = hdr.num_azi ;
    *nPulses = input_int(prompt, (char *)"nPulses", (char *)"Number of pulses to reconstruct in cphdFile", *nPulses);
    *nAzBeam = *nElBeam = 100 ;
    *nAzBeam = input_int((char *)"Azimuth rays in radar beam?", (char *)"nAzBeam",
                         (char *)"Number of azimuth rays to use to construct radar beam. More is better but slower",*nAzBeam);
    *nElBeam = input_int((char *)"Elevation rays in radar beam?", (char *)"nElBeam",
                         (char *)"Number of elevation rays to use to construct radar beam. More is better but slower",*nElBeam);
    

    if(*nPulses == 1){
        *useGPU = 0;
        createNewCPHD = 0;
        *bounceToShow = input_int((char *)"Which bounce number to show (1-8)", (char *)"bounceToShow",
                                  (char *)"Which radar bounce should be displayed? (<1 is do not show bounce info)", *bounceToShow);
        if (*bounceToShow > MAXBOUNCES || *bounceToShow < 1) *bounceToShow = 0 ;
    }
    
    // To save processing time only read in the CPHDFile PHD if we are going to create a new
    // CPHD file
    //
    if( createNewCPHD ){
        // change extension
        //
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
        *interogRad = input_dbl((char *)"Radius of region around point to interrogate", (char *)"interogRad",
                                (char *)"Radius of region around point to interrogate", *interogRad);
        
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
    
    if(*bounceToShow == 0 && *interrogate == 0)
        *useGPU = input_yesno((char *)"USE GPU?", (char *)"USEGPU",
                              (char *)"USE GPU if its available. If it isnt then use the CPU.",*useGPU);    
    
    free(prompt);

    destroy_cphd_header(&hdr, status);
    return (status->status) ;
    
}