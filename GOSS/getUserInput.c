/***************************************************************************
 *
 *       Module:    getUserInput.c
 *      Program:    GOSS
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

#include "GOSS.h"
static char *rootpath = "/local_storage/DGM" ;

int getUserInput(char **inCPHDFile, char **KdTreeFile, char **outCPHDFile,
                 int *startPulse, int *nPulses,
                 int *bounceToShow, int *nAzBeam, int *nElBeam, int *useGPU,
                 SPStatus *status){
    
    char * prompt;
    char * file ;
    SPStatus fileStat ;
    FILE *fp;
    
    *useGPU             =  0 ;
    *bounceToShow       =  0 ;
    
    prompt  = (char *)malloc(sizeof(char)*256);
//    file    = (char *)malloc(sizeof(char)*256);

    do {
        im_init_status(fileStat, 0) ;
        sprintf(prompt, "%s/KdTree.kdt",rootpath);
        *KdTreeFile = input_string((char *)"KdTree Filename", (char *)"KdTreeFile",
                                  (char *)"The name of a KdTree file containing the scene.",
                                  prompt);
        if( access( *KdTreeFile, F_OK ) == -1 ) {
            printf("Cannot access file %s\n",*KdTreeFile);
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
            printf("Cannot access file %s\n",*inCPHDFile);
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
    
    // change extension
    //
    char *pExt = strrchr(file, '.');
    if (pExt != NULL)
        strcpy(pExt, ".goss.cph");
    else
        strcat(prompt, ".goss.cph");
        
    do {
        im_init_status(fileStat, 0) ;
        *outCPHDFile = input_string((char *)"Output CPHD Filename", (char *)"CPHDOut",
                               (char *)"The name of a CPHD file to create.",
                               prompt);
        if ( (fp = fopen(*outCPHDFile, "w")) == NULL){
            printf("Cannot access file %s\n",*outCPHDFile);
            fileStat.status = BAD_FILE ;
        }else fclose(fp) ;
            
    } while(fileStat.status != NO_ERROR);
    
    // Read in CPHD Header data to give the user a hand in choosing correct parameters
    //
    CPHDHeader hdr;
    readCPHDHeader(*inCPHDFile, &hdr, status) ;
    
    *startPulse = 0;
    sprintf(prompt, "Start Pulse (0-%d)",hdr.num_azi);
    *startPulse = input_int(prompt, (char *)"startPulse", (char *)"Start pulse in CPHD file to reconstruct", *startPulse);

    sprintf(prompt, "Number of pulses (1-%d)",hdr.num_azi - *startPulse);
    *nPulses = hdr.num_azi ;
    *nPulses = input_int(prompt, (char *)"nPulses", (char *)"Number of pulses to reconstruct in cphdFile", *nPulses);
    if(*nPulses == 1){
        *useGPU = 0;
        *bounceToShow = input_int((char *)"Which bounce number to show (1-8)", (char *)"bounceToShow",
                                  (char *)"Which radar bounce should be displayed? (<1 is do not show bounce info)", *bounceToShow);
        if (*bounceToShow > MAXBOUNCES || *bounceToShow < 1) *bounceToShow = 0 ;
    }

    *nAzBeam = *nElBeam = 100 ;
    *nAzBeam = input_int((char *)"Azimuth rays in radar beam?", (char *)"nAzBeam",
                        (char *)"Number of azimuth rays to use to construct radar beam. More is better but slower",*nAzBeam);
    *nElBeam = input_int((char *)"Elevation rays in radar beam?", (char *)"nElBeam",
                        (char *)"Number of elevation rays to use to construct radar beam. More is better but slower",*nElBeam);
    if(*bounceToShow == 0)
        *useGPU = input_yesno((char *)"USE GPU?", (char *)"USEGPU",
                             (char *)"USE GPU if its available. If it isnt then use the CPU.",*useGPU);
    
    free(prompt);
//    free(file) ;

    destroy_cphd_header(&hdr, status);
    return (status->status) ;
    
}