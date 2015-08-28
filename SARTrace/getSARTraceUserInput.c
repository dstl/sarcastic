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

#include "SARTrace.h"
#include "colourCodes.h"

static char *rootpath = "/local_storage/DGM" ;

int getSARTraceUserInput(char **inCPHDFile, char **KdTreeFile, char **outDir, int *nRaysX, int *nRaysY, SPStatus *status){
    
    char * prompt;
    SPStatus fileStat ;
    
    prompt  = (char *)malloc(sizeof(char)*256);
    
    do {
        im_init_status(fileStat, 0) ;
        sprintf(prompt, "%s/KdTree.kdt",rootpath);
        *KdTreeFile = input_string((char *)"KdTree Filename", (char *)"KdTreeFile",
                                   (char *)"The name of a KdTree file containing the scene. The KdTree is a hierarchical structure containing all the scene information and can be created from a '.dae' file using the program 'buildKdTree'",
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
                                   (char *)"The name of a CPHD file to use. The CPHD file is used to position the radar sensor location and RF parameters relative to the scene",
                                   prompt);
        if( access( *inCPHDFile, R_OK ) == -1 ){
            printf(RED "Cannot access file %s\n" RESETCOLOR,*inCPHDFile);
            fileStat.status = BAD_FILE ;
        }
    } while(fileStat.status != NO_ERROR);
    
    do {
        im_init_status(fileStat, 0) ;
        sprintf(prompt, "%s/rayTraceOutput",rootpath);
        *outDir = input_string((char *)"Output directory", (char *)"OutDir",
                               (char *)"The name of a directory / folder to place bounce ray intersection information. Two sets of files are created in this folder, both in ',ply' format. The first contains the locations of each ray intersection with different files showing the locations after each bounce. The second contains the actual triangles that are hit on each bounce. This is useful as the Physical Optics properties of each triangle are calculated for the whole triangle and so only one hit per triangle is required. You can therefore use this to determine how many rays are required to fully illuminate the scene.",
                               prompt);
        if( access( *outDir, R_OK ) == -1 ){
            printf(RED "Cannot access file %s\n" RESETCOLOR,*outDir);
            fileStat.status = BAD_FILE ;
        }
    } while(fileStat.status != NO_ERROR);
    
    *nRaysX = 1024 ;
    *nRaysX = input_int((char *)"Number of rays in X", (char *)"nraysx",
                        (char *)"The number of rays to cast in the X or horizontal direction. A radar beam is simulated by generating a 2D gaussian random distribution of rays that illuminate the entire scene. Larger scenes therefore require more rays to fully sample the scene", *nRaysX);
    
    *nRaysY = 1024 ;
    *nRaysY = input_int((char *)"Number of rays in Y", (char *)"nraysy",
                        (char *)"The number of rays to cast in the Y or vertical direction. A radar beam is simulated by generating a 2D gaussian random distribution of rays that illuminate the entire scene. Larger scenes therefore require more rays to fully sample the scene", *nRaysY);
    
    free(prompt);
    
    return (status->status) ;
    
}
