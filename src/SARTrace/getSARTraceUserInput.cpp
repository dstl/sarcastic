/***************************************************************************
 * 
 *           Module :  getSARTraceUserInput.cpp
 *          Program :  SARTrace
 *       Created by :  Darren Muff on 30/07/2013
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *           Function to ask for user input before executing the
 *           rest of the program
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

#include "SARTrace.hpp"
#include "colourCodes.h"
#include "readMaterialFile.hpp"

int getSARTraceUserInput(char **inCPHDFile, char **meshFile, char **outDir, int *pulseToTrace, int *nRaysX, int *nRaysY,  SPStatus *status){
    
    char * prompt;
    SPStatus fileStat ;
    
    prompt  = (char *)malloc(sizeof(char)*256);
    
    do {
        im_init_status(fileStat, 0) ;
        sprintf(prompt, "delaunay.ply");
        *meshFile = input_string((char *)"Triangle mesh file", (char *)"plyfile",
                   (char *)"The name of a triangle mesh \'.ply\' file containing the scene. The file can be created using colldaToTriFile using a .dae file as input (you can create a .dae file from Google Sketchup (tm))",
                   prompt);
        if( access( *meshFile, F_OK ) == -1 ) {
            printf(RED "Cannot access file %s\n" RESETCOLOR,*meshFile);
            fileStat.status = BAD_FILE ;
        }
    } while(fileStat.status != NO_ERROR);
    
    do {
        im_init_status(fileStat, 0) ;
        sprintf(prompt, "cphdFile.cph");
        *inCPHDFile = input_string((char *)"CPHD Filename", (char *)"CPHDFile",
                                   (char *)"The name of a CPHD file to use. The CPHD file is used to position the radar sensor location and RF parameters relative to the scene",
                                   prompt);
        if( access( *inCPHDFile, R_OK ) == -1 ){
            printf(RED "Cannot access file %s\n" RESETCOLOR,*inCPHDFile);
            fileStat.status = BAD_FILE ;
        }
    } while(fileStat.status != NO_ERROR);
    
    CPHDHeader hdr;
    readCPHDHeader(*inCPHDFile, &hdr, status) ;
    CHECK_STATUS(status) ;
    printf("CPHDFile has %d pulses\n", hdr.num_azi);
    *pulseToTrace = hdr.num_azi / 2 ;
    do {
        *pulseToTrace = input_int("Whioh pulse would you like to trace?", "PULSETOTRACE", "Must be a pulse index of a valid pulse within the cphd dataset", *pulseToTrace) ;
    } while (*pulseToTrace < 0 || *pulseToTrace > hdr.num_azi-1) ;
    
    do {
        im_init_status(fileStat, 0) ;
        sprintf(prompt, "rayTraceOutput");
        *outDir = input_string((char *)"Output directory", (char *)"OutDir",
                               (char *)"The name of a directory / folder to place bounce ray intersection information. Two sets of files are created in this folder, both in ',ply' format. The first contains the locations of each ray intersection with different files showing the locations after each bounce. The second contains the actual triangles that are hit on each bounce. This is useful as the Physical Optics properties of each triangle are calculated for the whole triangle and so only one hit per triangle is required. You can therefore use this to determine how many rays are required to fully illuminate the scene.",
                               prompt);
        if( access( *outDir, R_OK ) == -1 ){
            printf(RED "Cannot access directory %s ! \n" RESETCOLOR,*outDir);
            printf(GREEN "I'll create it for you..." RESETCOLOR);
            const int dir_err = mkdir(*outDir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            if (-1 == dir_err){
                printf(RED "Error creating directory! \n" RESETCOLOR);
                fileStat.status = BAD_FILE ;
            }else{
                printf(GREEN "Directory successfully created !\n" RESETCOLOR);
            }
        }
    } while(fileStat.status != NO_ERROR);
    
    *nRaysX = 1024 ;
    *nRaysX = input_int((char *)"Number of rays in X", (char *)"nraysx",
                        (char *)"The number of rays to cast in the X or horizontal direction. A radar beam is simulated by generating a 2D gaussian random distribution of rays that illuminate the entire scene. Larger scenes therefore require more rays to fully sample the scene", *nRaysX);
    
    *nRaysY = 1024 ;
    *nRaysY = input_int((char *)"Number of rays in Y", (char *)"nraysy",
                        (char *)"The number of rays to cast in the Y or vertical direction. A radar beam is simulated by generating a 2D gaussian random distribution of rays that illuminate the entire scene. Larger scenes therefore require more rays to fully sample the scene", *nRaysY);
    
    free(prompt);
    
    // Read in the material properties file if required
    //
    char *matfile = input_string((char *)"Input materialfile filename", (char *)"materialfilename",
                                 (char *)"The name of a 'materialfile' or 'none' (defaults used)",
                                 (char *) MATERIALPROPS);
    initialiseMaterials(matfile, true);
    
    return (status->status) ;
    
}
