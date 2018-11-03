/***************************************************************************
 * 
 *           Module :  main.cpp
 *          Program :  colladaToPlyFile
 *       Created by :  Darren Muff on Wed Jul 11 08:12:36 2018 -0400
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *      Description :  Program to read in a collada (.dae) file and write 
 *                     it out as a Stanford Polygon File (.ply)
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

#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include "colladainterface.h"
#include <sarclib/sarclib.h>
#include "colourCodes.h"
#include "readMaterialFile.hpp"

extern "C" {
#include "matrixMultiplication.h"
}
#include "TriangleMesh.hpp"

void banner ();
char * tryReadFile(const char *prompt, const char * key, const char * help, const char *def) ;

std::vector<ColGeom> geom_vec;    // Vector containing COLLADA meshes
int num_objects;                  // Number of meshes in the vector

int main(int argc, char* argv[]) {
    
    SPStatus status ;
    im_init_status(status, 0) ;
    im_init_lib(&status, (char *)"colladaToPlyFile", argc, (char **)argv);
    CHECK_STATUS_NON_PTR(status);
    
    banner(); 
    char * inFile, * outFile ;
    
    inFile = tryReadFile("Collada .dae file to read", "inFile",
                         "The pathname of a Collada file containing the scene",
                         "scene.dae") ;
    
    outFile = input_string("Name of output file", "outFile",
                           "The pathname of the file to create. This will be a binary file containing triangle and material information of the scene",
                           "triangles.ply");
    
    // Read in the material properties file if required
    //
    char *matfile = input_string((char *)"Input materialfile filename", (char *)"materialfilename",
                                 (char *)"The name of a 'materialfile' or 'none' (defaults used)",
                                 (char *) MATERIALPROPS);
    initialiseMaterials(matfile, true);
    
    TriangleMesh mesh;
    mesh.readDAEFile(inFile) ;
    mesh.writePLYFile(outFile);
    
    im_close_lib(&status);

    std::cout << "Done" << std::endl;
    
    return 0;
}

char * tryReadFile(const char *prompt, const char * key, const char * help, const char *def)
    ///  prompt  : the text to display when asking the user for input (the default value will also be displayed)
    ///  key     : a unique bit of text (such az "AzPixelSize") which the input system will then use to store parameter values
    ///  help    : the text to display to the user when they just enter '?'
    ///  def     : the default, ie the value to take if the user just presses return
    ///
{
    FILE *fp;
    SPStatus fileStat ;
    char *prmpt, *fname;

    size_t len = strlen(def);
    prmpt = (char *)sp_malloc(sizeof(char) * len);
    
    do {
        im_init_status(fileStat, 0) ;
        strcpy(prmpt, def);
        fname = input_string(prompt, key, help, def);
        
        if ( (fp = fopen(fname, "r")) == NULL){
            printf(RED "Cannot access file %s\n" RESETCOLOR,fname);
            fileStat.status = BAD_FILE ;
        }
        
    } while(fileStat.status != NO_ERROR);
    
    fclose(fp);
    
    return(fname) ;
}
