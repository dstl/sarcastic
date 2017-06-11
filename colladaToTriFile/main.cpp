

#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include "colladainterface.h"
#include <SIlib2/SIlib2.h>
#include "materialProperties.h"
#include "colourCodes.h"

extern "C" {
#include "matrixMultiplication.h"
}
#include "TriangleMesh.hpp"

char * tryReadFile(const char *prompt, const char * key, const char * help, const char *def) ;

std::vector<ColGeom> geom_vec;    // Vector containing COLLADA meshes
int num_objects;                  // Number of meshes in the vector

int main(int argc, char* argv[]) {
    
    SPStatus status ;
    im_init_status(status, 0) ;
    im_init_lib(&status, (char *)"readCollada", argc, (char **)argv);
    CHECK_STATUS_NON_PTR(status);
    
    char * inFile, * outFile ;
    
    inFile = tryReadFile("Collada .dae file to read", "inFile",
                         "The pathname of a Collada file containing the scene",
                         "/Users/Darren/Development/Models/scene.dae") ;
    
    outFile = input_string("Name of output file", "outFile",
                           "The pathname of the file to create. This will be a binary file containing triangle and material information of the scene",
                           "/tmp/triangles.ply");
    
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
