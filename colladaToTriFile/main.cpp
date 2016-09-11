

#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include "colladainterface.h"
#include <SIlib2/SIlib2.h>
#include "trianglecpp.h"
#include "materialProperties.h"
#include "colourCodes.h"

extern "C" {
#include "matrixMultiplication.h"
}
#include "TriangleFile.h"

char * tryReadFile(const char *prompt, const char * key, const char * help, const char *def) ;

std::vector<ColGeom> geom_vec;    // Vector containing COLLADA meshes
int num_objects;                  // Number of meshes in the vector

int main(int argc, char* argv[]) {
    
    SPStatus status ;
    im_init_status(status, 0) ;
    im_init_lib(&status, (char *)"readCollada", argc, (char **)argv);
    CHECK_STATUS_NON_PTR(status);
    
    char * inFile, * outFile ;
    bool plyop=false;
    char *plyFname=NULL;
    
    inFile = tryReadFile("Collada .dae file to read", "inFile",
                         "The pathname of a Collada file containing the scene",
                         "/Users/Darren/Development/Models/scene.dae") ;
    
    outFile = input_string("Name of output file", "outFile",
                           "The pathname of the file to create. This will be a binary file containing triangle and material information of the scene",
                           "/tmp/triangles.tri");
    
    plyop = input_yesno("Write out collada file as PLY file", "plyop",
                        "The PLY format is a simple-to-understand representation ofteh transgles in the Collada scene",
                        plyop);
    if(plyop){
        plyFname = input_string("name of PLY file to create", "plyFname",
                               "Full path name of the PLY file to create", plyFname);
    }
    
    
    // Geom_vec is a C++ Vector containing an array of 'meshes' from the collada file
    // We are only interested in the triangles in the file so loop though the
    // Geom_vec array, find all the triangles then write them to a triangle vector
    //
    ColladaInterface::readGeometries(&geom_vec, inFile);
    
    std::vector<Triangle> tri_vec;
    num_objects = (int) geom_vec.size();
    for(int i=0; i<num_objects; i++){
        if(geom_vec[i].primitive == GL_TRIANGLES ){
            float *f      = (float *)geom_vec[i].map["POSITION"].data;
            int stride    = geom_vec[i].map["POSITION"].stride ;
            float scaling = geom_vec[i].scaling;
            double localPt[4], globalPt[4];

            for (int j=0; j<geom_vec[i].index_count/3; j++){
                short AAi,BBi,CCi;
                AAi = geom_vec[i].indices[j*3+0] ;
                BBi = geom_vec[i].indices[j*3+1] ;
                CCi = geom_vec[i].indices[j*3+2] ;
                float AAx,AAy,AAz, BBx,BBy,BBz, CCx,CCy,CCz;
                
                localPt[0] = f[AAi*stride+0] * scaling;
                localPt[1] = f[AAi*stride+1] * scaling;
                localPt[2] = f[AAi*stride+2] * scaling;
                localPt[3] = scaling ;
                matmul(geom_vec[i].transform, localPt, globalPt, 4, 4, 1, 4);
                AAx = globalPt[0];
                AAy = globalPt[1];
                AAz = globalPt[2];
                
                localPt[0] = f[BBi*stride+0] * scaling;
                localPt[1] = f[BBi*stride+1] * scaling;
                localPt[2] = f[BBi*stride+2] * scaling;
                localPt[3] = scaling ;
                matmul(geom_vec[i].transform, localPt, globalPt, 4, 4, 1, 4);
                BBx = globalPt[0];
                BBy = globalPt[1];
                BBz = globalPt[2];

                localPt[0] = f[CCi*stride+0] * scaling;
                localPt[1] = f[CCi*stride+1] * scaling;
                localPt[2] = f[CCi*stride+2] * scaling;
                localPt[3] = scaling ;
                matmul(geom_vec[i].transform, localPt, globalPt, 4, 4, 1, 4);
                CCx = globalPt[0];
                CCy = globalPt[1];
                CCz = globalPt[2];

                SPVector AA,BB,CC;
                VECT_CREATE(AAx, AAy, AAz, AA);
                VECT_CREATE(BBx, BBy, BBz, BB);
                VECT_CREATE(CCx, CCy, CCz, CC);
                
                std::string triMaterial ;
                if (geom_vec[i].materialSide1 == std::string("material") || geom_vec[i].materialSide1 == std::string("")) {
                    if (geom_vec[i].materialSide2 != std::string("material") && geom_vec[i].materialSide2 != std::string("")) {
                        triMaterial = geom_vec[i].materialSide2;
                    }else{
                        triMaterial = materialProperties[0].matname;
                    }
                }else if (geom_vec[i].materialSide2 == std::string("material") || geom_vec[i].materialSide2 == std::string("")){
                    if (geom_vec[i].materialSide1 != std::string("material") && geom_vec[i].materialSide1 != std::string("")) {
                        triMaterial = geom_vec[i].materialSide1;
                    }else{
                        triMaterial = materialProperties[0].matname;
                    }
                }
                
                Triangle t = Triangle(AA, BB, CC,triMaterial);
                if(t.area != 0.0){
                    printf("A: %5.2f,%5.2f,%5.2f B: %5.2f,%5.2f,%5.2f C: %5.2f,%5.2f.%5.2f N: %5.2f,%5.2f.%5.2f - %s\n",AA.x,AA.y,AA.z,BB.x,BB.y,BB.z,CC.x,CC.y,CC.z, t.NN.x,t.NN.y,t.NN.z,triMaterial.c_str());
                    tri_vec.push_back(t);
                }else{
                    printf("ERROR area is %f: A: %f,%f,%f  B: %f,%f,%f  C: %f,%f.%f - %s\n",t.area,AA.x,AA.y,AA.z,BB.x,BB.y,BB.z,CC.x,CC.y,CC.z,triMaterial.c_str());

                }
            }
        }
    }
    
    TriangleFile trifile(tri_vec);
    trifile.WriteFile(outFile);
    
    if(plyop){
        trifile.WritePLYFile(plyFname, false);
    }
    
    
    ColladaInterface::freeGeometries(&geom_vec);
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
