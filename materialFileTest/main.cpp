//
//  main.cpp
//  materialFileTest
//
//  Created by Darren Muff on 08/09/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#include <iostream>
#include <SIlib2/SIlib2.h>
#include "readMaterialFile.hpp"
#include "materialProperties.h"

using namespace std;
char * tryReadFile(const char *prompt, const char * key, const char * help, const char *def);

int main(int argc, const char * argv[]) {
    SPStatus status ;
    im_init_status(status, 0);
    im_init_lib(&status, "materialFileTest", argc, (char **)argv);
    CHECK_STATUS_NON_PTR(status) ;
    
    char *instr = input_string((char *)"Input materialfile filename", (char *)"infilename",
                              (char *)"The name of a materialfile'",
                              (char *) "materialProperties.txt");
    
    if(string(instr) == string("none") || string(instr) == string("None") || string(instr) == string("")){
        globalNumMats    = NMATERIALS ;
        globalMatProps   = new scatProps [globalNumMats] ;
        globalMatColours = new int [globalNumMats * 3] ;
        for(int i=0; i<globalNumMats; ++i){
            globalMatProps[i]  = materialProperties[i] ;
            globalMatColours[i*3+0] = materialColours[i][0] ;
            globalMatColours[i*3+1] = materialColours[i][1] ;
            globalMatColours[i*3+2] = materialColours[i][2] ;
        }
        
    }else{
        readMaterialFile(string(instr), true) ;
    }
    
    printf("global N materials is %d\n", globalNumMats) ;
    
    delete [] globalMatProps ;
    delete [] globalMatColours ;
    free(instr);
    
    im_close_lib(&status);

    return 0;
}

