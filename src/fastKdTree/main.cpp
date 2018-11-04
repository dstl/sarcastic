/***************************************************************************
 * 
 *           Module :  main.cpp
 *          Program :  fastKdTree
 *       Created by :  Darren Muff on 05/03/2017
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  04-Nov-2018
 *   Description:
 *      Programme to build a K-Dimensional tree quicker and more scalably than the pevious implementation
 *      This version uses the approach detailed in [1]. The previous approach uses that of [2]
 *
 *      1. Zhou, Kun, et al. "Real-time kd-tree construction on graphics hardware."
 *         ACM Transactions on Graphics (TOG) 27.5 (2008): 126.
 *
 *      2. Wald, Ingo, and Vlastimil Havran. "On building fast kd-trees for ray tracing,
 *         and on doing that in O (N log N)." Interactive Ray Tracing 2006, IEEE Symposium
 *         on. IEEE, 2006.
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

#include <iostream>
#include <sarclib/sarclib.h>
#include "colourCodes.h"
#include "buildTree.hpp"
#include "rayTrace.hpp"


#define USECOLOR 0
#define NAZRAYS 500
#define NELRAYS 500
#define NTESTRUNS 5

using namespace std;

char * tryReadFile(const char *prompt, const char * key, const char * help, const char *def);
//void indexTriangles( treeList **nodelist);
//void writeKdTreeToFile(string filename, std::vector<kdTreeNode *> *nodelist) ;


int main(int argc, const char * argv[]) {
    
    kdTree::KdData * tree;
    int treeSize;
    
    // Initialise SILib
    //
    SPStatus status ;
    im_init_status(status, 0);
    im_init_lib(&status, (char *)"fastKDTree", argc, (char **)argv);
    CHECK_STATUS_NON_PTR(status);
    
    // Read in triangle data
    //
    char *instr = tryReadFile((char *)"Input triangle filename", (char *)"infilename",
                              (char *)"The name of a triangle .ply file containing triangle data",
                              (char *) "delaunay.ply");
    

    for (int i=0; i<NTESTRUNS; ++i) {
        printf("=======================================================================\n");
        printf("\n");
        printf("                    R U N   %d\n",i);
        printf("\n");
        printf("=======================================================================\n");

        
        // Read in the triangle mesh from the input plyfile and check it's
        // integrity
        //
        TriangleMesh mesh;
        printf("Reading in file %s...\n",instr);
        mesh.readPLYFile(instr);
        printf("Done. Checking file Integrity...\n");
        mesh.checkIntegrityAndRepair();
        printf("Done \n");
        mesh.buildTriangleAABBs();
        
        // Initialise and start timer
        //
        Timer runTimer ;
        startTimer(&runTimer, &status) ;
        
        kdTree::buildTree(mesh, &tree, &treeSize, (kdTree::TREEOUTPUT)(/*kdTree::OUTPUTDATA |*/ kdTree::OUTPUTSUMM));
        printf("Done!\n");
        
        // end timer
        //
        endTimer(&runTimer, &status);
        
        
        if (USECOLOR==1) {
            printf("KdTree constructed in " BOLD BLINK GREEN " %f " RESETCOLOR "seconds \n",timeElapsedInSeconds(&runTimer, &status)) ;
        }else{
            printf("KdTree constructed in  %f seconds \n",timeElapsedInSeconds(&runTimer, &status)) ;
        }
        
        printf("Ray tracing through tree...\n");
        rayTrace(&mesh, tree, &treeSize, NAZRAYS, NELRAYS) ;
    
        free(tree) ;
        sleep(1);
        
    }
    
    im_close_lib(&status);
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
    free(prmpt) ;
    return(fname) ;
}



