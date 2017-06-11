/***************************************************************************
 *  main.m
 *  Sadilac
 *
 *  Created by Muff Darren on 19/05/2012.
 *  Copyright (c) 2014 [dstl]. All rights reserved.
 *
 *  CLASSIFICATION       :   UNCLASSIFIED
 *  Date of CLASSN       :   18/02/2014
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
 ***************************************************************************/

#import <Foundation/Foundation.h>
#import "MVectorObjc.h"
#import "KdTree.h"
#import "Triangle.h"
#include <SIlib2/SIlib2.h>
#include "materialProperties.h"
#define ROOTPATH "/Users/Darren/Development"
#define SHOWTRIANGLES 0
#define MAX_LINE_LENGTH 128

void banner() ;

int main (int argc, const char * argv[])
{
    SPStatus status ;
    im_init_status(status, 0);
	im_init_lib(&status, (char *)"buildKDTree", argc, (char **)argv);
	CHECK_STATUS_NON_PTR(status);

    @autoreleasepool {
        
        banner() ;
        char * instr = input_string((char *)"Triangle filename", (char *)"filename",
                                    (char *)"The name of a file containing triangle data. The file can be generated using colladaToTriFile",
                                    (char *) ROOTPATH"/DATA/triangles_Delauney.tri");
        
        FILE *fpin ;
        fpin = fopen(instr, "r");
        if (fpin == NULL) {
            printf("Error : failed to open input triangle file %s\n",instr);
            exit(1);
        }
        
         char *oustr  = input_string((char *)"KdTree Filename", (char *)"kdTreeFname",
                                     (char *)"Full path name where KdTree will be stored",
                                     (char *) ROOTPATH"/DATA/KdTree.kdt");
        
        FILE *fpout ;
        fpout = fopen(oustr, "w");
        if (fpout == NULL) {
            printf("Error : failed to open input triangle file %s\n",oustr);
            exit(1);
        }
        fclose(fpout) ;
        
        char line[MAX_LINE_LENGTH];
        char name1[MAX_LINE_LENGTH];
        char name2[MAX_LINE_LENGTH];
        char name3[MAX_LINE_LENGTH];
        char name4[MAX_LINE_LENGTH];
        char name5[MAX_LINE_LENGTH];
        char name6[MAX_LINE_LENGTH];
        char name7[MAX_LINE_LENGTH];
        char name8[MAX_LINE_LENGTH];
        int ntri=0;;
        int nvertices=0;
        int matId = 0;
        double x,y,z;
        do {
            fgets(line, MAX_LINE_LENGTH, fpin);
            sscanf(line, "%s%s%s", name1, name2, name3);
            if (strcmp(name1, "element") == 0) {
                if(strcmp(name2, "vertex") == 0){
                    nvertices = atoi(name3);
                }else if(strcmp(name2, "face") == 0){
                    ntri = atoi(name3);
                }
            }
        } while (strcmp( name1, "end_header") ) ;
        
        if (nvertices == 0 || ntri == 0) {
            printf("Error: failed to read number of triangles / vertices from .PLY header\n");
            exit(1);
        }
        
        float *vertices = (float *)sp_malloc(sizeof(float) * nvertices * 3);
        for (int i=0; i<nvertices; ++i){
            fgets(line, MAX_LINE_LENGTH, fpin);
            sscanf(line, "%s%s%s", name1, name2, name3);
            vertices[i*3+0] = atof(name1);
            vertices[i*3+1] = atof(name2);
            vertices[i*3+2] = atof(name3);
        }
        
        int *trinds = (int *)sp_malloc(sizeof(int) * ntri * 4);
        for (int i=0; i<ntri; ++i){
            fgets(line, MAX_LINE_LENGTH, fpin);
            sscanf(line, "%s%s%s%s%s%s%s%s",name1,name2, name3,name4,name5,name6,name7,name8);
            trinds[i*4+0] = atoi(name2);
            trinds[i*4+1] = atoi(name3);
            trinds[i*4+2] = atoi(name4);
            trinds[i*4+3] = atoi(name8);
        }
        
        NSMutableArray * triangles         = [[NSMutableArray alloc] init];
        for (int itri=0; itri < ntri; itri++ ) {
            x = vertices[trinds[itri*4+0]*3+0];
            y = vertices[trinds[itri*4+0]*3+1];
            z = vertices[trinds[itri*4+0]*3+2];
            MVector * AA = [MVector MVectorWithValuesX:x Y:y Z:z] ;
            x = vertices[trinds[itri*4+1]*3+0];
            y = vertices[trinds[itri*4+1]*3+1];
            z = vertices[trinds[itri*4+1]*3+2];
            MVector * BB = [MVector MVectorWithValuesX:x Y:y Z:z] ;
            x = vertices[trinds[itri*4+2]*3+0];
            y = vertices[trinds[itri*4+2]*3+1];
            z = vertices[trinds[itri*4+2]*3+2];
            MVector * CC = [MVector MVectorWithValuesX:x Y:y Z:z] ;
            matId = trinds[itri*4+3];
            Triangle * tri = [Triangle TriangleWithVerticesAa:AA Bb:BB Cc:CC andId:itri];
            [tri setMatId:matId];
            [triangles addObject:tri] ;
        }

        if (SHOWTRIANGLES) {
            printf("Triangle file contains the following triangles:\n");
            for (int itri=0; itri < ntri; itri++ ) {
                Triangle *t = [triangles objectAtIndex:itri];
                printf("Triangle %d\n", [t triId]);
                printf("  AA       : %3.6f,%3.6f,%3.6f\n", [[t vertexAt:0] X],[[t vertexAt:0] Y],[[t vertexAt:0] Z]);
                printf("  BB       : %3.6f,%3.6f,%3.6f\n", [[t vertexAt:1] X],[[t vertexAt:1] Y],[[t vertexAt:1] Z]);
                printf("  CC       : %3.6f,%3.6f,%3.6f\n", [[t vertexAt:2] X],[[t vertexAt:2] Y],[[t vertexAt:2] Z]);
                printf("  Normal   : %3.6f,%3.6f,%3.6f\n", [[t normal] X],[[t normal] Y],[[t normal] Z]);
                printf("  Material : %s \n", [[t materialName] UTF8String]);
            }
        }
        
        KdTree * kdtree = [[KdTree alloc] initWithTriangles:triangles] ;
        [kdtree printBoxes] ;
        Timer timer1,timer2 ;
        startTimer(&timer1, &status) ;
        printf("Building scene...\n");
        [kdtree Build] ;
        endTimer(&timer1, &status);
        startTimer(&timer2, &status);
        printf("Packing KDTree to file %s\n",oustr);
        [kdtree PacktoFile:[NSString stringWithUTF8String:oustr]];
        [kdtree printSummary];
        endTimer(&timer2, &status);
        printf("Scene built in %f secs\n",timeElapsedInSeconds(&timer1, &status));
        printf("Scene packed to file in %f secs\n",timeElapsedInSeconds(&timer2, &status));
        fclose( fpin ) ;

    }
    
    im_close_lib(&status);
    return 0;
}

