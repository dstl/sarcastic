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

        int ntri;
        int vertexBytes;
        int materialBytes;
        fread(&ntri, sizeof(int), 1, fpin);
        fread(&vertexBytes,sizeof(int), 1, fpin);
        fread(&materialBytes, sizeof(int), 1, fpin);
        
        if (vertexBytes != sizeof(double) ) {
            printf("Error : Triangle vertices of size %d bytes not yet supported\n",vertexBytes);
            exit(1);
        }
        if (materialBytes != MATBYTES ) {
            printf("Error : Material Bytes in file does not match the program\n");
            exit(1);
        }

        NSMutableArray * triangles         = [[NSMutableArray alloc] init];
        
        for (int itri=0; itri < ntri; itri++ ) {
            int val ;
            double AAx,AAy,AAz,BBx,BBy,BBz,CCx,CCy,CCz,NNx,NNy,NNz,area,x ;
            char mat[MATBYTES] ;
            fread(&val, sizeof(int),    1, fpin) ;
            fread(&AAx, sizeof(double), 1, fpin) ;
            fread(&AAy, sizeof(double), 1, fpin) ;
            fread(&AAz, sizeof(double), 1, fpin) ;
            MVector * AA = [MVector MVectorWithValuesX:AAx Y:AAy Z:AAz] ;
            fread(&BBx, sizeof(double), 1, fpin) ;
            fread(&BBy, sizeof(double), 1, fpin) ;
            fread(&BBz, sizeof(double), 1, fpin) ;
            MVector * BB = [MVector MVectorWithValuesX:BBx Y:BBy Z:BBz] ;
            fread(&CCx, sizeof(double), 1, fpin) ;
            fread(&CCy, sizeof(double), 1, fpin) ;
            fread(&CCz, sizeof(double), 1, fpin) ;
            MVector * CC = [MVector MVectorWithValuesX:CCx Y:CCy Z:CCz] ;
            fread(&NNx, sizeof(double), 1, fpin) ;
            fread(&NNy, sizeof(double), 1, fpin) ;
            fread(&NNz, sizeof(double), 1, fpin) ;
            MVector  * NN  = [MVector MVectorWithValuesX:NNx Y:NNy Z:NNz] ;
            Triangle * tri = [Triangle TriangleWithVerticesAa:AA Bb:BB Cc:CC andNormal:NN ];
            [tri setTriId:val] ;
            
            // Read Area - don't need it yet just read it in
            //
            fread(&area,sizeof(double), 1, fpin) ;
            
            // Same for the globalToLocalMatrix
            //
            for (int j=0; j<9; j++){
                fread(&x, sizeof(double), 1, fpin);
            }
            
            // Same for the LocalToGlobalMatrix
            //
            for (int j=0; j<9; j++){
                fread(&x, sizeof(double), 1, fpin);
            }

            fread(&mat, sizeof(char), materialBytes, fpin) ;
            NSString *strFromFile = [[NSString stringWithUTF8String: mat] uppercaseString];
            NSString *matStr ;
            for(int i=0; i<NMATERIALS; i++){
               matStr = [NSString stringWithUTF8String:materialProperties[i].matname];
                if ([strFromFile rangeOfString:matStr].location != NSNotFound) {
                    [tri setMaterialName:strFromFile] ;
                }
            }
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

