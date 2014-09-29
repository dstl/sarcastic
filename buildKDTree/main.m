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
#import "COLScene.h"
#import "KdTree.h"
#import "readCollada.h"
#import "Timer.h"
#include <mtlib/mtlib.h>
#define ROOTPATH "/Users/Darren/Development"
// #include "Sadilac.h"

void banner() ;

int main (int argc, const char * argv[])
{
    SPStatus mtlibstatus ;
    char * str;
    mtlibstatus.debug = 0 ;
	im_init_lib(&mtlibstatus, (char *)"Sadilac", argc, (char **)argv);
	CHECK_STATUS_NON_PTR(mtlibstatus);

    @autoreleasepool {
        
        banner() ;
        str =input_string((char *)"Model filename", (char *)"filename",
                                 (char *)"The name of a model in the collada DAE format. Export from Sketchup(R)",
                                 (char *) ROOTPATH"/Models/model.dae");
        
        NSString *filestr = [NSString stringWithUTF8String:str];
        if ([filestr rangeOfString:@"file://"].location == NSNotFound) {
            NSMutableString *newFileStr = [NSMutableString stringWithUTF8String:"file://"] ;
            filestr = [newFileStr stringByAppendingString:filestr];
        }
        NSURL * ipurl = [NSURL URLWithString:filestr];
        
        str =input_string((char *)"KdTree Filename", (char *)"kdTreeFname",
                          (char *)"Full path name where KdTree will be stored",
                          (char *) ROOTPATH"/DATA/KdTree.kdt");
        
        NSString * treeFname = [NSString stringWithUTF8String:str] ;

        
        Timer * timer = [[Timer alloc] init];
        [timer startTimer];
        COLScene * scene = [COLScene sceneWithURL:ipurl];
//        [scene printTree];
        NSLog(@"Scaling is %f\n",[scene sceneScaling]);
        NSLog(@"Scene has %lu triangles\n",(unsigned long)[[scene triangles] count]);
        KdTree * kdtree = [[KdTree alloc] initWithScene:scene];
        
        NSLog(@"Building Scene...");
        Timer * timer0 = [[Timer alloc] init];
        [timer0 startTimer];
        [kdtree Build];
        [timer0 stopTimer];
        NSLog(@"KdTree built in %f seconds",[timer0 timeElapsedInSeconds]);
        
        NSLog(@"Packing Kdtree...");
        Timer * timer1 = [[Timer alloc] init];
        [timer1 startTimer];
        [kdtree PacktoFile:treeFname];
        [timer1 stopTimer];
        NSLog(@"Packing completed in %f seconds",[timer1 timeElapsedInSeconds]);

        [timer stopTimer];
        [kdtree printSummary];
        printf("Completed in %f seconds\n",[timer timeElapsedInSeconds]);
        
    }
    
    im_close_lib(&mtlibstatus);
    return 0;
}

