//
//  main.c
//  POTriangleTest
//
//  Created by Darren on 30/08/2014.
//  Copyright (c) 2014 Dstl. All rights reserved.
//

#include <stdio.h>
#include <SIlib.h>
#include "sarcastic.h"

int test0();

int main(int argc, const char * argv[])
{
    int ret;
    
    ret = test0();
    
    if (ret != 0) {
        printf("Test Failed \n");
    }else{
        printf("Test Passed \n");
    }
    return ret;
}

int test0(){
    
    // Set up a triangle
    //
    SPVector AA,BB,CC ;
    VECT_CREATE(10, 0, 0, AA);
    VECT_CREATE(0, 10, 0, BB);
    VECT_CREATE(-10, -10, 0, CC);
    TriCoords triCoords ;
    triCoords.A=AA ;
    triCoords.B=BB ;
    triCoords.Cc=CC ;
    
    // Create a ray
    //
    Ray r ;
    SPVector hp, dir;
    VECT_CREATE(20.0, 20.0, 20.0, r.org);
    VECT_CREATE(0.0, 0.0, 0.0, hp);
    VECT_SUB(hp, r.org, dir);
    r.len = VECT_MAG(dir);
    VECT_NORM(dir, r.dir);
    r.pow = 1000000.0;
    
    printf("Transmit location %f,%f,%f\n", r.org.x,r.org.y,r.org.z);
    printf("Hipoint location  %f,%f,%f\n", hp.x,hp.y,hp.z);
    printf(" Facet:\n");
    printf("    %f%f%f\n",AA.x,AA.y,AA.z);
    printf("    %f%f%f\n",BB.x,BB.y,BB.z);
    printf("    %f%f%f\n",CC.x,CC.y,CC.z);
    printf("    %f%f%f\n",AA.x,AA.y,AA.z);
    
    return 0;
    
}
