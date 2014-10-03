//
//  main.c
//  sarcastic
//
//  Created by Darren Muff on 03/10/2014.
//  Copyright (c) 2014 Dstl. All rights reserved.
//

#include <stdio.h>
#include <SIlib.h>
#include "POTriangle.h"
#include "matrixMultiplication.h"

void buildTriangle(SPVector AA, SPVector BB, SPVector CC, triangle * tri) ;

int main(int argc, const char * argv[])
{
    
    triangle tri ;
    SPVector AA,BB,CC ;
    VECT_CREATE( 10,   0,  0, AA);
    VECT_CREATE(  0,  10,  0, BB);
    VECT_CREATE(-10, -10,  0, CC);
    buildTriangle(AA, BB, CC, &tri) ;
    
    // Create a ray
    //
    Ray r ;
    SPVector hp, dir;
    VECT_CREATE(20.0, 20.0, 20.0, r.org);
    VECT_CREATE(0.0, 0.0, 0.0, hp);
    VECT_SUB(hp, r.org, dir);
    r.len = 0.0;
    VECT_NORM(dir, r.dir);
    r.pow = 1000000.0;
    SPVector zhat, Hpol, Vpol ;
    VECT_CREATE(0, 0, 1, zhat);
    VECT_CROSS(r.dir, zhat, Hpol) ;
    VECT_CROSS(Hpol, r.dir, Vpol) ;
    r.pol = Vpol ;
    
    printf("Transmit location %f,%f,%f\n", r.org.x,r.org.y,r.org.z);
    printf("Hipoint location  %f,%f,%f\n", hp.x,hp.y,hp.z);
    printf(" Facet:\n");
    printf("    %f%f%f\n",AA.x,AA.y,AA.z);
    printf("    %f%f%f\n",BB.x,BB.y,BB.z);
    printf("    %f%f%f\n",CC.x,CC.y,CC.z);
    printf("    %f%f%f\n",AA.x,AA.y,AA.z);

    POTriangle(tri, ray) ;
    
    return 0;
}


void buildTriangle(SPVector AA, SPVector BB, SPVector CC, triangle * tri){
    

    SPVector NN, l1,l2, zhat ;
    double alpha, beta;
    double T_dash[9];
    double T_dashdash[9];
    
    tri->id = 0;
    tri->AA = AA ;
    tri->BB = BB ;
    tri->CC = CC ;
    VECT_SUB(BB, AA, l1) ;
    VECT_SUB(CC, AA, l2) ;
    VECT_CROSS(l1, l2, NN) ;
    tri->area = VECT_MAG(NN) * 0.5 ;
    VECT_NORM(NN, tri->NN) ;
    strcpy (tri->mat, "MATERIAL") ;
    
    VECT_CREATE(0, 0, 1, zhat);
    alpha = atan2(NN.y, NN.x);
    beta  = acos(VECT_DOT(zhat, NN));
    T_dash[0] = cos(alpha);
    T_dash[1] = sin(alpha);
    T_dash[2] = 0;
    T_dash[3] = -sin(alpha);
    T_dash[4] = cos(alpha);
    T_dash[5] = 0;
    T_dash[6] = 0;
    T_dash[7] = 0;
    T_dash[8] = 1;
    T_dashdash[0] = cos(beta);
    T_dashdash[1] = 0;
    T_dashdash[2] = -sin(beta);
    T_dashdash[3] = 0;
    T_dashdash[4] = 1;
    T_dashdash[5] = 0;
    T_dashdash[6] = sin(beta);
    T_dashdash[7] = 0;
    T_dashdash[8] = cos(beta);
    
    
    
    matmul(T_dashdash, T_dash, tri->globalToLocalMat, 3, 3, 3, 3);
    mat3by3inv(tri->globalToLocalMat, tri->localToGlobalMat);

    return ;
}