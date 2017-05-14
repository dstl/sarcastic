//
//  threadCore.hpp
//  sarcastic
//
//  Created by Darren Muff on 17/03/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#ifndef threadCore_hpp
#define threadCore_hpp

#include <stdio.h>
#include <SIlib2/SIlib2.h>
#include "AABB.hpp"
extern "C" {
#include "OpenCLUtils.h"
}
#include "kdTreeNode.hpp"
#include "accelerateTriangles.hpp"
#include "rayTrace.hpp"
#include "reflect.hpp"
#include "shadowRays.hpp"
#include "buildRays.hpp"

#define NPOINTS (32)
#define OVERSAMP (512)
#define MAXBOUNCES  10  // Maximum number of ray bounces

typedef struct SPCmplxD{
    double r ;
    double i ;
} SPCmplxD;

typedef struct rangeAndPower {
    double  range ;
    SPCmplxD Es ;
} rangeAndPower ;

typedef struct rnpData_t {
    SPCmplxD Es;
    double  samplingOffset;
    int     samplingOffsetInt;
    double  indexOffset ;
    double  rdiff;
}rnpData_t ;

typedef struct threadData {
    int tid;
    int nThreads ;
    AABB SceneBoundingBox ;
    int startPulse ;
    int nPulses ;
    int nAzBeam ;
    int nElBeam ;
    CPHDHeader * cphdhdr ;
    double gainRx ;
    SPImage * phd ;             // Pointer to beginning of cphd data
    int bounceToShow ;
    SPStatus status ;
    double PowPerRay ;          // Ray Power (Pp = (Pt * Gtx ))
    int interrogate ;           // Do we want to write out details about an interrogation point in the scene?
    SPVector interogPt ;        // Position in scene coordinates of a point to be interogated
    double interogRad ;         // Radius in metres around interrogation point to calculate scattering for
    FILE ** interogFP ;         // File pointer to dump out interrogation data
    TriangleMesh *sceneMesh ;    // Mesh of triangles in scene
    TriangleMesh *moverMesh ;
} threadData ;

void * devPulseBlock ( void * threadArg ) ;


void packSinc(SPCmplxD point, SPCmplx *outData, double rdiff, double sampleSpacing, long long nxInData, double * ikernel);
void ham1dx(double * data, int nx) ;
void sinc_kernel(int oversample_factor, int num_of_points, double resolution, double sampleSpacing, double *ikernel);
void generate_gaussian_ray_distribution(double xStdev,double yStdev, SPVector txPos, SPVector GRP, int nx,int ny, Ray *rayArray);


#endif /* threadCore_hpp */
