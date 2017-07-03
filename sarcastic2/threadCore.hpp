/***************************************************************************
 *
 *       Module:    threadCore.hpp
 *      Program:    SARCASTIC
 *   Created by:    Darren on 15/03/2017.
 *                  Copyright (c) 2017 Dstl. All rights reserved.
 *
 *   Description:
 *      Core routine to be run on each thread
 *
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  15/03/2017
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
 *
 ***************************************************************************/

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
    kdTree::KdData ** tree ;
    int treesize; 
    ATS **accelTriangles ;
    int polarisation ;
} threadData ;

enum  POLARISATION { VV, VH, HV, HH, V_, H_ } ;

void * devPulseBlock ( void * threadArg ) ;
void packSinc(SPCmplxD point, SPCmplx *outData, double rdiff, double sampleSpacing, long long nxInData, double * ikernel);
void ham1dx(double * data, int nx) ;
void sinc_kernel(int oversample_factor, int num_of_points, double resolution, double sampleSpacing, double *ikernel);
void generate_gaussian_ray_distribution(double xStdev,double yStdev, SPVector txPos, SPVector GRP, int nx,int ny, Ray *rayArray);


#endif /* threadCore_hpp */
