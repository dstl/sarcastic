/***************************************************************************
 * 
 *           Module :  threadCore.hpp
 *          Program :  bircs
 *       Created by :  Darren Muff on 12/06/2017
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *   Description:
 *     Core thread routine for BIRCS
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
#ifndef threadCore_hpp
#define threadCore_hpp

#include <stdio.h>
#include <sarclib/sarclib.h>
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

typedef struct bircsThreadData {
    int tid ;
    int nThreads ;
    AABB SceneBoundingBox ;
    int startPulse ;
    int nPulses ;
    int nAzBeam ;
    int nElBeam ;
    SPVector * TxPositions ;    // Pointer to beginning of TxPos data
    SPVector * RxPositions ;    // Pointer to beginning of RxPos data
    double gainRx ;
    double freq_centre ;        // Centre frequency
    int bounceToShow ;
    SPStatus status ;
    double PowPerRay ;          // Ray Power (Pp = (Pt * Gtx ))
    int interrogate ;           // Do we want to write out details about an interrogation point in the scene?
    SPVector interogPt ;        // Position in scene coordinates of a point to be interogated
    double interogRad ;         // Radius in metres around interrogation point to calculate scattering for
    FILE ** interogFP ;         // File pointer to dump out interrogation data
    int polarisation ;
    TriangleMesh *sceneMesh ;    // Mesh of triangles in scene
    kdTree::KdData ** tree ;
    int treesize;
    ATS **accelTriangles ;
    SPVector * results ;
    int rayGenMethod ;
    
} bircsThreadData ;

enum  POLARISATION { VV, VH, HV, HH, V_, H_ } ;

void * devPulseBlock ( void * threadArg ) ;

void oclReflect(cl_context          context,            // OpenCL context - alrready built
                cl_command_queue    Q,                  // OpenCl command Q - already instatiated
                cl_kernel           kernel,             // OpenCl kernel for this routine to call
                cl_mem              dTriangles,
                int                 nRays,              // Number of rays to reflect
                size_t              localWorkSize,      // Local workgroupsize to use for OpenCL Kernel
                Ray                 *rays,              // Array of rays to consider
                Hit                 *hits,              // Array of hit points for each ray
                Ray                 *reflectedRays      // output array of reflected rays
);

void oclBuildShadowRays(cl_context          context,            // OpenCL context - alrready built
                        cl_command_queue    Q,                  // OpenCl command Q - already instatiated
                        cl_kernel           kernel,             // OpenCl kernel for this routine to call
                        size_t              localWorkSize,      // Local workgroupsize to use for OpenCL Kernel
                        int                 nRays,              // The number of reflected rays being considered
                        SPVector            RxPos,              // The Receiver location in x,y,z
                        Ray                 *reflectedRays,     // Array of reflected rays - used for their origin as its the reflection point to Rx
                        Ray                 *shadowRays,        // Output - shadow rays to be tested for occlusion later
                        double              *ranges             // Output - the range to the receiver for each shadowray. Calculated here to save calcs later
);

void oclRandomRays(cl_context context,          // OpenCL context - alrready built
                   cl_command_queue Q,          // OpenCl command Q - already instatiated
                   cl_kernel RRkernel,          // OpenCl kernel for this routine to call
                   int nAzBeam,                 // Number of azimuth rays in beam
                   int nElBeam,                 // Number of elevation rays in beam
                   size_t localWorkSize[2],     // local work size for this device - prealculated
                   double azStdDev,             // Az standard deviation of rays
                   double elStdDev,             // El standard deviation of rays
                   SPVector origin,             // Location of origin or rays
                   SPVector aimpoint,           // Mean aimpoint for rays
                   double raypow,               // Ray power (Pp = (Pt * Gtx / (4*PI)) * dAz * dEl)
                   Ray *rayArray                // Array that will be returned
);

void packSinc(SPCmplxD point, SPCmplx *outData, double rdiff, double sampleSpacing, long long nxInData, double * ikernel);
void ham1dx(double * data, int nx) ;
void sinc_kernel(int oversample_factor, int num_of_points, double resolution, double sampleSpacing, double *ikernel);
void generate_gaussian_ray_distribution(double xStdev,double yStdev, SPVector txPos, SPVector GRP, int nx,int ny, Ray *rayArray);


#endif /* threadCore_hpp */
