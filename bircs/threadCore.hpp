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
    int devIndex ;
    int nThreads ;
    OCLPlatform platform ;
    AABB SceneBoundingBox ;
    int startPulse ;
    int nPulses ;
    int nAzBeam ;
    int nElBeam ;
    SPVector * TxPositions ;    // Pointer to beginning of TxPos data
    SPVector * RxPositions ;    // Pointer to beginning of RxPos data
    //    double * Fx0s ;             // Pointer to beginning of Fx0s data
    //    double * FxSteps ;          // Pointer to beginning of FxSteps data
    //    double * amp_sf0 ;          // Pointer to beginning of amp_sf0 data
    double gainRx ;
    //    double chirpRate ;
    //    double ADRate ;
    //    double pulseDuration ;      // Pulse duration in seconds
    //    double oneOverLambda ;
    double freq_centre ;        // Centre frequency
    //    double StartFrequency ;
    int bounceToShow ;
    SPStatus status ;
    double PowPerRay ;          // Ray Power (Pp = (Pt * Gtx ))
    int interrogate ;           // Do we want to write out details about an interrogation point in the scene?
    SPVector interogPt ;        // Position in scene coordinates of a point to be interogated
    double interogRad ;         // Radius in metres around interrogation point to calculate scattering for
    FILE ** interogFP ;         // File pointer to dump out interrogation data
    TriangleMesh *sceneMesh ;    // Mesh of triangles in scene
} bircsThreadData ;

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
void buildRays(Ray **rayArray, int *nRays, int nAzRays, int nElRays, TriangleMesh *mesh, SPVector TxPos,
               double PowPerRay, AABB SceneBoundingBox,
               cl_context context, cl_command_queue commandQ, cl_kernel  randRaysKL, size_t randRaysLWS[2],
               SPVector **rayAimPoints);

#endif /* threadCore_hpp */
