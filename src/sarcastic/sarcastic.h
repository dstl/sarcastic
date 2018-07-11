/***************************************************************************
 *
 *       Module:    sarcastic.h
 *      Program:    SARCASTIC
 *   Created by:    Darren on 27/07/2013.
 *                  Copyright (c) 2013 Dstl. All rights reserved.
 *
 *   Description:
 *   <ENTER DESCRIPTION HERE>
 *
 *
 *   CLASSIFICATION        :  <PENDING>
 *   Date of CLASSN        :  14/03/2013
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

#ifndef SARCASTIC_SARCASTIC_h
#define SARCASTIC_SARCASTIC_h

#include <sarclib/sarclib.h>
#include "SarcasticVersion.h"
#include <string.h>
#include <sys/time.h>
#if defined (__APPLE__) || defined(MACOSX)
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif
#include "OpenCLUtils.h"
#include "materialProperties.h"
#include "RCS.h"

#define EPSILON ((double)(0.00000001))  // Small number
#define NILROPE ((int) -666666 )
#define NOINTERSECTION -1
#define KDT_ISLEAF(n)    (n->leaf.flagDimAndOffset & (unsigned int)(1<<31))
#define KDT_DIMENSION(n) (n->leaf.flagDimAndOffset & 0x3)
#define KDT_OFFSET(n)    ((n->leaf.flagDimAndOffset & (0x7FFFFFFC))>>2)
#define BEAMMARGIN 2.0
#define GPUCAPABILITY_MAJOR 2
#define GPUCAPABILITY_MINOR 0
#define MAXDEVICES 8
#define RXLNA 30        // Receiver Low-Noise amplifier in dB
#define RXImp 50        // Receiver Impedence
#define MAXBOUNCES  10  // Maximum number of ray bounces
#define BOUNCEINDEXERROR    (-1)

typedef struct AABB {
    SPVector AA;
    SPVector BB;
} AABB;

typedef struct SPCmplxD{
    double r ;
    double i ;
} SPCmplxD;

typedef union {
    struct KdTreeLeaf {
        unsigned int flagDimAndOffset;
        // bits 0..1        : splitting dimension
        // bits 2..30       : offset bits
        // bits 31 (sign)   : flag whether node is a leaf
        float splitPosition;
        int Ropes[6];
        AABB aabb;
    } leaf;
    
    struct KdTreeBranch {
        unsigned int flagDimAndOffset;
        // bits 0..30       : offset to first child
        // bit 31 (sign)    : flag whether node is a leaf
        float splitPosition;
        int Ropes[6];
        AABB aabb;
    } branch;
    
} KdData ;

//typedef struct Texture {
//    float ka;    // Constant for ambient reflection
//    float kd;    // Constant for diffuse scattering
//    float ks;    // Constant for specular scattering
//    float n;     // Shininess constant
//} Texture;


//typedef struct TriCoords {
//    SPVector A ;      // Cartesian coordinates of triangle
//    SPVector B ;
//    SPVector Cc ;
//} TriCoords ;

typedef struct Hit {
    double dist;
    int trinum;
    double u;
    double v;
} Hit;

typedef struct Ray {
    SPVector org;    // Origin
    SPVector dir;    // Direction
    double   pow;    // Power for this ray
    double   len;    // Distance travelled to this ray's origin from transmission
    SPVector pol ;   // unit vector of direction of E field of ray
} Ray;

typedef struct rangeAndPower {
    double  range ;
    SPCmplxD Es ;
} rangeAndPower ;

typedef struct cplxf {
    float r;
    float i;
} cplxf;

typedef struct ATS {    // Accelerated Triangle Structure
    int  triNum;        // Triangle ID
    double d;           // Constant of plane equation
    double nd_u;        // Normal.u / normal.k
    double nd_v;        // normal.v / normal.k
    int k;              // projection dimension
    double kbu;
    double kbv;
    double kbd;
    double kcu;
    double kcv;
    double kcd;
    int    matInd;      // Material Index
} ATS;

typedef struct Triangle{
    int        id ;
    SPVector   AA ;
    SPVector   BB ;
    SPVector   CC ;
    SPVector   NN ;
    double     area ;
    double     globalToLocalMat[9];
    double     localToGlobalMat[9];
    int        matId;
} Triangle;

typedef struct KdTreeStruct {
    int                nTriangles;         // number of triangles in array 'triangles'
    ATS *              accelTriangles;     // Array of triangles of size nTriangles
    int                nTreeNodes;         // number of nodes in KdTree
    KdData *           KdTree;             // SAH - KdTree to optimise ray traversal through volume
    int                triListDataSize;    // size of trianglelist data
    int *              triangleListData;   // array of triangle indices into Triangles
    int                nLeaves;            // number of leaves (nodes with triangles) in KdTree
    int *              triangleListPtrs;   // array of pointers into triangleListData for each KdTree node
} KdTreeStruct ;



/*typedef struct RadarParams {
    int nPulses;                    // radar parameter - pulses in output pulses
    int nSamps;                     // radar parameter - samples in output pulses
    float oneOverLambda;            // radar parameter - one over wavelength
    double TxPowPerRay;             // radar parameter - transmitter power per ray in beam
    double TxRaySolidAngle;         // radar parameter - Solid angle of a single ray cast by transmitter
    float chirpRate;                // radar parameter - chirp rate Hz/Sec
    float ADRate;                   // radar parameter - ad sample rate
    float dAz;                      // radar parameter - az slice in radians
    float dEl;                      // radar parameter - el slice in radians
    int nAzBeam;                    // num of rays in azimuth in radar beam
    int nElBeam;                    // num of rays in elevation in radar beam
    int bounceToShow ;              // Which bounce to display -1= none, 0=1st bounce etc...
    double Aeff;                    // Radar Parameter - The effective area of the Receive Antenna
} RadarParams;*/

typedef struct threadData {
    int devIndex ;
    int nThreads ;
    OCLPlatform platform ;
    KdTreeStruct       KDT ;                // Structure containing all KDTree info
    AABB SceneBoundingBox ;
    int startPulse ;
    int nPulses ;
    int nAzBeam ;
    int nElBeam ;
    SPVector * TxPositions ;    // Pointer to beginning of TxPos data
    SPVector * RxPositions ;    // Pointer to beginning of RxPos data
    double * Fx0s ;             // Pointer to beginning of Fx0s data
    double * FxSteps ;          // Pointer to beginning of FxSteps data
    double * amp_sf0 ;          // Pointer to beginning of amp_sf0 data
    double gainRx ;
    SPImage * phd ;             // Pointer to beginning of cphd data
    double chirpRate ;
    double ADRate ;
    double pulseDuration ;      // Pulse duration in seconds
    double oneOverLambda ;
    double freq_centre ;        // Centre frequency
    double StartFrequency ;
    int bounceToShow ;
    SPStatus status ;
    double PowPerRay ;          // Ray Power (Pp = (Pt * Gtx ))
    int interrogate ;           // Do we want to write out details about an interrogation point in the scene?
    SPVector interogPt ;        // Position in scene coordinates of a point to be interogated
    double interogRad ;         // Radius in metres around interrogation point to calculate scattering for
    FILE ** interogFP ;         // File pointer to dump out interrogation data
    Triangle *triangles;        // Array of triangle coordinates
} threadData ;

typedef struct rnpData_t {
    SPCmplxD Es;
    double  samplingOffset;
    int     samplingOffsetInt;
    double  indexOffset ;
    double  rdiff;
}rnpData_t ;

typedef struct threadDataBF {
    int nx ;
    int nrnpItems ;
    int startSamp ;
    int nSamp ;
    rnpData_t *rnpData ;
    int pulseIndex ;
    SPImage * phd ;
    double A ;
    double B ;
} threadDataBF ;


int getUserInput(char **inCPHDFile, char **KdTreeFile, char **outCPHDFile,
                 int *startPulse, int *nPulses,
                 int *bounceToShow, int *nAzBeam, int *nElBeam, int *useGPU,
                 int *interrogate, SPVector *interogPt, double *interograd,
                 FILE **interogateFP, int * pulseUndersampleFactor, SPStatus *status) ;

void * devPulseBlock ( void * threadArg ) ;
void * beamForm ( void * threadArg ) ;
void oclRayTrace(cl_context         context,            // OpenCL context - already instantiated
                 cl_command_queue   Q,                  // OpenCl command Q - already instatiated
                 cl_kernel          RTkernel,           // Ray Tracing kernel, already compiled and created from a program
                 size_t             globalWorkSize[2],  // Total number of rays in each dimension to cast
                 size_t             localWorkSize[2],   // Work dimensions for this device
                 
                 KdTreeStruct       KDT,                // Structure containing all KDTree info
                 
                 SPVector           RxPos,              // Receive position of radar (used to calculate range)
                 double             gainRx,             // receive antenna gain
                 double             TxPowPerRay,        // TxPowerPerRay
                 AABB               SceneBoundingBox,   // Bounding box of scene - required for ray traversal optimisation
                 
                 int                bounceToShow,       // Useful for debugging
                 int                pulseIndex,          // Useful for debugging
                 
                 Ray *              rays,               // Array of rays (size nAzBeam*nElBeam). Each ray will be cast through KdTree
                 rangeAndPower      *rnp                // output array of ranges and powers for each ray intersection
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

void oclKdTreeHits(cl_context         context,            // OpenCL context - already instantiated
                   cl_command_queue   Q,                  // OpenCl command Q - already instatiated
                   cl_kernel          STkernel,           // stacklessTraverse kernel, already compiled and created from a program
                   int                nRays,              // Total number of rays to cast
                   size_t             localWorkSize,      // Work dimensions for this device
                   cl_mem             dTriangles,
                   cl_mem             dKdTree,
                   cl_mem             dtriListData,
                   cl_mem             dtriListPtrs,
                   AABB               SceneBoundingBox,   // Bounding box of scene - required for ray traversal optimisation
                   Ray *              rays,               // Array of rays (size nRays). Each ray will be cast through KdTree
                   Hit *              hits                // output array of hit locations
);

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

void oclReflectPower(cl_context          context,            // OpenCL context - alrready built
                     cl_command_queue    Q,                  // OpenCl command Q - already instatiated
                     cl_kernel           kernel,             // OpenCl kernel for this routine to call
                     size_t              localWorkSize,      // Local workgroupsize to use for OpenCL Kernel
                     cl_mem              dTriangles,
                     cl_mem              dTextures,
                     Hit                 *hits,              // Array of hit locations to x-ref with triangles (and then Textures) for material props
                     SPVector            RxPos,              // Location of Receiver in x,y,z
                     double              GrxOverFourPi,      // Receiver antenna gain / 4Pi.
                     int                 nRays,              // The number of reflected rays being considered
                     Ray                 *rays,              // unit vector rays arriving at hitpoint
                     Ray                 *reflectedRays,     // Array of reflected rays - used for their origin as its the reflection point to Rx
                     Ray                 *shadowRays,        // Viewpoint unit vector. The vector direction to calculate power for (usually to receiver)
                     double              *ranges,            // Range to receiver for each shadow ray (precalculated in shadowRay generation)
                     rangeAndPower       *rnp                // Output array of ray power at, and range to reciever
);

void oclPOField(cl_context          context,            // OpenCL context - already built
                cl_command_queue    Q,                  // OpenCl command Q - already instatiated
                cl_kernel           kernel,             // OpenCl kernel for this routine to call
                size_t              localWorkSize,      // Local workgroupsize to use for OpenCL Kernel
                Triangle            *triangles,         // Array of triangles
                int                 ntris,              // Number of triangles
                Hit                 *hits,              // Array of hit locations to x-ref with triangles for material props
                int                 nRays,              // The number of reflected rays being considered
                Ray                 *rays,              // unit vector rays arriving at hitpoint
                Ray                 *shadowRays,        // Array of reflected rays - used for their origin as its the reflection point to Rx
                SPVector            RxPos,              // Location of Receiver in x,y,z
                double              k,                  // Wavenumber constant k = 2 * PI / Lambda
                double              *ranges,            // Range to receiver for each shadow ray (precalculated in shadowRay generation)
                double              gainRX,             // Receiver gain
                int                 firstBounce,        // If 1 then origin is used for phase calculation in Filed equations causing triangles to correctly combine
                rangeAndPower       *rnp                // Output array of ray power at, and range to reciever
);

int buildKernel(cl_context context,             // OpenCL Context. Already created
                const char *kernelCodePath,     // String defining FQPN for kernel code
                const char *kernelCodeName,     // OpenCL kernel name for kernel code
                cl_device_id devID,             // Device Id on this platform to compile for
                int        workDims,            // WorkDimensions to optimise localWorkSize for
                cl_program *program,            // Output - Compiled OCL programme
                cl_kernel  *kernel,             // Output - OpenCL Kernel
                size_t     *localWorkSize       // Output - Prefered local worksize for kernel
);

void banner () ;

void cpuPOField(Triangle            *triangles,         // Array of triangles
                int                 ntris,              // Number of triangles
                Hit                 *hits,              // Array of hit locations to x-ref with triangles for material props
                int                 nRays,              // The number of reflected rays being considered
                Ray                 *rays,              // unit vector rays arriving at hitpoint
                Ray                 *shadowRays,        // Array of reflected rays - used for their origin as its the reflection point to Rx
                SPVector            RxPos,              // Location of Receiver in x,y,z
                double              k,                  // Wavenumber constant k = 2 * PI / Lambda
                double              *ranges,            // Range to receiver for each shadow ray (precalculated in shadowRay generation)
                double              gainRx,             // Receiver gain used for power calculations
                int                 firstBounce,        // if 1 then PO calcs use origin for field calculations
                rangeAndPower       *rnp                // Output array of ray power at, and range to reciever
);

#endif
