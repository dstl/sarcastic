/***************************************************************************
 *
 *       Module:    GOSS.h
 *      Program:    GOSS
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

#ifndef GOSS_GOSS_h
#define GOSS_GOSS_h

#include <SIlib.h>

#include <string.h>
#include <sys/time.h>
#if defined (__APPLE__) || defined(MACOSX)
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

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
#define MAXBOUNCES  8   // Maximum number of ray bounces

#define BOUNCEINDEXERROR    (-1)

typedef struct AABB {
    SPVector AA;
    SPVector BB;
} AABB;

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

typedef struct Texture {
    float ka;    // Constant for ambient reflection
    float kd;    // Constant for diffuse scattering
    float ks;    // Constant for specular scattering
    float n;     // Shininess constant
} Texture;


typedef struct TriCoords {
    SPVector A ;      // Cartesian coordinates of triangle
    SPVector B ;
    SPVector Cc ;
} TriCoords ;

typedef struct Hit {
    double dist;
    int trinum;
    double u;
    double v;
} Hit;

typedef struct Ray {
    SPVector org;    // Origin
    SPVector dir;    // Direction
    float  san;      // Solid Angle for Ray
    
} Ray;

typedef struct rangeAndPower {
    double range ;
    double power ;
} rangeAndPower ;

typedef struct cplxf {
    float r;
    float i;
} cplxf;

typedef struct Triangle {
    int  triNum;    // Triangle ID
    double d;         // Constant of plane equation
    double nd_u;      // Normal.u / normal.k
    double nd_v;      // normal.v / normal.k
    int k;          // projection dimension
    double kbu;
    double kbv;
    double kbd;
    double kcu;
    double kcv;
    double kcd;
    int textureInd;
} Triangle;


void banner () ;

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
    OCLPlatform platform ;
    int nTriangles ;
    Triangle * Triangles ;
    int nTextures ;
    Texture * Textures ;
    AABB SceneBoundingBox ;
    int nLeaves ;
    int * triPtrs ;
    int triListDataSize ;
    int * triListData ;
    int nTreeNodes ;
    KdData * KdTree ;
    int startPulse ;
    int nPulses ;
    int nAzBeam ;
    int nElBeam ;
    double dAz ;
    double dEl ;
    double Aeff ;               // The effective area of the Receive Antenna
    SPVector * TxPositions ;    // Pointer to beginning of TxPos data
    SPVector * RxPositions ;    // Pointer to beginning of RxPos data
    double * Fx0s ;             // Pointer to beginning of Fx0s data
    double * FxSteps ;          // Pointer to beginning of FxSteps data
    double raySolidAng ;
    double TxPowPerRay ;
    SPImage * phd ;             // Pointer to beginning of cphd data
    double chirpRate ;
    double ADRate ;
    double oneOverLambda ;
    double StartFrequency ;
    int bounceToShow ;
    SPStatus status ;

} threadData ;


int getUserInput(char **inCPHDFile, char **KdTreeFile, char **outCPHDFile,
                 int *startPulse, int *nPulses,
                 int *bounceToShow, int *nAzBeam, int *nElBeam, int *useGPU,
                 SPStatus *status) ;
char * loadProgramSource(const char *filename);
/*void radarRayBouncer(int iPulse,              // Pulse to process
                     RadarParams radar,         // struct for radar params
                     SPVector * RxPos,          // array - nPulses long
                     SPVector * TxPos,          // array - nPulses long
                     Triangle * Triangles,      // array - triangle data
                     Texture * textureData,     // array - texture data
                     AABB SceneBoundingBox,     // Scene bounding box
                     KdData * KdTree,           // array containing KdTree
                     int **triangleLists,       // array of triangle list arrays
                     cplxf * pulses             // output array - nPulses x nSamps
                     );*/
void * devPulseBlock ( void * threadArg ) ;


#endif
