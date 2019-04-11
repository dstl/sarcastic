/***************************************************************************
 * 
 *           Module :  tdp.h
 *          Program :  tdpocl
 *       Created by :  Darren Muff on 12/03/2013
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *    Main routine for Time Domain Processor.
 *    This module is responsible for collecting user input, determining
 *    memory requirements and then farming out the work to key tdp_core
 *    processing modules.
 *
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
 
#ifndef __TDP_H
#define __TDP_H

#include <sys/time.h>
#include <unistd.h>
#include <fftw3.h>
#include "tdpoclVersion.h"
#include <sarclib/sarclib.h>
#include "tdpcore.h"
#include "OpenCLUtils.h"

/* These define the number of points in the sinc interpolate and the
   over sample factor */

#define PULSESPERPULSEBLOCK (8192)
#define RCOMPOVERSAMPLE ((float)1.3)

/* define minimum GPU capability required to run the programme */
#define GPUCAPABILITY_MAJOR 2
#define GPUCAPABILITY_MINOR 0
#define MAXDEVICES 4

/* Specify a defualt EIRP. Required for RCS correction if used */
#define DEFAULT_EIRP ((double)1.567417e+08)

struct oclLWS { char *devName; size_t x; size_t y; } ;


/* This structure is used when we are running in K-space trimming mode. The
   structure holds the data for an individual surface/output pixel and the 
   bounds of the K-space (in pulses & samples) that should be used when
   forming the image */
/* lowpulse[lowsamp, highsamp] - highpulse[lowsamp, highsamp] */

typedef struct {
  int lowpulse;
  int highpulse;
  int lowpulse_samples[2];
  int highpulse_samples[2];
} KSpaceSupport;
#define TTM_LINE_SIZE (512)


/* A quick'n'dirty structure to hold all the KSpaceSupport data - should be moved into
   image.h at somepoint......*/
typedef struct {
  int nx;
  int ny;
  KSpaceSupport * data;
} KSpaceSupportImage;

/* Structure that holds all the information that is required by a thread
   for the thread to do its work. Each thread gets a unique instance of this
   structure */

typedef struct {
  struct timeval start;   /* The time we started to process */
  SPImage * image;        /* The output image */
  SPImage * surface;      /* The focal surface */
  SPImage * data;         /* The range compressed radar data */
  SPStatus * status;      /* Status */
  CPHDHeader *hdr;        /* CPHD header information */
  KSpaceSupportImage * kbounds;  /* The K-space bounds */
  double * kernel;        /* Sinc interpolate kernel */
  double rng_off;
  int starty;             /* The first row (of output image) that this thread should work on */
  int stopy;              /* The last row this thread should work on */
  int pulse;              /* What pulse are we working on (in CPHD file numbering scheme) */
  int first_pulse;        /* The pulse we start on (CPHD num scheme) */
  int last_pulse;         /* The last pulse (CPHD num scheme) */
  int id;                 /* Unique ID number of this thread */
  int nPulsesPerBlock;           /* The number of pulses in a block to process */
  int dataStart;          /* Start pulse od data chunk in data image*/
  int gpuLabel;           // If using GPU this give the selected card id
} ThreadInfo;

/* function prototypes */
void topotrim_setup(FILE **topotrim_fp,
                    char **input_filename,
                    char **surface_fname,
                    int *NX, int *NY, int *xstart, int *ystart);

int create_surface(SPVector sc_pos,
                   SPVector sensor_pos,
                   char * surface_fname,
                   int nx,
                   int ny,
                   double sx,
                   double sy,
                   SPStatus * status);

void getUserInput(int *USEGPU,              // Whether to use the GPU if its available or run on CPU
                  char **cphdFilename,      // Name of CPHD file to process
                  char **surfaceFilename,   // Name of surface file - defining each pixel in ECEF coords
                  int *Nx,                  // Size of output image in X
                  int *Ny,                  // Size of output image in Y
                  int *xstart,              // First x sample in surface file to use
                  int *ystart,              // First y sample in surface file to use
                  int *nPulses,             // Number of pulses in CPHD file to process
                  int *startPulse,          // First pulse to process
                  double *PtGt,             // If non zero then use this value of Pt*Gt to output image as RCS (m^2)
                  char **outputFilename,    // Name of file to store the output image
                  SPStatus *status );       // mtlib status flag
size_t getMemorySize( );
void oclSRPRanges(OCLPlatform platform,     // OpenCl Platform definition structure
                  int nPulses,              // Number of pulses to convert
                  SPVector *TXs,            // Array of transmit locations. One for each pulse
                  SPVector *RXs,            // Array of receive locations. One for each pulse
                  SPVector *SRP,            // Array of SRP locations. One for each pulse
                  double **ranges,          // Pointer to array of ranges to be filled. One for each pulse
                  SPStatus *status          // Status array
);
void banner ();


//#ifdef __CUFILE__
//extern "C" {
//#endif
//void GPUCheck(int **gpuLabel,int *gpuCards);
//void * gpuWork(void * info);
//void * gpuKspaceWork(void * info);
//#ifdef __CUFILE__
//}
//#endif

#endif
