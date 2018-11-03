/***************************************************************************
 * 
 *           Module :  ranf.h
 *          Program :  Sarcastic
 *       Created by :  Darren Muff on 27/09/2014
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *     Description : 
 *|The routines below implement an interface based on the drand48 family of
 *|random number routines, but with seed handling and default starting seed
 *|compatible with the CRI MATHLIB random number generator (RNG) family "ranf".
 *|
 *|Original source code  :
 *|    http://w3.pppl.gov/~hammett/comp/python/LLNLDistribution11/RNG/Src/
 *|
 *|-----------------------------------------------------------------------
 *|
 *|User-callable routines defined:
 *|   Seedranf - Set seed from 32-bit integer
 *|   Mixranf  - Set seed, with options; return seed
 *|   Getranf  - Get 48-bit seed in integer array
 *|   Setranf  - Set seed from 48-bit integer
 *|   Getmult  - Get 48-bit multiplier in integer array
 *|   Setmult  - Set multiplier from 48-bit integer
 *|   Ranf     - The generator itself
 *|
 *|Currently responsible person:
 *|  Fred N. Fritsch
 *|   Computer Applications Organization
 *|   LCPD, ICF Group
 *|   Lawrence Livermore National Laboratory
 *|   fnf@llnl.gov
 *|   Darren Muff, Defence Science and Technology Laboratory, UK
 *| Revision history:
 *|  (yymmdd)
 *|   91????  DATE WRITTEN
 *|           This started with the ranf.c that was checked out of the
 *|           Basis repository in August 1996.  It was written by Lee
 *|           Busby to provide an interface for Basis to the drand48
 *|           family of generators that is compatible with the CRI
 *|           MATHLIB ranf.  Following are the relevant cvs log entries.
 *|   950511  Added back seedranf, setranf, getranf, mixranf from older
 *|           Basis version. (LB)
 *|   950530  Changed type of u32 from long to int, which is currently ok
 *|           on all our supported platforms. (LB)
 *|   960823  Revised to use the portable PMATH generator instead. (FNF)
 *|   960903  Added new routines Getmult and Setmult. (FNF)
 *|   960904  Changed declarations for 48-bit quantities from pointers to
 *|           two-element arrays, to be consistent with actual usage, and
 *|           brought internal documentation up to LDOC standard. (FNF)
 *|   960905  Moved definitions of internal procedures PM_16to24 and
 *|           PM_24to16 to pmath_rng.c (see ranf.h).  (FNF)
 *|   960911  Corrected some problems with return values from Mixranf.
 *|           Also improved calling sequence descriptions for Seedranf and
 *|           Mixranf. (FNF)
 *|   960916  Eliminated state variable ranf_started.  Since the PMATH
 *|           family has the Cray default starting values built in, this
 *|           is no longer needed. (FNF)
 *|   961011  Modifed Setranf and Setmult to work OK on Cray's when given
 *|           values saved on a 32-bit workstation. (FNF)
 *|   961212  Corrected dimension error in Getmult.  (FNF per E. Brooks)
 *| 27-09-14  Adapted to work with OpenCK routines and scan primitives (DGM)
 * -----------------------------------------------------------------------
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

#ifndef sarcastic_ranf_h
#define sarcastic_ranf_h

#include <stdlib.h>
#include <math.h>
#ifdef RAN_DEBUG
#include <stdio.h>
#endif

#ifdef __MWERKS__
/*#include <utime.h>*/
#include <time.h>
#include <Timer.h>
#else
#if defined(_WIN32)
#include <time.h>
#else
#include <sys/time.h>
#endif
#endif

typedef unsigned int u32;
typedef unsigned short int u16;
typedef double f64;

/*  Prototypes for routines defined in ranf.c  */
#ifdef __STDC__
void Seedranf(u32 *s);            /* Set seed from 32-bit integer */
void Mixranf(int *s, u32 s48[2]); /* Set seed, with options; return seed */
void Getranf(u32 s48[2]);         /* Get 48-bit seed in integer array */
void Setranf(u32 s48[2]);         /* Set seed from 48-bit integer */
void Getmult(u32 m48[2]);         /* Get 48-bit multiplier in integer array */
void Setmult(u32 m48[2]);         /* Set multiplier from 48-bit integer */
f64  Ranf();                      /* The generator itself */
#else
void Seedranf();
void Mixranf();
void Getranf();
void Setranf();
void Getmult();
void Setmult();
f64  Ranf();
#endif

/*  Prototypes for routines defined in pmath_rng.c  */
#ifdef __STDC__
void PM_16to24(u16 x16[3], double x24[2]);
/* Convert 3 16-bit shorts to 2 24-bit doubles */
void PM_24to16(double x24[2], u16 x16[3]);
/* Convert 2 24-bit doubles to 3 16-bit shorts */
void PM_GSeed(double seedout[2]);  /* Get the current seed */
void PM_SSeed(double seedin[2]);   /* Reset the seed (unsafe) */
void PM_GMult(double multout[2]);  /* Get the current multiplier */
void PM_SMult(double multin[2]);   /* Reset the multiplier (unsafe) */
f64  PM_RANF();                    /* The generator itself */
#else
void PM_16to24();
void PM_24to16();
void PM_GSeed();
void PM_SSeed();
void PM_GMult();
void PM_SMult();
f64  PM_RANF();
#endif

#endif
