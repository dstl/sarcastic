/***************************************************************************
 *
 *       Module:    image_randim.c
 *      Program:    sarclib
 *   Created by:    Matt Notingham on 20/08/2007.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains functions associated with random numbers
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  19/07/2013
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

#include "sarclib.h"
#include <time.h>

#define LOOP(dtype, vtype) {                                                                \
                                                                                            \
    for(i = 0; i < s; i++) {                                                                \
        inp->data. dtype [i] = (vtype) (rand_gauss() * *(vtype *)sd + *(vtype *)mean);      \
    }                                                                                       \
}

#define CLOOP(dtype, vtype) {                                                               \
                                                                                            \
    for(i = 0; i < s; i++) {                                                                \
        inp->data. dtype [i].r = (rand_gauss() * ((vtype *)sd)->r) + ((vtype *)mean)->r;    \
        inp->data. dtype [i].i = (rand_gauss() * ((vtype *)sd)->i) + ((vtype *)mean)->i;    \
    }                                                                                       \
}

// Fill an array with gaussian noise - the mean & SD are assumed to be the same type as the image,
// hence all the scary casting in the macros above.
//
SPStatus *
im_gauss(SPImage * inp, void * mean, void * sd, SPStatus * status)
{
    int64_t s;
    int64_t i;
    
    CHECK_STATUS(status);
    CHECK_IMINTEG(inp, status);
    CHECK_PTR(mean, status);
    CHECK_PTR(sd, status);
    
    s = inp->nx * inp->ny;
    
    switch (inp->image_type)
    {
        case ITYPE_DOUBLE:
            LOOP(d, double);
            break;
        case ITYPE_FLOAT:
            LOOP(f, float);
            break;
        case ITYPE_INT64:
            LOOP(i64, int64_t);
            break;
        case ITYPE_INT32:
            LOOP(i32, int32_t);
            break;
        case ITYPE_INT16:
            LOOP(i16, int16_t);
            break;
        case ITYPE_INT8:
            LOOP(i8, int8_t);
            break;
            
        case ITYPE_CMPL_FLOAT:
            CLOOP(cmpl_f, SPCmplx);
            break;
        case ITYPE_CMPL_INT64:
            CLOOP(cmpl_i64, SPCmplxInt64);
            break;
        case ITYPE_CMPL_INT32:
            CLOOP(cmpl_i32, SPCmplxInt32);
            break;
        case ITYPE_CMPL_INT16:
            CLOOP(cmpl_i16, SPCmplxInt16);
            break;
        case ITYPE_CMPL_INT8:
            CLOOP(cmpl_i8, SPCmplxInt8);
            break;
            
        default:
            printf("Sorry, type %d is not supported\n", inp->image_type);
            status->status = INVALID_TYPE;
            return (status);
    }
    
    return(status);
}

// Initialise the seed of the random number generator - if the user supplies 0 as the
// seed then set the seed to be the current time.
//
SPStatus *
im_init_rand(long seed, SPStatus * status)
{
    CHECK_STATUS(status);
    
    if (seed == 0) {
        seed = time(NULL);
    }
    
    if (status->debug > 0) {
        printf("The seed has been set to %ld\n", seed);
    }
    
    srand48(seed);
    
    return(status);
}

// This uses the polar form of the Box-Muller transformation to generate
// gaussian random numbers with a mean of zero and a SD of 1.0
//
double
rand_gauss()
{
    double x1, x2, w, y1;
    do {
        x1 = 2.0 * drand48() - 1.0;
        x2 = 2.0 * drand48() - 1.0;
        w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);
    
    w = sqrt((-2.0 * log(w)) / w);
    
    y1 = x1 * w;
    
    return y1;
}

