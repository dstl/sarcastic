/***************************************************************************
 *  boxMuller_Random.c
 *  triDecomp
 *
 *  Created by Muff Darren on 27/09/2014.
 *  Copyright (c) 2014 [dstl]. All rights reserved.
 *
 * Original source code
 *  (c) Copyright 1994, Everett F. Carter Jr.
 *  Permission is granted by the author to use
 *  this software for any application provided this
 *  copyright notice is preserved.
 *
 *   Description:
 *      Implements the Polar form of the Box-Muller Transformation
 *      This module is taken from http://www.taygeta.com/random/gaussian.html
 *      There is a useful discussion on generating Gaussian random numbers
 *      at http://www.design.caltech.edu/erik/Misc/Gaussian.html
 *
 *
 *  CLASSIFICATION       :   UNCLASSIFIED
 *  Date of CLASSN       :   25/09/2014
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
 ***************************************************************************/

#ifndef sarcastic_boxMullerRandom_h
#define sarcastic_boxMullerRandom_h

#include <math.h>
#include "ranf.h"

float box_muller(float mean, float standardDeviation);	/* normal random variate generator */

#endif