/***************************************************************************
 * 
 *           Module :  boxMullerRandom.h
 *          Program :  Sarcastic
 *       Created by :  Darren Muff on 27/09/2014
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *   Description:
 *      Implements the Polar form of the Box-Muller Transformation
 *      This module is taken from http://www.taygeta.com/random/gaussian.html
 *      There is a useful discussion on generating Gaussian random numbers
 *      at http://www.design.caltech.edu/erik/Misc/Gaussian.html
 *     Original source code (c) Copyright 1994, Everett F. Carter Jr.
 *       Permission is granted by the author to use
 *       this software for any application provided this
 *       copyright notice is preserved.
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

#ifndef sarcastic_boxMullerRandom_h
#define sarcastic_boxMullerRandom_h

#include <math.h>
#include "ranf.h"

float box_muller(float mean, float standardDeviation);	/* normal random variate generator */

#endif
