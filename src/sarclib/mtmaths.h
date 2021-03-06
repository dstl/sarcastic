/***************************************************************************
 * 
 *           Module :  mtmaths.h
 *          Program :  sarclib
 *       Created by :  Darren Muff on Sat Jun 30 10:34:25 2018 +0100
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      This file contains information for the maths functions used in the time
 *      domain processor
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

#ifndef sarclib_MTMATHS_H__
#define sarclib_MTMATHS_H__

#define SP_MIN(a, b) ((a) > (b) ? (b) : (a))   ///< calculates the minimum of 2 numbers 
#define SP_MAX(a, b) ((a) > (b) ? (a) : (b))   ///< calculates the maximum of 2 numbers 
#define SP_SIGN(a) ((a) > 0.0 ? 1.0 : -1.0)    ///< calculates the sign of a number, useful for working out the geometry 
#define SP_SQR(a) (a*a)                        ///< squares what is in the brackets, useful if "a" is a long expression 

#endif
