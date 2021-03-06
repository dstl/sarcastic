/***************************************************************************
 * 
 *           Module :  interpolate.h
 *          Program :  sarclib
 *       Created by :  Alan Blake on 23/08/2005
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *       Sinc interpolation routines.
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

#ifndef sarclib_INTERPOLATE_H__
#define sarclib_INTERPOLATE_H__

int interpolate8pt(SPImage * data, SPImage * new_image);
int image_interpolate(SPImage * data, SPImage *map_data, double * k_data, 
                      int kernel_nx, int kernel_oversample, 
                      int use_kernel_x, int use_kernel_y, SPImage *o_data);
SPCmplx * interpolate_rotate(SPCmplx * data, int nx, int ny, double theta);

double * generate_sinc_kernel(int oversample_factor, int num_of_points);

#endif

