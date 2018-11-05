/***************************************************************************
 * 
 *           Module :  image_interp.h
 *          Program :  sarclib
 *       Created by :  Mark Ashe on 14/08/2007
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      This file contains functions needed to interpolate images
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

#ifndef sarclib_IMAGE_INTERP_H__
#define sarclib_IMAGE_INTERP_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "sarclib.h"

// Function prototypes
//
int generate_axis(double* axis, double axis_start, double axis_stop, long int np, double* output_resolution);
void simple_tpose(double* in, long int nx, long int ny, double* out);
long int index_below(double value, double* axis, long int np);
long int index_above(double value, double* axis, long int np);
int index_span(double minv, double maxv, double *axis, long int np, long int *below, long int* above);
int index_bracketing(double value, double* axis_in, long int axis_length, long int* low_bound, long int* high_bound);
int spline(double* in_x, double* in_y, long int np, double *out_y2);
double splint(double out_x, double* in_x, double* in_y, long int np, double* in_y2);
int reg_spline(double *in_x, double *in_y, long int np, double *out_y2);
int reg_splint(double *out_x, double *out_y,long int np, double x0, double x1, double y0, double y1, double dy0, double dy1);
int bilinterp_image(SPImage* inp_img, SPImage* ret_img, double* in_ties, double* out_ties, long int out_nx, long int out_ny, SPStatus* status);
int splinterp_image(SPImage* inp_img, SPImage* ret_img, double* in_ties, double* out_ties, long int out_nx, long int out_ny, SPStatus* status);
double blint_dbl(double x_out, double y_out, double x0, double x1, double y0, double y1, double z00, double z01, double z10, double z11);
int bilinterp_vecim(SPImage *inp_img, SPImage *ret_img, double *out_ties, long int out_nx, long int out_ny, SPStatus *status);
int splinterp_vecim(SPImage *inp_img, SPImage *ret_img, double *out_ties, long int out_nx, long int out_ny, SPStatus *status);

#endif
