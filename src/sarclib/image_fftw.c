/***************************************************************************
 *
 *       Module:    image_fftw.c
 *      Program:    sarclib
 *   Created by:    Matt Nottingham in 2012.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains functions that perform talks to the FFTW library
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
#include <fftw3.h>

extern pthread_mutex_t fftw_plan_lock;

SPStatus *
im_fftw_r2c(SPImage * inp, FFTMODE mode, SPImage * out, SPStatus * status)
{
  fftwf_plan p = NULL;
  int flags;
  SPImage save_inp;

  CHECK_STATUS(status);

  if (out->image_type != ITYPE_CMPL_FLOAT && inp->image_type != ITYPE_FLOAT) {
    printf("im_fftw_r2c: Only ITYPE_FLOAT input (found %s) & ITYPE_CMPL_FLOAT output (found %s) are supported.\n",
	   itype2string(inp->image_type), itype2string(out->image_type));
    status->status = INVALID_TYPE;
  }
  
  if (inp->ny != 1) {
    printf("Can only handle 1D real to complex FFTs at the moment\n");
    status->status = INPUT_NY_MISMATCHED;
  }

  if (out->data.v != NULL) {
    if (out->nx != (inp->nx / 2 + 1)) {
      printf("NX of output data is not correct size (it is %"PRId64", it should be %"PRId64")\n",
	     out->nx, (inp->nx / 2 + 1));
      status->status = OUTPUT_NX_MISMATCHED;
    }
    if (out->ny != inp->ny) {
      printf("NY output data is not correct size, (it is %"PRId64", it should be %"PRId64")\n",
	     out->ny, inp->ny);
      status->status = OUTPUT_NY_MISMATCHED;
    }
  } else {
    im_create(out, ITYPE_CMPL_FLOAT, (inp->nx / 2 + 1), 1, 1.0, 1.0, status);
  }

  CHECK_STATUS(status);
 
  if (mode & FFT_FFTW_MEASURE) {
    flags = FFTW_MEASURE;
  } else if (mode & FFT_FFTW_PATIENT) {
    flags = FFTW_PATIENT;
  } else {
    flags = FFTW_ESTIMATE;
  }

  if (flags != FFTW_ESTIMATE) {
    im_init(&save_inp, status);
    im_copy(inp, &save_inp, status);
  }

  pthread_mutex_lock(&fftw_plan_lock);

  p = fftwf_plan_dft_r2c_1d((int)inp->nx, inp->data.f, (fftwf_complex *) out->data.cmpl_f,
			    flags);

  pthread_mutex_unlock(&fftw_plan_lock);

  if (p != NULL) {
    if (flags != FFTW_ESTIMATE) {
      im_copy(&save_inp, inp, status);
      im_destroy(&save_inp, status);
    }
    fftwf_execute(p);
  } else {
    status->status = INVALID_FFT_SIZE;
    fprintf(stderr, "Failed to set a FFT plan!\n");
  }

  pthread_mutex_lock(&fftw_plan_lock);
  fftwf_destroy_plan(p);
  pthread_mutex_unlock(&fftw_plan_lock);
  
  return status;
}


SPStatus* im_fftw(SPImage *a, FFTMODE mode, SPStatus * status)                     // fft an image
{
    fftwf_plan p = NULL;
    fftwf_r2r_kind kind[2];
    int sz[2];
    int sgn;
    double num_scale = 1.0;
    int flags;
    SPImage save_inp;

    CHECK_STATUS(status);
    CHECK_IMINTEG(a,status);
    
    if (a->image_type != ITYPE_CMPL_FLOAT && a->image_type != ITYPE_FLOAT) {
        printf("im_fftw: Only ITYPE_FLOAT & ITYPE_CMPL_FLOAT are supported - data type is %s\n",
               itype2string(a->image_type));
        status->status = INVALID_TYPE;
    }
    CHECK_STATUS(status);

    if (mode & FFT_FFTW_MEASURE) {
      flags = FFTW_MEASURE;
    } else if (mode & FFT_FFTW_PATIENT) {
      flags = FFTW_PATIENT;
    } else {
      flags = FFTW_ESTIMATE;
    }

    if (flags != FFTW_ESTIMATE) {
      im_init(&save_inp, status);
      im_copy(a, &save_inp, status);
    }
    
    if (a->image_type == ITYPE_CMPL_FLOAT) {
        if (mode & REV) {
            sgn = FFTW_BACKWARD;
        } else {
            sgn = FFTW_FORWARD;
        }
        
        if (mode & FFT_X_ONLY && !(mode & FFT_Y_ONLY)) {
            sz[0] = (int)a->nx;
            sz[1] = (int)a->ny;
            pthread_mutex_lock(&fftw_plan_lock);
            
            p = fftwf_plan_many_dft(1, sz, (int)a->ny, (fftwf_complex *) a->data.cmpl_f, NULL, 1, (int)a->nx, 
				    (fftwf_complex *) a->data.cmpl_f, NULL, 1, (int)a->nx, sgn, flags);
            pthread_mutex_unlock(&fftw_plan_lock);
            num_scale = a->nx;
        }
        if (mode & FFT_Y_ONLY && !(mode & FFT_X_ONLY)) {
            im_transp(a, status);
            sz[0] = (int)a->nx;
            sz[1] = (int)a->ny;
            
            pthread_mutex_lock(&fftw_plan_lock);
            p = fftwf_plan_many_dft(1, sz, (int)a->ny, (fftwf_complex *)a->data.cmpl_f, NULL, 1, (int)a->nx, 
				    (fftwf_complex *) a->data.cmpl_f, NULL, 1, (int)a->nx, sgn, flags);
            pthread_mutex_unlock(&fftw_plan_lock);
            
            num_scale = a->nx;
            
        }
        if (mode & FFT_X_ONLY && mode & FFT_Y_ONLY) {
            pthread_mutex_lock(&fftw_plan_lock);
            p = fftwf_plan_dft_2d((int)a->ny, (int)a->nx, (fftwf_complex *) a->data.cmpl_f, (fftwf_complex *)a->data.cmpl_f, sgn, FFTW_ESTIMATE);
            pthread_mutex_unlock(&fftw_plan_lock);
            num_scale = a->nx * a->ny;
            
        }
        
    } else if (a->image_type == ITYPE_FLOAT) {
        if (mode & REV) {
            kind[0] = FFTW_REDFT01;
        } else {
            kind[0] = FFTW_REDFT10;
        }
        
        if (mode & FFT_X_ONLY && !(mode & FFT_Y_ONLY)) {
            sz[0] = (int)a->nx;
            sz[1] = (int)a->ny;
            
            pthread_mutex_lock(&fftw_plan_lock);
            p = fftwf_plan_many_r2r(1, sz, (int)a->ny, a->data.f, NULL, 1, (int)a->nx, a->data.f, NULL, 1, (int)a->nx, kind, flags);
            pthread_mutex_unlock(&fftw_plan_lock);
            
            num_scale = a->nx;
        }
        
        if (mode & FFT_Y_ONLY && !(mode & FFT_X_ONLY)) {
            im_transp(a, status);
            sz[0] = (int)a->nx;
            sz[1] = (int)a->ny;
            
            pthread_mutex_lock(&fftw_plan_lock);
            p = fftwf_plan_many_r2r(1, sz, (int)a->ny, a->data.f, NULL, 1, (int)a->nx, a->data.f, NULL, 1, (int)a->nx, kind, FFTW_ESTIMATE);
            pthread_mutex_unlock(&fftw_plan_lock);
            
            num_scale = a->nx;
        }
        
        if (mode & FFT_X_ONLY && mode & FFT_Y_ONLY) {
            pthread_mutex_lock(&fftw_plan_lock);
            
            p = fftwf_plan_r2r_2d((int)a->ny, (int)a->nx, a->data.f, a->data.f, kind[0], kind[0], flags);
            pthread_mutex_unlock(&fftw_plan_lock);
            num_scale = a->nx * a->ny;
        }
    }

    CHECK_STATUS(status);
    
    if (p != NULL) {
      if (flags != FFTW_ESTIMATE) {
	im_copy(&save_inp, a, status);
	im_destroy(&save_inp, status);
      }
      fftwf_execute(p);
    } else {
        status->status = INVALID_FFT_SIZE;
        fprintf(stderr, "Failed to set a FFT plan!\n");
        return status;
    }
    
    if (mode & FFT_Y_ONLY && !(mode & FFT_X_ONLY)) {
        im_transp(a, status);
    }
    
    if((mode & (FFT_X_ONLY+FFT_Y_ONLY)) && (mode & (SCALE_N+SCALE_R))) {		// need to scale
        double scale=1.0;
        
        // Calculate the scale factor from the number of FFT points and if its scale_R or _N
        //
        scale/= (double) ((mode & SCALE_R) ? sqrt((double)num_scale) : num_scale);
        
        im_mult_scalar(a, ITYPE_DOUBLE, &scale, status);
    }
    
    CHECK_STATUS(status);

    pthread_mutex_lock(&fftw_plan_lock);
    fftwf_destroy_plan(p);
    pthread_mutex_unlock(&fftw_plan_lock);
    
    return(status);
}



