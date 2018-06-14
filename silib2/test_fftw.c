/***************************************************************************
 *
 *       Module:    image_fftw.c
 *      Program:    SIlib2
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

#include "SIlib2.h"
#include <fftw3.h>

int
main(int argc, char *argv[])
{
  SPStatus status;
  SPImage test_im_meas;
  SPImage test_im_est;
  SPImage test_im_def;
  SPImage diff_est;
  SPImage diff_meas;
  int64_t x;

  status.debug = 0;
  im_init_lib(&status, "test_fftw", argc, argv);

  srandom(10182651);
  im_create(&test_im_def, ITYPE_CMPL_FLOAT, 16384 * 8, 1, 1.0, 1.0, &status);

  for(x = 0; x < test_im_def.nx; x++) {
    test_im_def.data.cmpl_f[x].r = (random() - RAND_MAX/2) * 100.0 / RAND_MAX;
    test_im_def.data.cmpl_f[x].i = (random() - RAND_MAX/2) * 100.0 / RAND_MAX;
  }

  im_init(&test_im_est, &status);
  im_init(&test_im_meas, &status);

  im_copy(&test_im_def, &test_im_est, &status);
  im_copy(&test_im_def, &test_im_meas, &status);

  im_fft(&test_im_def, FFT_X_ONLY + FWD + SCALE_N, &status);
  test_im_def.xspc = 1.0;
  test_im_def.yspc = 1.0;
  im_save(&test_im_def, "test_fft_def.dat", &status);

  im_fftw(&test_im_est, FFT_FFTW_ESTIMATE + FFT_X_ONLY + FWD + SCALE_N, &status);
  im_save(&test_im_est, "test_fft_est.dat", &status);

  im_fftw(&test_im_meas, FFT_FFTW_MEASURE + FFT_X_ONLY + FWD + SCALE_N, &status);
  im_save(&test_im_meas, "test_fft_meas.dat", &status);

  im_init(&diff_est, &status);
  im_init(&diff_meas, &status);

  im_sub(&test_im_est, &test_im_def, &diff_est, &status);
  im_sub(&test_im_meas, &test_im_def, &diff_meas, &status);
  
  im_save(&diff_est, "diff_fftw_est.dat", &status);
  im_save(&diff_meas, "diff_fftw_meas.dat", &status);

  im_destroy(&diff_est, &status);
  im_destroy(&diff_meas, &status);

  im_destroy(&test_im_def, &status);
  im_destroy(&test_im_est, &status);
  im_destroy(&test_im_meas, &status);

  printf("All done.\n");
  im_close_lib(&status);

  return 0;
}
