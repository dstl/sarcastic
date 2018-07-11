/***************************************************************************
 *
 *       Module:    tdpKernel.h
 *      Program:    tdpocl
 *   Created by:    Darren on 21/03/2013.
 *                  Copyright (c) 2013 [Dstl]. All rights reserved.
 *
 *   Description:
 *   <ENTER DESCRIPTION HERE>
 *
 *
 *   CLASSIFICATION        :  <PENDING>
 *   Date of CLASSN        :  14/03/2013
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

#ifndef tdpocl_tdpKernel_h
#define tdpocl_tdpKernel_h

static const char *kernelCode     = " \n \
 \n \
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n \
#define CMPLX_MULT(a,b,out) {out.r = a.r * b.r - a.i * b.i; out.i = a.i * b.r + a.r * b.i;} \n \
#define CMPLX_ADD(a,b,out) {out.r = a.r + b.r; out .i = a.i + b.i;} \n \
#define VECTMAG(v,m){m=sqrt(v.x*v.x+v.y*v.y+v.z*v.z);} \n \
 \n \
typedef struct { double x, y, z; } VectorH ; \n \
 \n \
typedef struct CMPLX { \n \
    float r; \n \
    float i; \n \
} CMPLX ; \n \
 \n \
void gpu_sinc_interp(__global CMPLX *f, CMPLX *g, float rangeLabel, __global float * ikernel) ; \n \
 \n \
__kernel void tdpocl(__global CMPLX *dev_Data,             // Raw pulse data. x range compressed samples; y pulse index \n \
                     __global CMPLX *dev_Image,            // Output image. \n \
                     const int nx, const int ny,           // image size in x & y \n \
                     __global VectorH *dev_satPos_rx,      // Sensor receiver positions \n \
                     __global VectorH *dev_satPos_tx,      // Sensor transmitter positions \n \
                     __global double *dev_srpRanges,       // Range from pulse to SRP and back \n \
                     __global VectorH *dev_surf,           // Surface array (2D) defining location of pixels in output image \n \
                     const int dev_phaseSgn,               // +1 or -1 : defines direction of increasing range in receiver IQ \n \
                     __global float *dev_fx_step,          // Frequency step per pulse sample for a given pulse \n \
                     const int dev_range_n,                // number of range samples in dev_Data \n \
                     const double dev_freq_centre,         // centre frequency of transmitted waveform \n \
                     __global float * dev_ikernel,         // Sinc interpolation kernel array of (OVERSAMP(512) * NPOINTS(8) +1) length \n \
                     const int dev_nPulses)                // Number of pulses in dev_Data \n \
{ \n \
 \n \
    int pulse; \n \
    double rng_to_srp, rng_to_pxl, range_diff; \n \
    double r1,r2; \n \
    double sx; \n \
    float phs; \n \
    double CC = 299792458.0 ; \n \
    double MPI = 3.141592653589793 ; \n \
    float rangeLabel; \n \
    VectorH diff1; \n \
    VectorH diff2; \n \
    VectorH rxp, txp, gpt; \n \
    double fxp; \n \
    CMPLX pcorr; \n \
    CMPLX intData; \n \
    CMPLX tmp; \n \
 \n \
    const int pxlx = get_global_id(0); \n \
    const int pxly = get_global_id(1); \n \
    const int pxl  = pxly * get_global_size(0) + pxlx ; \n \
 \n \
    if(pxlx>=0 && pxlx<nx && pxly>=0 && pxly<ny){ \n \
         \n \
        const int rgn = dev_range_n ; \n \
         \n \
        for (pulse = 0; pulse < dev_nPulses; pulse++){ \n \
             \n \
            intData.r = 0; \n \
            intData.i = 0; \n \
            fxp = dev_fx_step[pulse] ; \n \
             \n \
            rxp = dev_satPos_rx[pulse] ; \n \
            txp = dev_satPos_tx[pulse] ; \n \
            gpt = dev_surf[pxl]; \n \
             \n \
            rng_to_srp = dev_srpRanges[pulse]; \n \
             \n \
//            VECTSUB(rxp,gpt,diff1); \n \
            diff1.x = rxp.x - gpt.x ; \n \
            diff1.y = rxp.y - gpt.y ; \n \
            diff1.z = rxp.z - gpt.z ; \n \
             \n \
//            VECTSUB(txp,gpt,diff2); \n \
            diff2.x = txp.x - gpt.x ; \n \
            diff2.y = txp.y - gpt.y ; \n \
            diff2.z = txp.z - gpt.z ; \n \
             \n \
            VECTMAG(diff1,r1); \n \
            VECTMAG(diff2,r2); \n \
             \n \
            rng_to_pxl = r1+r2; \n \
             \n \
            range_diff = dev_phaseSgn * (rng_to_srp - rng_to_pxl) * 0.5 ; \n \
            sx = CC/(2 * fxp * rgn); \n \
            rangeLabel = range_diff/sx + (rgn / 2); \n \
             \n \
            // Calculating phase corrections: \n \
            // \n \
            phs = 4.0 * range_diff * MPI * dev_freq_centre / CC; \n \
            pcorr.i = sincos(phs, &(pcorr.r)); \n \
             \n \
            if (rangeLabel > 4 && rangeLabel < dev_range_n - 8) { \n \
                // Interpolate onto range profile: \n \
                gpu_sinc_interp(&dev_Data[pulse * rgn], &intData, rangeLabel, dev_ikernel); \n \
                 \n \
                // Applying phase correction to interpolated data \n \
                CMPLX_MULT(pcorr,intData,tmp); \n \
                CMPLX_ADD(tmp,dev_Image[pxl],dev_Image[pxl]); \n \
                 \n \
            } \n \
        } \n \
         \n \
    } \n \
} \n \
 \n \
void gpu_sinc_interp(__global CMPLX *f, CMPLX *g, float rangeLabel, __global float * ikernel) \n \
{ \n \
    int iidx; \n \
     \n \
    int offset; \n \
    const int nint_rl = (int) (rangeLabel+0.5);  // We assume that rangeLabel > 0 \n \
     \n \
    iidx = (int)((nint_rl-rangeLabel) * 512); \n \
    if (iidx < 0) { \n \
        iidx += 512; \n \
        offset = 1; \n \
    } else { \n \
        offset = 0; \n \
    } \n \
     \n \
    // Assume NPOINTS is 8 and loop unroll \n \
    // Assume OVERSAMP is 512; \n \
    // \n \
 \n \
    const float4 r_tmpL = (float4)(f[nint_rl-4+offset].r * ikernel[iidx], \n \
                                   f[nint_rl-3+offset].r * ikernel[iidx+512], \n \
                                   f[nint_rl-2+offset].r * ikernel[iidx+1024], \n \
                                   f[nint_rl-1+offset].r * ikernel[iidx+1536]); \n \
    const float4 r_tmpR = (float4)(f[nint_rl+offset].r * ikernel[iidx+2048], \n \
                                   f[nint_rl+1+offset].r * ikernel[iidx+2560], \n \
                                   f[nint_rl+2+offset].r * ikernel[iidx+3072], \n \
                                   f[nint_rl+3+offset].r * ikernel[iidx+3584]); \n \
     \n \
    const float4 r_tmp = r_tmpL + r_tmpR ; \n \
     \n \
    const float4 i_tmpL = (float4)(f[nint_rl-4+offset].i * ikernel[iidx], \n \
                                   f[nint_rl-3+offset].i * ikernel[iidx+512], \n \
                                   f[nint_rl-2+offset].i * ikernel[iidx+1024], \n \
                                   f[nint_rl-1+offset].i * ikernel[iidx+1536]); \n \
    const float4 i_tmpR = (float4)(f[nint_rl+offset].i * ikernel[iidx+2048], \n \
                                   f[nint_rl+1+offset].i * ikernel[iidx+2560], \n \
                                   f[nint_rl+2+offset].i * ikernel[iidx+3072], \n \
                                   f[nint_rl+3+offset].i * ikernel[iidx+3584]); \n \
     \n \
    const float4 i_tmp = i_tmpL + i_tmpR ; \n \
     \n \
    const float4 ones = (float4)(1.0f,1.0f,1.0f,1.0f); \n \
     \n \
    g->r = dot(r_tmp, ones) ; \n \
    g->i = dot(i_tmp, ones) ; \n \
} \n \
 \n \
" ;

#endif
