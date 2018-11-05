/***************************************************************************
 * 
 *           Module :  tdpKernel.cl
 *          Program :  kernels
 *       Created by :  Darren Muff on 21/03/2013
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :  OpenCL kernel for the time domain SAR processor.
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

#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define CMPLX_MULT(a,b,out) {out.r = a.r * b.r - a.i * b.i; out.i = a.i * b.r + a.r * b.i;}
#define CMPLX_ADD(a,b,out) {out.r = a.r + b.r; out .i = a.i + b.i;}
#define VECTCREATE(a,b,c,vOut){vOut.x=a; vOut.y=b; vOut.y=c;}
#define VECTSUB(a,b,vOut){vOut.x=a.x-b.x; vOut.y=a.y-b.y; vOut.z=a.z-b.z;}
#define VECTMAG(v,m){m=sqrt(v.x*v.x+v.y*v.y+v.z*v.z);}
#define NPOINTS (8)
#define OVERSAMP (512)

typedef struct { double x, y, z; } VectorH ;

typedef struct CMPLX {
    float r;
    float i;
} CMPLX ;
    
void gpu_sinc_interp(__global CMPLX *f, CMPLX *g, float rangeLabel, __global float * ikernel) ;

__kernel void tdpocl(__global CMPLX *dev_Data,             // Raw pulse data. x range compressed samples; y pulse index
                     __global CMPLX *dev_Image,            // Output image.
                     const int nx, const int ny,           // image size in x & y
                     __global VectorH *dev_satPos_rx,      // Sensor receiver positions
                     __global VectorH *dev_satPos_tx,      // Sensor transmitter positions
                     __global double *dev_srpRanges,       // Range from pulse to SRP and back
                     __global VectorH *dev_surf,           // Surface array (2D) defining location of pixels in output image
                     const int dev_phaseSgn,               // +1 or -1 : defines direction of increasing range in receiver IQ
                     __global float *dev_fx_step,          // Frequency step per pulse sample for a given pulse
                     const int dev_range_n,                // number of range samples in dev_Data
                     const double dev_freq_centre,         // centre frequency of transmitted waveform
                     __global float * dev_ikernel,         // Sinc interpolation kernel array of (OVERSAMP * NPOINTS +1) length
                     const int dev_nPulses)                // Number of pulses in dev_Data
{
    
    int pulse;
    double rng_to_srp, rng_to_pxl, range_diff;
    double r1,r2;
    double sx;
    float phs;
    double CC = 299792458.0 ;
    double MPI = 3.141592653589793 ;
    float rangeLabel;
    VectorH diff1;
    VectorH diff2;
    VectorH rxp, txp, gpt;
    double fxp;
    CMPLX pcorr;
    CMPLX intData;
    CMPLX tmp;
    
    const int pxlx = get_global_id(0);
    const int pxly = get_global_id(1);
    const int pxl  = pxly * get_global_size(0) + pxlx ;
    
    if(pxlx>=0 && pxlx<nx && pxly>=0 && pxly<ny){

        const int rgn = dev_range_n ;

        for (pulse = 0; pulse < dev_nPulses; pulse++){

            intData.r = 0;
            intData.i = 0;
            fxp = dev_fx_step[pulse] ;
            
            // Unpack variables from global scope to local registers
            // Required for some GPU cards (eg ATI Radeon HD 5870)
            // Only required for vectorH calcs as their vars do not specify memory scope
            //
            rxp = dev_satPos_rx[pulse] ;
            txp = dev_satPos_tx[pulse] ;
            gpt = dev_surf[pxl];
           
            rng_to_srp = dev_srpRanges[pulse];
            
//            VECTSUB(rxp,gpt,diff1);
            diff1.x = rxp.x - gpt.x ;
            diff1.y = rxp.y - gpt.y ;
            diff1.z = rxp.z - gpt.z ;
            
//            VECTSUB(txp,gpt,diff2);
            diff2.x = txp.x - gpt.x ;
            diff2.y = txp.y - gpt.y ;
            diff2.z = txp.z - gpt.z ;
            
            VECTMAG(diff1,r1);
            VECTMAG(diff2,r2);

            rng_to_pxl = r1+r2;

            range_diff = dev_phaseSgn * (rng_to_srp - rng_to_pxl) * 0.5 ;
            sx = CC/(2 * fxp * rgn);
            rangeLabel = range_diff/sx + (rgn / 2);
            
            // Calculating phase corrections:
            //
            phs = 4.0 * range_diff * MPI * dev_freq_centre / CC;
            pcorr.i = sincos(phs, &(pcorr.r));
                      
            if (rangeLabel > NPOINTS/2 && rangeLabel < dev_range_n - NPOINTS) {
                // Interpolate onto range profile:
                gpu_sinc_interp(&dev_Data[pulse * rgn], &intData, rangeLabel, dev_ikernel);
                
                // Applying phase correction to interpolated data
                CMPLX_MULT(pcorr,intData,tmp);
                CMPLX_ADD(tmp,dev_Image[pxl],dev_Image[pxl]);

            }
        }

    }
}

void gpu_sinc_interp(__global CMPLX *f, CMPLX *g, float rangeLabel, __global float * ikernel)
{
//    float i_tmp = 0.0;
//    float r_tmp = 0.0;
    int iidx;
//    int i;
    
    int offset;
    const int nint_rl = (int) (rangeLabel+0.5);  // We assume that rangeLabel > 0
    
    iidx = (int)((nint_rl-rangeLabel) * OVERSAMP);
    if (iidx < 0) {
        iidx += OVERSAMP;
        offset = 1;
    } else {
        offset = 0;
    }
    
    // Assume NPOINTS is 8 and loop unroll
    //
    // Assume OVERSAMP is 512;
    
    const float4 r_tmpL = (float4)(f[nint_rl-4+offset].r * ikernel[iidx],
                    f[nint_rl-3+offset].r * ikernel[iidx+512],
                    f[nint_rl-2+offset].r * ikernel[iidx+1024],
                    f[nint_rl-1+offset].r * ikernel[iidx+1536]);
    const float4 r_tmpR = (float4)(f[nint_rl+offset].r * ikernel[iidx+2048],
                     f[nint_rl+1+offset].r * ikernel[iidx+2560],
                     f[nint_rl+2+offset].r * ikernel[iidx+3072],
                     f[nint_rl+3+offset].r * ikernel[iidx+3584]);
    
    const float4 r_tmp = r_tmpL + r_tmpR ;
    
    const float4 i_tmpL = (float4)(f[nint_rl-NPOINTS/2+offset].i * ikernel[iidx],
                     f[nint_rl-3+offset].i * ikernel[iidx+512],
                     f[nint_rl-2+offset].i * ikernel[iidx+1024],
                     f[nint_rl-1+offset].i * ikernel[iidx+1536]);
    const float4 i_tmpR = (float4)(f[nint_rl+offset].i * ikernel[iidx+2048],
                     f[nint_rl+1+offset].i * ikernel[iidx+2560],
                     f[nint_rl+2+offset].i * ikernel[iidx+3072],
                     f[nint_rl+3+offset].i * ikernel[iidx+3584]);
    
    const float4 i_tmp = i_tmpL + i_tmpR ;


//    for(i = 0; i < NPOINTS; i++) {
//
//        r_tmp += f[nint_rl+i-NPOINTS/2+offset].r * ikernel[iidx];
//        i_tmp += f[nint_rl+i-NPOINTS/2+offset].i * ikernel[iidx];
//        
//        iidx += OVERSAMP;
//    }
    
//    g->r = r_tmp;
//    g->i = i_tmp;
    
    
    const float4 ones = (float4)(1.0f,1.0f,1.0f,1.0f);

    g->r = dot(r_tmp, ones) ;
    g->i = dot(i_tmp, ones) ;
}

