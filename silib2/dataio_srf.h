/** @file********************************************************************
 *
 *       Module:    dataio_srf.h
 *      Program:    SIlib2
 *   Created by:    Matt Nottingham on 18/08/2008.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      Description: Files for reading/writing SRF format files.
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

#ifndef SILib2_DATAIO_SRF_H__
#define SILib2_DATAIO_SRF_H__

#define SRF_RECORD_LENGTH  (128)

typedef struct
{
    char datetime[SRF_RECORD_LENGTH + 1];       ///< Date and time stamp.
    char datum[16];                             ///< Geographic datum (e.g. WGS-84)
    FILE *  fp;                                 ///< The file pointer
    SPImageType  data_type;                     ///< Image data tyoe
    
    SPVector ll1;                               ///< Corner Coordinates; lat,long; UL
    SPVector ll2;                               ///< Corner Coordinates; lat,long; UR
    SPVector ll3;                               ///< Corner Coordinates; lat,long; LR
    SPVector ll4;                               ///< Corner Coordinates; lat,long; LL
    
    SPVector grp;                               ///< GRP (ECEF, m)
    
    SPVector ipn;                               ///< Image Plane Normal (ECEF, m)
    SPVector ip_rng;                            ///< Image Plane Range (ECEF, m)
    SPVector ip_cr;                             ///< Image Plane CrsRng (ECEF, m)
    SPVector fpn;                               ///< Focus Plane Normal (ECEF, m)
    SPVector apc;                               ///< Aperture Center (ECEF, m)
    SPVector apc_tx;                            ///< Aperture Center, bistatic xmit (ECEF, m)
    SPVector apc_rx;                            ///< Aperture Center, bistatic recr (ECEF, m)
    SPVector syn;                               ///< Synthetic Aperture (ECEF, m)
    SPVector vel;                               ///< Velocity (ECEF, m/s)
    double dof[2];                              ///< Depth of focus, (m)
    
    double lambda[2];                           ///< Effective Wavelength (m)
    int img_grp[2];                             ///< Image GRP pixel, X,Y
    int num_pulse;                              ///< Number of pulses
    int pulse_offset;                           ///< Pulse Offset
    int num_samp;                               ///< Number of samples
    int samp_offset;                            ///< Sample offset
    int fft_dim[2];                             ///< Dimensions of 2-D FFT, X,Y
    double sl_ratio[2];                         ///< Window function sidelobe ratio, dB
    int rect_dim[2];                            ///< Dimensions of dispersed image, X,Y
    double rect0[2];                            ///< Dispersed image center freq, X,Y (m^-1)
    double rect_sp[2];                          ///< Dispersed image sample spacing, X,Y (m^-1/samp)
    double pulse_duration;                      ///< Pulse duration (sec)
    double tx_bw;                               ///< tx BW (Hz)
    double rx_bw;                               ///< rx BW (Hz)
    
    int     nrows;                              ///< The number of rows in the range dimension (ny)
    int     ncols;                              ///< The number of elements per row (crossrange) (nx)
    double  sx;                                 ///< Image scale factor in X (m/pixel)
    double  sy;                                 ///< Image scale factor in Y (m/pixel)
    double  resX;                               ///< Image resolution in X (m)
    double  resY;                               ///< Image resolution in Y (m)
    int     do_byte_swap;                       ///< Whether we need to swap bytes due to little/big endian issues
    off_t  data_start;                          ///< Position in the file where the actual data starts
} SRFheader;

SPStatus * im_open_srf(SRFheader * hdr, const char * filename, SPStatus * status);
SPStatus * im_load_srf_subset(SPImage * data, SRFheader * hdr, int startx, int starty, SPStatus * status);
SPStatus * im_load_srf(SPImage * data, SRFheader * hdr, SPStatus * status);
SPStatus * im_close_srf(SRFheader * hdr, SPStatus * status);
SPStatus * im_destroy_srf(SRFheader * hdr, SPStatus * status);

#endif
