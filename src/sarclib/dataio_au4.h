/***************************************************************************
 * 
 *           Module :  dataio_au4.h
 *          Program :  sarclib
 *       Created by :  Darren Muff on Sat Jun 30 10:34:25 2018 +0100
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      Reads in the au4 file produced by etpm2au4 which is part of IFP4
 *      Originally written by Paul Eichel, Sandia National Labs
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

#ifndef sarclib_DATAIO_AU4_H__
#define sarclib_DATAIO_AU4_H__

#include "sarclib.h"

typedef struct
{
    int pulse_number;
    double tx_x;
    double tx_y;
    double tx_z;
    double rx_x;
    double rx_y;
    double rx_z;
    double time;
    double delta_r;
    double fscale;
    double c0;
    double c1;
    double c2;
} AuxInfo;

void read_au4(const char *fname, AuxInfo * au_info, int num_azi);

#endif

