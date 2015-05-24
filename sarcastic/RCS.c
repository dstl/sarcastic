/***************************************************************************
 *
 *       Module:    RCS.c
 *      Program:    SARCASTIC
 *   Created by:    Darren on 07/05/2015.
 *                  Copyright (c) 2013 Dstl. All rights reserved.
 *
 *   Description:
 *      Function to calculate the radar cross section (RCS) of a target
 *      at a given range and with a provided Efield magnitude
 *
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  7/05/2015
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

#include "RCS.h"
#include <SIlib/SIlib.h>

double RCS(double PtGt, double EMagAtRx, double TxRange, double RxRange){

    
    // first calculation the isotropically radiated power from a point at a range
    // from the transitter
    //
    double iso = PtGt / (4 * SIPC_pi * TxRange * TxRange) ; // power in Watts / m^2
    
    // Correct the received radar field magnitude to be the power at the reflection point
    // An interesting point here. The E field measured at the receiver is the mean
    // amplitude of a sinusoidal wave function. The peak amplitude is therefore twice as large
    // (hence the factor of 4 below when calculating the power (amplitude squared)
    //
    double Ppt = EMagAtRx * EMagAtRx * 4 * ( 4.0 * SIPC_pi * RxRange * RxRange); // Watts / m^2

    // RCS is defined as the ratio of radar reflected power in the direction of a receiver
    // compared the the power reflected isotropically
    //
    double rcs = Ppt / iso ; // unitless
    
    return rcs;
}

