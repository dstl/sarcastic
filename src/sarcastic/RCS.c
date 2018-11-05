/***************************************************************************
 * 
 *           Module :  RCS.c
 *          Program :  sarcastic
 *       Created by :  Darren Muff on 07/05/2015
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      Function to calculate the radar cross section (RCS) of a target
 *      at a given range and with a provided Efield magnitude
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

#include "RCS.h"
#include <sarclib/sarclib.h>

double RCS(double PtGt, double EMagAtRx, double TxRange, double RxRange){

    
    // first calculation the isotropically radiated power from a point at a range
    // from the transmitter
    //
    double iso = PtGt / (4 * SIPC_pi * TxRange * TxRange) ; // power in Watts / m^2
    
    // Correct the received radar field magnitude to be the power at the reflection point
    //
    double Ppt = EMagAtRx * EMagAtRx * ( 4.0 * SIPC_pi * RxRange * RxRange); // Watts

    // RCS is defined as the ratio of radar reflected power in the direction of a receiver
    // compared the the power reflected isotropically
    //
    double rcs = Ppt / iso ; // unitless
    
    return rcs;
}

