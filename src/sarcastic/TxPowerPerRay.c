/***************************************************************************
 *
 *       Module:    TxPowerPerRay.c
 *      Program:    Sadilac
 *   Created by:    Darren on 16/06/2013.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *       Function to calculate the power density of a given ray.
 *   The Transmit and receive antennae are defined in the header file.
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  5/04/2014
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
#include <sarclib/sarclib.h>
// #include "TxPowerPerRay.h"

// double TxPowerPerRay(double rayWidthRadians, double rayHeightRadians, double *receiverGain){
    
//     double txl,txh,rxl,rxh;
//     double gainTx, gainRx ;
//     double antAreaTx, antAreaRx ;
//     double Pray ;

//     if (TxDishAntenna) {
//         txh = txl = TxAntLen ;
//         antAreaTx = SIPC_pi * (txl/2) * (txl/2) ;
//     }else{
//         txh = TxAntHei ;
//         txl = TxAntLen ;
//         antAreaTx = txh * txl ;
//     }
//     gainTx   = 4.0 * SIPC_pi * antAreaTx * (ApEffTx/100.0) / (lambda * lambda);

//     if (monoStatic) {
// //        rxl = txl ;
// //        rxh = txh ;
//         gainRx = gainTx ;
//     }else{
//         if (RxDishAntenna) {
//             rxh = rxl = RxAntLen ;
//             antAreaRx = SIPC_pi * (rxl/2) * (rxl/2) ;
//         }else{
//             rxh = RxAntHei ;
//             rxl = RxAntLen ;
//             antAreaRx = rxh * rxl ;
//         }
//         gainRx   = 4.0 * SIPC_pi * antAreaRx * (ApEffRx/100.0) / (lambda * lambda);

//     }

// //    Pray = Pt * gainTx * rayWidthRadians * rayHeightRadians / (4.0*SIPC_pi) ;
//     Pray = Pt * gainTx ;
    
//     *receiverGain = gainRx ;
//     return (Pray) ;
   
// }

double TxPowerPerRay(CPHDHeader *hdr, double *receiverGain){
    double Pt, lambda, antLenTx, antHeiTx, antAreaTx, gainTx, gainRx, Pray, AppEff ;

    Pt        = 1000.0 ;       // Peak Transmit power assumed to be 1000 w (until its put in the CPHDFile)
    AppEff    = 0.7 ;          // Aperture efficiency - assumed to be 70 %
    lambda    = SIPC_c / hdr->freq_centre ;
    antLenTx  = lambda / hdr->antenna_width_az ;
    antHeiTx  = lambda / hdr->antenna_width_el ;
    antAreaTx = antLenTx * antHeiTx ;
    gainTx    = 4.0 * SIPC_pi * antAreaTx * AppEff / (lambda * lambda);
    gainRx    = gainTx ; // Assume transmit antenna is same as receive antenna
    Pray      = Pt * gainTx ;
    *receiverGain = gainRx ;
    return (Pray) ;
}
