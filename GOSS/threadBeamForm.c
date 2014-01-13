/***************************************************************************
 *
 *       Module:    threadBeamForm.c
 *      Program:    GOSS
 *   Created by:    Darren on 14/09/2013.
 *                  Copyright (c) 2013 Dstl. All rights reserved.
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
#include <stdio.h>
#include "GOSS.h"

void * beamForm ( void * threadArg ) {
    
    struct threadDataBF *td;
    td = (struct threadDataBF *) threadArg;
    
    int nx             = td->nx ;
    int nrnpItems      = td->nrnpItems;
    int startSamp      = td->startSamp;
    int nSamp          = td->nSamp ;
    rnpData_t *rnpData = td->rnpData ;
    int pulseIndex     = td->pulseIndex ;
    SPImage * phd      = td->phd ;
    double A           = td->A ;
    double B           = td->B ;
    
    double currentReal, currentImag, power, phse ;
    int samplingOffsetInt;

    // Create the uncompressed (but deramped) range profile for this pulse, storing it
    // back in its original location
    //
    for (int j=0; j<nrnpItems; j++ ){
        for (int ksamp=0; ksamp<nSamp; ksamp++) {
            samplingOffsetInt = rnpData[j].samplingOffsetInt ;
            if ((ksamp+samplingOffsetInt >=0) && (startSamp+ksamp+samplingOffsetInt < nx)) {
                currentReal = phd->data.cmpl_f[(pulseIndex)*phd->nx + (startSamp+ksamp+samplingOffsetInt)].r ;
                currentImag = phd->data.cmpl_f[(pulseIndex)*phd->nx + (startSamp+ksamp+samplingOffsetInt)].i ;
                
                power = rnpData[j].power ;
                // Reminder
                //  double A = 4.0*SIPC_pi*td->chirpRate / (SIPC_c * td->ADRate);
                //  double B = 4.0*SIPC_pi*td->oneOverLambda ;
                phse  = -1 * rnpData[j].rdiff * ( A * (startSamp+ksamp+rnpData[j].indexOffset) + B);
                td->phd->data.cmpl_f[(pulseIndex)*td->phd->nx + (startSamp+ksamp+samplingOffsetInt)].r = currentReal+power*cos(phse) ;
                td->phd->data.cmpl_f[(pulseIndex)*td->phd->nx + (startSamp+ksamp+samplingOffsetInt)].i = currentImag+power*sin(phse) ;
            }
        }
    }
    
    // return to parent thread
    //
    pthread_exit(NULL) ;
}
