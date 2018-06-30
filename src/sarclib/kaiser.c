/***************************************************************************
 *
 *       Module:    kaiser.c
 *      Program:    sarclib
 *   Created by:    Emma Griffiths on 10/08/2005.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This is a short piece of code for creating a kaiser-bessel window
 *      function to use for things like the sinc interpolation.
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

// alpha is the windowing parameter
//
SPStatus * kaiser (SPImage *in, SPImage *window, double alpha, SPStatus *status)
{
    CHECK_STATUS(status);
    
    int i;
    
    double tmp;
    
    SPImage temp;
    SPImage tempc;
    SPImage I0a;
    SPImage I0b;
    
    im_create(&temp, ITYPE_DOUBLE, in->nx, in->ny, in->xspc, in->yspc, status);
    im_create(&I0a, ITYPE_DOUBLE, in->nx, in->ny, in->xspc, in->yspc, status);
    
    im_create(&tempc, ITYPE_DOUBLE, 1, in->ny, in->xspc, in->yspc, status);
    im_create(&I0b, ITYPE_DOUBLE, 1, in->ny, in->xspc, in->yspc, status);
    
    for (i = 0; i < in->nx; i++)
    {
        // calculate the input for the bessel fn
        //
        temp.data.d[i] = alpha*M_PI*sqrt(1.0 - SP_SQR((in->data.d[i]/in->data.d[in->nx-1]))); 
        if (status->debug >= 30)
        {
            printf("\n%s:%d\n", __FILE__, __LINE__);
            printf("value of input to bessel fn %g alpha %g in %g\n", temp.data.d[i], alpha,in->data.d[i]);
        }
    }
    
    // alpha * pi is also a bessel function input, the kaiser equation is = I0( alpha*pi (1-(f(n)/f(n-1))^2)^(1/2) ) / IO(alpha*pi)
    //
    tempc.data.d[0] = alpha*M_PI;  
    
    I0(&temp, &I0a, status);
    I0(&tempc, &I0b, status);
    
    tmp = 1.0/I0b.data.d[0];
    
    im_mult_scalar(&I0a, ITYPE_DOUBLE, &tmp, status);
    im_copy (&I0a, window, status);
    
    im_destroy(&temp, status);
    im_destroy(&tempc, status);
    im_destroy(&I0a, status);
    im_destroy(&I0b, status);
    
    if (status->debug >= 30)
    {
        for (i = 0; i < window->nx*window->ny; i++)
        {
            printf("input to k-b fn %g output of k-b fn %g\n", in->data.d[i], window->data.d[i]);
        }
    }
    
    return(status);
}


SPStatus* I0 (SPImage *in, SPImage *out, SPStatus *status)
{
    // zero order modified Bessel function, I don't really know enough about them to add oodles of comments
    //
    int k, i;
    
    double fact;
    double p;
    double t;
    double eps = 1.0e-6;   // accuracy parameter 
    
    if (status->debug >= 30)
    {
        printf("\n%s:%d\n", __FILE__, __LINE__);
    }
    
    for (i = 0 ; i < in->nx * in->ny; i++)
    {
        fact = 1.0;
        p = in->data.d[i] * 0.5;
        t = SP_SQR(p);
        out->data.d[i] = 1.0 + t;
        
        for (k = 2; t > eps; k++)
        {
            p *= in->data.d[i] * 0.5;
            fact *= k;
            t = SP_SQR(p/fact);
            out->data.d[i] += t;
        }
        
        if (status->debug >= 30)
        {
            printf("Output of bessel function %g\n", out->data.d[i]);
        }
    }
    return(status);
}

