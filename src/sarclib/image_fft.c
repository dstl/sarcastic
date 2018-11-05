/***************************************************************************
 * 
 *           Module :  image_fft.c
 *          Program :  sarclib
 *       Created by :  Emma Griffiths on 04/07/2005
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      This file contains functions that perform FFTs
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

#include "sarclib.h"

// The private (local) function that actualy calls the nightmare code that is fft_c8. 
// ptr is a pointer to the first element in an array that we need to FFT 
// xwidth is the width of the data. (nx) 
// ywidth is the height (number of rows) of the data (ny) 
// stride is the number of elements to skip between FFT's - always = xwdith (= nx) in our simple usage
//
static int fft_x_c8(SPCmplx *ptr, int xwidth, int ywidth, int stride, FFTMODE mode)
{
    int error = NO_ERROR;
    int j;
    
    // So, we do an FFT of a row of our image, hence we need to do ywidth (ny) of them
    //
    for (j=0; j < ywidth; j++)
    {
        // This calls the dreaded fft_c8 function. By this point, the
        // only bits of mode that it cares about are if its a forward or
        // reverse FFT. Due to the strangeness (flexiability) of fft_c8,
        // we need to pass it the number of elements in the data
        // it needs to FFT 3 times - I'm sure its API made sense to
        //  someone.
        //
        if((error=fft_c8(ptr, xwidth, xwidth, xwidth, ((mode & REV) ? 2 : -2) )) > NO_ERROR) break;
        
        // Increment the pointer by the xwidth of the data, so that ptr now points to the start
        // of the next row
        //
        ptr += stride;
    }
    return error;
}

// This is the private (local) FFT routine that is called by im_fft() after making sure the image
// is complex cartesian
//
static SPStatus* im_cart_fft(SPImage *a, FFTMODE mode, SPStatus * status)
{
    double num_scale = 1.0;
    CHECK_STATUS(status);
    
    // If the FFT_X_ONLY flag is set in mode (which is also set if FFT_2D is selected!) then do an FFT in
    // the X direction
    //
    if(mode & FFT_X_ONLY)
    {
        // Calculate the new spacing
        //
        a->xspc = (float)1.0/((float) a->nx *a->xspc);	// reciprocal units 
        
        // Call the private FFT function above 
        //
        status->status = fft_x_c8(a->data.cmpl_f, (int)a->nx, (int)a->ny, (int)a->nx, mode);
        
        // Store the number of elements in FFT for later scaling 
        //
        num_scale = a->nx;
        
        // Fall out if we've already created an error 
        //
        if (status->status > NO_ERROR)  return (status);
    }
    
    
    // If the FFT_Y_ONLY flag is set in mode (which is also set if FFT_2D is selected!) then do an FFT in
    // the Y direction
    //
    if(mode & FFT_Y_ONLY)
    {
        SPImage img_t;
        
        // The FFT routine we use can only do FFTs in the X direction, so we have to transpose the data,
        // do an FFT in the X direction, and then transpose it back again
        //
        // Calculate what the new spacing should be
        //
        a->yspc = (float)1.0/((float) a->ny *a->yspc);	// reciprocal units
        
        // Create a copy of the supplied data and transpose it 
        //
        im_init(&img_t, status);
        im_copy(a, &img_t, status);
        im_transp(&img_t, status);
        if(status->status == NO_ERROR)
        {
            // Clear away the original image 
            //
            im_destroy(a, status);
            
            // Perform the X-direction FFT on the transposed data 
            //
            status->status = fft_x_c8(img_t.data.cmpl_f, (int)img_t.nx, (int)img_t.ny, (int)img_t.nx, mode);
            if (status->status == NO_ERROR)
            {
                // Increment the number of operation count by the length of the FFT we just did 
                //
                num_scale *= img_t.nx;
                
                // Transpose the data back again 
                //
                im_transp(&img_t, status);
                
                // Copy it back to the original structure 
                //
                im_copy(&img_t, a, status);
            }
        }
        
        // Destroy the temporary image that contained the transposed data 
        //
        im_destroy(&img_t, status);
        if (status->status > NO_ERROR)  return (status);
    }
    
    // If we did an FFT (ie. either X_ONLY and/or Y_ONLY were set) and the user requested scaling
    // (i.e. either SCALE_N or SCALE_R were set), the scale the data appropriatly
    //
    if((mode & (FFT_X_ONLY+FFT_Y_ONLY)) && (mode & (SCALE_N+SCALE_R))) {		// need to scale
        double scale=1.0;
        
        // Calculate the scale factor from the number of FFT points and if its scale_R or _N 
        //
        scale/= (double) ((mode & SCALE_R) ? sqrt((double)num_scale) : num_scale);
        
        im_mult_scalar(a, ITYPE_DOUBLE, &scale, status);
    }
    
    return (status);
}

// This function is the main entry point for the user FFT 
// a  is the image to FFT, and the output from the FFT 
// mode is a bitfield that is created from or'ing (or adding ) together the fields
// specifed in fft.h, for example:
//  FWD + SCALE_N + FFT_Y_ONLY
//      or
//  REV + NOSCALE + FFT_2D
//
SPStatus* im_fft(SPImage *a, FFTMODE mode, SPStatus * status)                     // fft an image 
{
    SPImage b;
    
    CHECK_STATUS(status);
    CHECK_IMINTEG(a,status);
    
    // If we have a polar format type complex number, then convert it to cartesian, otherwise
    //just call the private FFT routine (im_cart_fft)
    //
    if (a->image_type == ITYPE_POLAR)
    {
        im_init(&b, status);
        // Convert to cart. 
        //
        im_polar_to_cart(a, &b, status);
        CHECK_STATUS(status);
        
        // Do the normal FFT 
        //
        im_cart_fft(&b, mode, status);
        
        CHECK_STATUS(status);
        
        // Convert back to polar 
        //
        im_cart_to_polar(&b, a, status);
        CHECK_STATUS(status);
        
        im_destroy(&b, status);
    }
    else if (a->image_type == ITYPE_CMPL_FLOAT)
    {
        im_cart_fft(a, mode, status);
    }
    else if (a->image_type == ITYPE_FLOAT)
    {
        SPImage tmp;
        im_create(&tmp, ITYPE_CMPL_FLOAT, a->nx, a->ny, a->xspc, a->yspc, status);
        CHECK_STATUS(status);
        im_cast(a, SET_REAL, &tmp, status);
        im_destroy(a, status);
        CHECK_STATUS(status);
        im_cart_fft(&tmp, mode, status);
        CHECK_STATUS(status);
        im_clone(&tmp, a, status);
    }
    else
    {
        printf("im_fft: Only ITYPE_FLOAT & ITYPE_CMPL_FLOAT are supported - data type is %s\n",
               itype2string(a->image_type));
        status->status = INVALID_TYPE;
    }
    return(status);
}
