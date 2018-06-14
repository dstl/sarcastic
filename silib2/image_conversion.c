/***************************************************************************
 *
 *       Module:    image_conversion.c
 *      Program:    SIlib2
 *   Created by:    Emma Griffiths on 30/06/2005.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains functions needed to operate on complex polar arrays.
 *
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

#include "SIlib2.h"

#define CAST_SCALER( type ) {                                                                               \
    switch (dst->image_type) {                                                                              \
        case ITYPE_CMPL_FLOAT:                                                                              \
        {                                                                                                   \
            switch (mode) {                                                                                 \
                                                                                                            \
                /* casting from float containing phase info to output the phase as a cartesian image */     \
                case SET_REAL:                                                                              \
                {                                                                                           \
                    for (y = 0; y < src->ny; y++) {                                                         \
                        for (x = 0; x < src->nx; x++) {                                                     \
                            dst->data.cmpl_f[x+src->nx * y].r = src->data. type [x + src->nx * y];          \
                            dst->data.cmpl_f[x+src->nx * y].i = 0.0;                                        \
                        }                                                                                   \
                    }                                                                                       \
                }                                                                                           \
                break;                                                                                      \
                                                                                                            \
                /* casting from float containing phase info to output the phase as a cartesian image */     \
                case SET_IMAG:                                                                              \
                {                                                                                           \
                    for (y = 0; y < src->ny; y++) {                                                         \
                        for (x = 0; x < src->nx; x++) {                                                     \
                            dst->data.cmpl_f[x+src->nx * y].r = 0.0;                                        \
                            dst->data.cmpl_f[x+src->nx * y].i = src->data. type [x + src->nx * y];          \
                        }                                                                                   \
                    }                                                                                       \
                }                                                                                           \
                break;                                                                                      \
                                                                                                            \
                default:                                                                                    \
                    fprintf(stderr, "Cast mode (%s) not supported in im_cast\n", castmode2string(mode));    \
                    status->status = INVALID_TYPE;                                                          \
                    break;                                                                                  \
            }                                                                                               \
            break;                                                                                          \
            default:                                                                                        \
                fprintf(stderr, "Cast mode (%s) not supported in im_cast\n", castmode2string(mode));        \
                status->status = INVALID_TYPE;                                                              \
                break;                                                                                      \
        }                                                                                                   \
    }                                                                                                       \
}

// Converts a cartesian type complex number image to a polar format number image.
// Will probably be deprecated and its functionality moved into im_cast.
//
SPStatus* im_cart_to_polar (SPImage *a, SPImage *out, SPStatus *status)                      
{
    int64_t n;
    CHECK_STATUS(status);
    CHECK_PTR(out,status);
    CHECK_IMINTEG(a,status);
    CHECK_TYPE(a, ITYPE_CMPL_FLOAT, status);
    
    if ((out->data.v == NULL) && (out->nx == 0) && (out->ny ==0))
    {
        im_create (out, ITYPE_POLAR, a->nx, a->ny, a->xspc, a->yspc, status);
    }
    
    CHECK_IMSIZEOUT(a, out, status);
    CHECK_TYPE(out, ITYPE_POLAR, status);
    
    for (n=0; n<(a->nx*a->ny); n++)
    {
        CMPLX_CART2POL(a->data.cmpl_f[n], out->data.pol[n]);
    }
    return(status);
}

// Converts a polar type complex number image to a cartesian format number image.
// Will probably be deprecated and its functionality moved into im_cast.
//
SPStatus* im_polar_to_cart (SPImage *a, SPImage *out, SPStatus *status) 
{
    int64_t n;
    CHECK_STATUS(status);
    CHECK_PTR(out,status);
    CHECK_IMINTEG(a,status);
    CHECK_TYPE(a, ITYPE_POLAR, status);
    
    if ((out->data.v == NULL) && (out->nx == 0) && (out->ny ==0))
    {
        im_create (out, ITYPE_CMPL_FLOAT, a->nx, a->ny, a->xspc, a->yspc, status);
    }
    
    CHECK_IMSIZEOUT(a, out, status);
    CHECK_TYPE(out, ITYPE_CMPL_FLOAT, status);
    
    for (n=0; n<(a->nx*a->ny); n++)
    {
        CMPLX_POL2CART(a->data.pol[n], out->data.cmpl_f[n]);
    }
    return(status);
}

// Converts an image of one type into another type of image. The exact mode of
// the trasnfer (eg, when going from complex to float - do you want the magnitude
// or the phase)
//
SPStatus* im_cast(SPImage *src, SPCastMode mode, SPImage *dst, SPStatus * status)
{
    int64_t x;
    int64_t y;
    
    CHECK_STATUS(status);
    CHECK_PTR(src,status);
    CHECK_PTR(dst,status);
    CHECK_IMCONT(src,status);
    CHECK_IMINTEG(src,status);
    
    CHECK_TYPE_VALID(dst, status);
    
    if (dst->data.v == NULL)
    {
        im_create(dst, dst->image_type, src->nx, src->ny, src->xspc, src->yspc, status);
    }
    else
    {
        CHECK_IMCONT(dst,status);
        CHECK_IMINTEG(dst,status);
    }
    CHECK_IMSIZEOUT(src, dst, status);
    
    if (status->debug >= 20)
    {
        printf("Casting from %s to %s\n", itype2string(src->image_type), itype2string(dst->image_type));
    }
    
    switch(src->image_type)
    {
        case ITYPE_CMPL_FLOAT:
        {
            switch (dst->image_type)
            {
                case ITYPE_FLOAT:
                {
                    switch (mode)
                    {
                        case CAST_MODULUS:  // casting from complex float to output the modulus as a float 
                            for(y = 0; y < src->ny; y++)
                            {
                                for(x = 0; x < src->nx; x++)
                                {
                                    dst->data.f[x + src->nx * y] = CMPLX_MAG(src->data.cmpl_f[x + src->nx * y]);
                                }
                            }
                            break;
                            
                        case CAST_PHASE:  // casting from complex float to output the phase as a float
                            for(y = 0; y < src->ny; y++)
                            {
                                for(x = 0; x < src->nx; x++)
                                {
                                    dst->data.f[x + src->nx * y] = CMPLX_PHASE(src->data.cmpl_f[x + src->nx * y]);
                                }
                            }
                            break;
                            
                        case CAST_REAL:  // casting from complex float to output the real as a float 
                            for(y = 0; y < src->ny; y++)
                            {
                                for(x = 0; x < src->nx; x++)
                                {
                                    dst->data.f[x + src->nx * y] = src->data.cmpl_f[x + src->nx * y].r;
                                }
                            }
                            break;
                            
                        case CAST_IMAG:  // casting from complex float to output the real as a float
                            for(y = 0; y < src->ny; y++)
                            {
                                for(x = 0; x < src->nx; x++)
                                {
                                    dst->data.f[x + src->nx * y] = src->data.cmpl_f[x + src->nx * y].i;
                                }
                            }
                            break;
                            
                        default:
                            fprintf(stderr, "Cast mode (%s) not supported in im_cast\n", castmode2string(mode));
                            status->status = INVALID_TYPE;
                            break;
                    }
                }
                    break;
                default:
                    fprintf(stderr, "Dest image type (%s) not supported in im_cast\n", itype2string(dst->image_type));
                    status->status = INVALID_TYPE;
                    break;
            }
        }
            break;
            
        case ITYPE_UINT16:
        {
            CAST_SCALER( ui16 );
        }
            break;
        case ITYPE_UINT8:
        {
            CAST_SCALER( ui8 );
        }
            break;
        case ITYPE_FLOAT:
        {
            switch (dst->image_type)
            {
                case ITYPE_POLAR:
                {
                    switch (mode)
                    {
                        case SET_MODULUS:  // casting from float containing modulus info to output the modulus as a polar image with no phase
                            for (y = 0; y < src->ny; y++)
                            {
                                for (x = 0; x < src->nx; x++)
                                {
                                    dst->data.pol[x+src->nx * y].pabs = src->data.f[x + src->nx * y];
                                    dst->data.pol[x+src->nx * y].parg = 0.0;
                                }
                            }
                            break;
                            
                        case SET_PHASE:  // casting from float containing phase info to output the phase as a polar image with magnitude 1
                            for (y = 0; y < src->ny; y++)
                            {
                                for (x = 0; x < src->nx; x++)
                                {
                                    dst->data.pol[x+src->nx * y].pabs = 1.0;
                                    dst->data.pol[x+src->nx * y].parg = src->data.f[x + src->nx * y];
                                }
                            }
                            break;
                            
                        default:
                            fprintf(stderr, "Cast mode (%s) not supported in im_cast\n", castmode2string(mode));
                            status->status = INVALID_TYPE;
                            break;
                    }
                }
                    break;
                    
                case ITYPE_CMPL_FLOAT:
                {
                    switch (mode)
                    {
                        case SET_PHASE:  // casting from float containing phase info to output the phase as a cartesian image 
                        {
                            for (y = 0; y < src->ny; y++)
                            {
                                for (x = 0; x < src->nx; x++)
                                {
                                    dst->data.cmpl_f[x+src->nx * y].r = cos(src->data.f[x + src->nx * y]);
                                    dst->data.cmpl_f[x+src->nx * y].i = sin(src->data.f[x + src->nx * y]);
                                }
                            }
                        }
                            break;
                            
                        case SET_REAL:  // casting from float containing phase info to output the phase as a cartesian image 
                        {
                            for (y = 0; y < src->ny; y++)
                            {
                                for (x = 0; x < src->nx; x++)
                                {
                                    dst->data.cmpl_f[x+src->nx * y].r = src->data.f[x + src->nx * y];
                                    dst->data.cmpl_f[x+src->nx * y].i = 0.0;
                                }
                            }
                        }
                            break;
                            
                        case SET_IMAG:  // casting from float containing phase info to output the phase as a cartesian image
                        {
                            for (y = 0; y < src->ny; y++)
                            {
                                for (x = 0; x < src->nx; x++)
                                {
                                    dst->data.cmpl_f[x+src->nx * y].r = 0.0;
                                    dst->data.cmpl_f[x+src->nx * y].i = src->data.f[x + src->nx * y];
                                }
                            }
                        }
                            break;
                            
                        default:
                            fprintf(stderr, "Cast mode (%s) not supported in im_cast\n", castmode2string(mode));
                            status->status = INVALID_TYPE;
                            break;
                    }
                }
                    break;
                    
                default:
                    fprintf(stderr, "Dest image type (%s) not supported in im_cast\n", itype2string(dst->image_type));
                    status->status = INVALID_TYPE;
                    break;
            }
            break;
        }
            
        case ITYPE_VECTOR:
        {
            switch (dst->image_type)
            {
                case ITYPE_VECTOR:
                {
                    switch (mode)
                    {
                        case CAST_UNIT:
                        {
                            for (y = 0; y < src->ny; y++)
                            {
                                for (x = 0; x < src->nx; x++)
                                {
                                    VECT_UNIT(src->data.vect[x+src->nx*y], dst->data.vect[x+src->nx*y]);
                                }
                            }
                        }
                            break;
                        default:
                            fprintf(stderr, "Cast mode (%s) not supported in im_cast\n", castmode2string(mode));
                            status->status = INVALID_TYPE;
                            break;
                    }
                }
                    break;
                case ITYPE_DOUBLE:
                {
                    switch (mode)
                    {
                        case CAST_MODULUS:
                        {
                            for (y = 0; y < src->ny; y++)
                            {
                                for (x = 0; x < src->nx; x++)
                                {
                                    dst->data.d[x+src->nx*y] = VECT_MOD(src->data.vect[x+src->nx*y]);
                                }
                            }
                        }
                            break;
                        default:
                            fprintf(stderr, "Cast mode (%s) not supported in im_cast\n", castmode2string(mode));
                            status->status = INVALID_TYPE;
                            break;
                    }
                }
                    break;
                case ITYPE_FLOAT:
                {
                    switch (mode)
                    {
                        case CAST_MODULUS:
                        {
                            for (y = 0; y < src->ny; y++)
                            {
                                for (x = 0; x < src->nx; x++)
                                {
                                    dst->data.f[x+src->nx*y] = VECT_MOD(src->data.vect[x+src->nx*y]);
                                }
                            }
                        }
                            break;
                        default:
                            fprintf(stderr, "Cast mode (%s) not supported in im_cast\n", castmode2string(mode));
                            status->status = INVALID_TYPE;
                            break;
                    }
                }
                    break;
                default:
                    fprintf(stderr, "Dest image type (%s) not supported in im_cast\n", itype2string(dst->image_type));
                    status->status = INVALID_TYPE;
                    break;
            }
            break;
        }
            
        default:
            fprintf(stderr, "Src image type (%s) not supported in im_cast\n", itype2string(src->image_type));
            status->status = INVALID_TYPE;
            break;
    }
    
    return (status);
}


// Converts two images (re & im) to a complex float
//
SPStatus *im_cmplx_f_make(SPImage *real, SPImage *imaginary, SPImage *out, SPStatus *status)
{
    int x, y;
    
    CHECK_IMINTEG(real, status);
    CHECK_IMINTEG(imaginary, status);
    CHECK_IMSIZEIN(real, imaginary, status);
    CHECK_TYPES_EQUAL(real, imaginary, status);
    if ((out->data.v == NULL) && (out->nx == 0) && (out->ny ==0))
    {
        im_create (out, ITYPE_CMPL_FLOAT, real->nx, real->ny, real->xspc, real->yspc, status);
    }
    
    CHECK_IMSIZEOUT(real, out, status);
    CHECK_TYPE(out, ITYPE_CMPL_FLOAT, status);
    CHECK_STATUS(status);
    
    switch(real->image_type)
    {
        case ITYPE_DOUBLE:
            for (y = 0; y < real->ny; y++)
            {
                for (x = 0; x < real->nx; x++)
                {
                    CMPLX_F_MAKE(real->data.d[x + y*real->nx], imaginary->data.d[x + y*real->nx], out->data.cmpl_f[x + y*real->nx]);
                }
            }
            break;
        default:
            fprintf(stderr, "Unfortunately im_cmpl_f_make only works with double inputs at present\n");
            status->status = INVALID_TYPE;
            break;
    }
    CHECK_STATUS(status);
    return(status);
}
