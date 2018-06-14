/***************************************************************************
 *
 *       Module:    dataio.h
 *      Program:    SIlib2
 *   Created by:    Emma Griffiths on 01/02/2005.
 *                  and Darren Muff 7th Jun 2017
 *                  Copyright (c) 2017 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains functions for reading and writing to and from files to
 *      the variable types used in the time-domain SAR processor
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  7/06/2017
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

size_t
fread_byte_swap(void * ptr, size_t size, size_t nitem, FILE * stream, int do_swap, SPStatus *status)
{
    size_t f;
    
    CHECK_STATUS(status);
    
    f = fread(ptr, size, nitem, stream);
    CHECK_BYTES(f, nitem, status);
    
    if (do_swap)
    {
        im_swap_bytes(ptr, size, nitem);
    }
    
    return(f);
}

size_t
fwrite_byte_swap(void * ptr, size_t size, size_t nitem, FILE * stream, int do_swap, SPStatus *status)
{
    size_t f;
    
    CHECK_STATUS(status);
    
    unsigned char * buff = (unsigned char *)malloc(sizeof(char) * size * nitem);
    memcpy(buff, ptr, size * nitem);
    
    if (do_swap) {
        im_swap_bytes(buff, size, nitem);
    }
    
    f = fwrite(buff, size, nitem, stream);
    CHECK_BYTES(f, nitem, status);
    
    free(buff) ;
    
    return(f);
}

void
im_swap_bytes(void * ptr, size_t size, size_t nitem)
{
    unsigned char * p;
    size_t n;
    size_t b;
    unsigned char tmp;
    p = (unsigned char *)ptr;
    
    for(n = 0; n < nitem; n++)
    {
        for(b = 0; b < size/2; b++)
        {
            tmp = *(p+b);
            *(p+b) = *(p+size-b-1);
            *(p+size-b-1) = tmp;
        }
        p+=size;
    }
}
