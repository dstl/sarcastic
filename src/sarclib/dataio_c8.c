/***************************************************************************
 *
 *       Module:    dataio_c8.c
 *      Program:    sarclib
 *   Created by:    Matt Nottingham on 05/10/2006.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains functions for reading c8 files made by IFP4
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

#include "dataio_c8.h"

// Parse the c8 image filename
//
SPStatus*
im_c8_parse_fname(const char * fname, int * nx, int * ny, SPStatus *status)
{
    int p;
    char * fn = strdup(fname);
    /* The filename must be of the form xxxxxx.nx.ny.c8 */
    
    if (strncmp(&fname[strlen(fname)-3], ".c8", 3) != 0) {
        fprintf(stderr, "Badly formed c8 filename - %s\n", fname);
        status->status = BAD_FILE;
        return status;
    }
    fn[strlen(fname)-3] = 0;
    
    for(p = (int)strlen(fname)-4; p > 0 && fname[p] != '.'; p--); /* look for the next '.' */
    fn[p] = 0;
    
    sscanf(&fn[p+1], "%d", ny);
    p--;
    while(fname[p] != '.') {
        p--;
    }
    fn[p] = 0;
    sscanf(&fn[p+1], "%d", nx);
    
    return(status);
}

// Read in a complete c8 file from the given file name
//
SPStatus*
im_load_c8(SPImage * im, const char * filename, SPStatus * status) 
{
    
    int nx;
    int ny;
    
    CHECK_STATUS(status);
    CHECK_PTR(filename, status);
    
    im_c8_parse_fname(filename, &nx, &ny, status);
    
    im_create(im, ITYPE_CMPL_FLOAT, nx, ny, 1.0, 1.0, status);
    CHECK_STATUS(status);
    
    im_load_c8_subset(im, filename, 0, 0, status);
    
    CHECK_STATUS(status);
    
    return(status);
}

// reads in a im->nx, im->ny subset of image from offset ox, oy */
//
SPStatus*
im_load_c8_subset(SPImage *im, const char * fname, long ox, long oy, SPStatus * status)
{
    FILE * ci;
    long y;
    size_t f;
    int do_swap = 0;
    int nx;
    int ny;
    
    CHECK_STATUS(status);
    
    CHECK_PTR(fname, status);
    
    CHECK_PTR(im, status);
    
    im_c8_parse_fname(fname, &nx, &ny, status);
    CHECK_STATUS(status);
    
    if (im->nx <= 0 || ox < 0 || (im->nx + ox) > nx)
    {
        fprintf(stderr, "im_load_c8_subset: Input image nx (%ld) + ox (%ld) is either < 0 or > %d\n",
                (long) im->nx, ox, nx);
        status->status = INPUT_NX_MISMATCHED;
        return(status);
    }
    
    if (im->ny <= 0 || oy < 0 || (im->ny + oy) > ny)
    {
        fprintf(stderr, "im_load_c8_subset: Input image ny (%ld) + oy (%ld) is either < 0 or > %d\n",
                (long) im->ny, oy, ny);
        status->status = INPUT_NY_MISMATCHED;
        return(status);
    }
    
    
    ci = fopen(fname, "r");
    CHECK_FP(ci, status);
    
    fseeko(ci, (oy * nx) * sizeof(SPCmplx), SEEK_SET);
    
    for(y = 0; y < im->ny; y++)
    {
        fseeko(ci, ox * sizeof(SPCmplx), SEEK_CUR);
        f = fread_byte_swap(&im->data.cmpl_f[y * im->nx], sizeof(float), 2 * im->nx, ci, do_swap, status);
        
        fseeko(ci, (nx - ox - im->nx) * sizeof(SPCmplx), SEEK_CUR);
    }
    
    fclose(ci);
    
    return(status);
}

