/***************************************************************************
 *
 *       Module:    dataio_ce.c
 *      Program:    sarclib
 *   Created by:    Matt Nottingham on 06/01/2006.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains functions for reading from CaseExec files
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

#include "dataio_ce.h"

// Reads in the metadata from a case exec dataset
//
SPStatus*
im_info_CEdataset(SPCEinfo * ce_info, const char * dir, SPStatus * status) {
    char * temp;
    FILE * sd;
    FILE * ci;
    
    CHECK_STATUS(status);
    CHECK_PTR(ce_info, status);
    CHECK_PTR(dir, status);
    
    switch (im_machine_type())
    {
        case IM_BIG_ENDIAN:
            ce_info->do_swap = 0;
            break;
            
        case IM_LITTLE_ENDIAN:
            ce_info->do_swap = 1;
            break;
            
        default:
            fprintf(stderr, "Unknown endian in im_load_CEdataset!!!\n");
            status->status = ENDIAN_ERROR;
            return(status);
            
    }
    
    ce_info->sd_name = calloc(strlen(dir)+30, sizeof(char));
    CHECK_PTR(ce_info->sd_name, status);
    
    snprintf(ce_info->sd_name, strlen(dir)+30, "%s/_Dataset.sd", dir);
    sd = fopen(ce_info->sd_name, "r");
    CHECK_FP(sd, status);
    
    ce_info->ci_name = calloc(strlen(dir)+30, sizeof(char));
    CHECK_PTR(ce_info->ci_name, status);
    
    snprintf(ce_info->ci_name, strlen(dir)+30, "%s/_Dataset.ci", dir);
    ci = fopen(ce_info->sd_name, "r");
    CHECK_FP(ci, status);
    fclose(ci);
    
    temp = calloc(1024, sizeof(char));
    CHECK_PTR(temp, status);
    
    while (!feof(sd)) {
        fgets(temp, 1024, sd);
        if (strncmp("ImageWidth=", temp, 11) == 0) {
            ce_info->nx = atol(&temp[11]);          // Width in pixels
        }
        if (strncmp("ImageHeight=", temp, 12) == 0) {
            ce_info->ny = atol(&temp[12]);          // Height in pixels
        }
        if (strncmp("SD_CR=", temp, 6) == 0) {
            ce_info->sx = 0.3048 / atof(&temp[6]);  // Convert into metres
        }
        if (strncmp("SD_RG=", temp, 6) == 0) {
            ce_info->sy = 0.3048 / atof(&temp[6]);  // Convert into metres
        }
    }
    free(temp);
    fclose(sd);
    
    if (ce_info->nx == 0 || ce_info->ny == 0) {
        fprintf(stderr, "Problem reading image size from %s/_Dataset.sd\n", dir);
        status->status = BAD_IMAGE;
        return (status);
    }
    
    return(status);
}

// Read in a complete case exec dataset from the given directory
//
SPStatus*
im_load_CEdataset(SPImage * im, const char * dir, SPStatus * status) {
    SPCEinfo info;
    CHECK_STATUS(status);
    im_info_CEdataset(&info, dir, status);
    
    CHECK_STATUS(status);
    im_create(im, ITYPE_CMPL_FLOAT, info.nx, info.ny, info.sx, info.sy, status);
    CHECK_STATUS(status);
    
    im_load_CEdataset_subset(im, &info, 0, 0, status);
    
    CHECK_STATUS(status);
    
    im_destroy_SPCEinfo(&info, status);
    
    return(status);
}

// reads in a im->nx, im->ny subset of image from offset ox, oy
// using the already filled in ce_info structure
//
SPStatus*
im_load_CEdataset_subset(SPImage *im, SPCEinfo * ce_info, long ox, long oy, SPStatus * status) {
    FILE * ci;
    long x, y;
    size_t f;
    SPCmplxInt16 * buffer;
    
    CHECK_STATUS(status);
    
    CHECK_PTR(ce_info, status);
    
    CHECK_PTR(im, status);
    
    if (im->nx <= 0 || ox < 0 || (im->nx + ox) > ce_info->nx)
    {
        fprintf(stderr, "im_load_CEdataset_subset: Input image nx (%ld) + ox (%ld) is either < 0 or > %ld\n",
                (long) im->nx, ox, ce_info->nx);
        status->status = INPUT_NX_MISMATCHED;
        return(status);
    }
    
    if (im->ny <= 0 || oy < 0 || (im->ny + oy) > ce_info->ny)
    {
        fprintf(stderr, "im_load_CEdataset_subset: Input image ny (%ld) + oy (%ld) is either < 0 or > %ld\n",
                (long) im->ny, oy, ce_info->ny);
        status->status = INPUT_NY_MISMATCHED;
        return(status);
    }
    
    buffer = calloc(im->nx, sizeof(SPCmplxInt16));
    CHECK_PTR(buffer, status);
    
    ci = fopen(ce_info->ci_name, "r");
    CHECK_FP(ci, status);
    
    fseek(ci, (oy * ce_info->nx) * sizeof(SPCmplxInt16), SEEK_SET);
    
    for(y = 0; y < im->ny; y++)
    {
        fseek(ci, ox * sizeof(SPCmplxInt16), SEEK_CUR);
        f = fread_byte_swap(buffer, sizeof(int16_t), 2 * im->nx, ci, ce_info->do_swap, status);
        
        CHECK_BYTES(f, 2 * im->nx, status);
        for(x = 0; x < im->nx; x++)
        {
            im->data.cmpl_f[x + y * im->nx].r = buffer[x].r;
            im->data.cmpl_f[x + y * im->nx].i = buffer[x].i;
        }
        
        fseek(ci, (ce_info->nx - ox - im->nx) * sizeof(SPCmplxInt16), SEEK_CUR);
    }
    free(buffer);
    fclose(ci);
    
    return(status);
}

// Clean up the SPCEinfo structure
//
SPStatus*
im_destroy_SPCEinfo(SPCEinfo * info, SPStatus * status)
{
    CHECK_STATUS(status);
    if (info->ci_name) free(info->ci_name);
    if (info->sd_name) free(info->sd_name);
    
    return(status);
}

