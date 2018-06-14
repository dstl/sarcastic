/***************************************************************************
 *
 *       Module:    dataio_envi.c
 *      Program:    SIlib2
 *   Created by:    Matt Nottingham on 21/09/2007.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains the prototypes for the functions (def'ed in dataio_envi.c)
 *      for reading and writing ENVI files
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

static void load_envi_meta(SPImage *a, FILE * fp, int * do_swap, off_t * offset, SPStatus * status);

// Loads up a complete ENVI image from file
//
SPStatus* im_load_envi (SPImage *a, const char *fname, SPStatus *status) 
{
    FILE *fp_data;
    FILE *fp_hdr;
    int64_t f;
    char * hdr_name;
    int do_swap = 0;
    off_t offset;
    
    CHECK_STATUS(status);
    CHECK_PTR(a,status);
    CHECK_PTR(fname,status);
    
    hdr_name = calloc(strlen(fname)+10, sizeof(char));
    CHECK_PTR(hdr_name,status);
    
    snprintf(hdr_name, strlen(fname)+10, "%s.hdr", fname);
    fp_hdr = fopen(hdr_name, "r");
    if (fp_hdr == NULL) {
        snprintf(hdr_name, strlen(fname)+10, "%s.HDR", fname);
        fp_hdr = fopen(hdr_name, "r");
        if (fp_hdr == NULL) {
            CHECK_FP(fp_hdr,status);
            CHECK_STATUS(status);
        }
    }
    
    im_init(a, status);
    
    fp_data = fopen(fname, "r");
    CHECK_FP(fp_data,status);
    CHECK_STATUS(status);
    
    load_envi_meta(a, fp_hdr, &do_swap, &offset, status);
    fclose(fp_hdr);
    CHECK_STATUS(status);
    
    im_create (a, a->image_type, a->nx, a->ny, a->xspc, a->yspc, status);
    
    CHECK_IMCONT(a,status);
    
    fseeko(fp_data, offset, SEEK_SET);
    
    f = fread(a->data.v, im_getsizeoftype(a->image_type), a->nx*a->ny, fp_data);
    CHECK_BYTES(a->nx*a->ny, f, status);
    
    swap_data(a, do_swap);
    
    CHECK_IMINTEG(a, status);
    
    fclose(fp_data);
    
    return(status);
}

static void load_envi_meta(SPImage *a, FILE * fp, int * do_swap, off_t * offset, SPStatus * status)
{
    char line[256];
    int tmp;
    
    while (!feof(fp))
    {
        fgets(line, 256, fp);
        if (!feof(fp))
        {
            if (strncmp("samples", line, 7) == 0) {
                sscanf(&line[10], "%d", &tmp);
                a->nx = tmp;
                printf("nx tmp = %d\n", tmp);
            }
            if (strncmp("lines", line, 5) == 0) {
                sscanf(&line[10], "%d", &tmp);
                a->ny = tmp;
                printf("ny tmp = %d\n", tmp);
            }
            if (strncmp("bands", line, 5) == 0) {
                sscanf(&line[10], "%d", &tmp);
                printf("bands tmp = %d\n", tmp);
                
                if (tmp != 1) {
                    status->status = INVALID_TYPE;
                    fprintf(stderr, "Too many bands (%d)\n", tmp);
                    return;
                }
            }
            if (strncmp("data type", line, 9) == 0) {
                sscanf(&line[12], "%d", &tmp);
                printf("dt tmp = %d\n", tmp);
                
                a->image_type = im_conv_from_envi_type(tmp);
                
                
                if (a->image_type == ITYPE_UNKNOWN) {
                    fprintf(stderr, "Unknown datatype in ENVI hdr (%d)\n", tmp);
                    status->status = INVALID_TYPE;
                    return;
                }
            }
            
            if (strncmp("byte order", line, 10) == 0) {
                sscanf(&line[13], "%d", &tmp);
                printf("bo tmp = %d\n", tmp);
                
                if ((im_machine_type() == IM_LITTLE_ENDIAN && tmp == 0) ||
                    (im_machine_type() == IM_BIG_ENDIAN && tmp == 1))
                {
                    *do_swap = 0;
                }
                else
                {
                    *do_swap = 1;
                }
            }
            if (strncmp("header offset", line, 13) == 0) {
                sscanf(&line[16], "%d", &tmp);
                printf("oset tmp = %d\n", tmp);
                
                *offset = tmp;
            }
            
            if (strncmp("file type", line, 9) == 0) {
                if (strncmp("ENVI Standard", &line[12], 13) != 0) {
                    status->status = BAD_FILE;
                    fprintf(stderr, "Bad file type in ENVI hdr (%s)\n", &line[12]);
                    return;
                }
            }
        }
    }
}

// Loads up just the metadata from file
//
SPStatus* im_load_envi_metadata (SPImage *a, const char *fname, SPStatus *status)  
{
    FILE *fp;
    int do_swap = 0;
    off_t offset;
    
    CHECK_STATUS(status);
    CHECK_PTR(a,status);
    im_init(a, status);
    
    fp = fopen(fname, "r");
    CHECK_FP(fp,status);
    CHECK_STATUS(status);
    
    load_envi_meta(a, fp, &do_swap, &offset, status);
    fclose(fp);
    
    CHECK_STATUS(status);
    
    return(status);
}

// Loads up a subset of an image from a file
//
SPStatus * im_load_envi_subset (SPImage *a, const char *fname, int64_t ox, int64_t oy, SPStatus *status) 
{
    FILE *fp;
    int64_t f;
    int64_t y;
    int do_swap = 0;
    SPImage orig;
    off_t offset;
    
    CHECK_STATUS(status);
    CHECK_PTR(a,status);
    
    fp = fopen(fname, "r");
    CHECK_FP(fp,status);
    CHECK_STATUS(status);
    
    load_envi_meta(&orig, fp, &do_swap, &offset, status);
    
    CHECK_IMCONT(a,status);
    if (a->image_type != orig.image_type) {
        fclose(fp);
        fprintf(stderr, "Original image and subimage are different types!");
        status->status = TYPES_NOT_SAME;
        return status;
    }
    
    fseeko(fp, (ox + oy * orig.nx) * im_getsizeoftype(a->image_type) + offset, SEEK_SET);
    
    for(y = 0; y < a->ny; y++) {
        f = fread(&a->data.i8[im_getsizeoftype(a->image_type)*a->nx*y], im_getsizeoftype(a->image_type), a->nx, fp);
        CHECK_BYTES(a->nx, f, status);
        fseeko(fp, (orig.nx - (ox + a->nx) + ox) * im_getsizeoftype(a->image_type), SEEK_CUR);
    }
    
    swap_data(a, do_swap);
        
    CHECK_IMINTEG(a, status);
    
    fclose(fp);
    
    return(status);
}

SPImageType
im_conv_from_envi_type(EnviImageType itype)
{
    SPImageType image_type= ITYPE_UNKNOWN;
    
    switch(itype) {
        case ENVI_ITYPE_UINT8:
            image_type = ITYPE_UINT8;
            break;
            
        case ENVI_ITYPE_INT16:
            image_type = ITYPE_INT16;
            break;
            
        case ENVI_ITYPE_INT32:
            image_type = ITYPE_INT32;
            break;
            
        case ENVI_ITYPE_FLOAT:
            image_type = ITYPE_FLOAT;
            break;
            
        case ENVI_ITYPE_DOUBLE:
            image_type = ITYPE_DOUBLE;
            break;
            
        case ENVI_ITYPE_CMPL_FLOAT:
            image_type = ITYPE_CMPL_FLOAT;
            break;
            
        case ENVI_ITYPE_CMPL_DOUBLE:
            image_type = ITYPE_UNKNOWN;
            break;
            
        case ENVI_ITYPE_UINT16:
            image_type = ITYPE_UINT16;
            break;
            
        case ENVI_ITYPE_UINT32:
            image_type = ITYPE_UINT32;
            break;
            
        case ENVI_ITYPE_INT64:
            image_type = ITYPE_INT64;
            break;
            
        case ENVI_ITYPE_UINT64:
            image_type = ITYPE_UINT64;
            break;
            
        default:
            image_type = ITYPE_UNKNOWN;
            break;
    }
    
    return (image_type);
}

EnviImageType
im_conv_to_envi_type(SPImageType itype)
{
    EnviImageType image_type = ENVI_ITYPE_UNKNOWN;
    
    switch (itype) {
        case ITYPE_POLAR:
            image_type = ENVI_ITYPE_UNKNOWN;
            break;
            
        case ITYPE_CMPL_FLOAT:
            image_type = ENVI_ITYPE_CMPL_FLOAT;
            break;
            
        case ITYPE_FLOAT:
            image_type = ENVI_ITYPE_FLOAT;
            break;
            
        case ITYPE_DOUBLE:
            image_type = ENVI_ITYPE_DOUBLE;
            break;
            
        case ITYPE_INT64:
            image_type = ENVI_ITYPE_INT64;
            break;
            
        case ITYPE_INT32:
            image_type = ENVI_ITYPE_INT32;
            break;
            
        case ITYPE_INT16:
            image_type = ENVI_ITYPE_INT16;
            break;
            
        case ITYPE_INT8:
            image_type = ENVI_ITYPE_UNKNOWN;
            break;
            
        case ITYPE_UINT64:
            image_type = ENVI_ITYPE_UINT64;
            break;
            
        case ITYPE_UINT32:
            image_type = ENVI_ITYPE_UINT32;
            break;
            
        case ITYPE_UINT16:
            image_type = ENVI_ITYPE_UINT16;
            break;
            
        case ITYPE_UINT8:    
            image_type = ENVI_ITYPE_UINT8;
            break;
            
        case ITYPE_VECTOR:
            image_type = ENVI_ITYPE_UNKNOWN;
            break;
            
        case ITYPE_CMPL_INT64:
            image_type = ENVI_ITYPE_UNKNOWN;
            break;
            
        case ITYPE_CMPL_INT32:
            image_type = ENVI_ITYPE_UNKNOWN;
            break;
            
        case ITYPE_CMPL_INT16:
            image_type = ENVI_ITYPE_UNKNOWN;
            break;
            
        case ITYPE_CMPL_INT8:
            image_type = ENVI_ITYPE_UNKNOWN;
            break;
            
        default:
            image_type = ENVI_ITYPE_UNKNOWN;
            break;
    }
    
    return (image_type);
}

