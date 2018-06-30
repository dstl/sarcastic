/***************************************************************************
 *
 *       Module:    dataio_im.c
 *      Program:    sarclib
 *   Created by:    Emma Griffiths on 19/07/2006.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains functions for reading and writing to and from files.
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

#include "dataio.h"
#include "error_defn.h"
#include "dataio_im.h"

#define N 256
#define N2 200

#define IMAGE_MAGIC 0x87654321

static SPStatus * im_save_envi_hdr(SPImage * a, const char * fname, SPStatus * status);
static SPStatus * im_save_envi_hdr_with_geo_ll(SPImage * a, SPVector tl_ll, SPVector tr_ll, SPVector br_ll, SPVector bl_ll, const char * fname, SPStatus * status);
static SPStatus* load_meta(SPImage * a, FILE * fp, int * do_swap, SPStatus * status);

// Saves the header for an image
//
// The header conists of:
//
//  - a uint32  magic number 0x87654321
//  - a uint64 pair of numbers describing the size of the image (x then y axis)
//  - a pair of doubles describing the pixels spacing of each axis (x then y axis)
//  - an int32 describing the image type (eg. complex float = 2, double = 4,
//
//  look in image.h for the declaration of SPImageType for the values of all the
//     other types
//
SPStatus* im_save_header(SPImage * a, FILE * fp, SPStatus * status)
{
    int magic = IMAGE_MAGIC;
    int64_t f;
    
    CHECK_STATUS(status);
    CHECK_PTR(a, status);
    CHECK_PTR(fp, status);
    
    f = fwrite(&magic, sizeof(magic), 1, fp);
    CHECK_BYTES(1, f, status);
    f = fwrite(&a->nx, sizeof(a->nx), 1, fp);
    CHECK_BYTES(1, f, status);
    f = fwrite(&a->ny, sizeof(a->ny), 1, fp);
    CHECK_BYTES(1, f, status);
    f = fwrite(&a->xspc, sizeof(a->xspc), 1, fp);
    CHECK_BYTES(1, f, status);
    f = fwrite(&a->yspc, sizeof(a->yspc), 1, fp);
    CHECK_BYTES(1, f, status);
    f = fwrite(&a->image_type, sizeof(a->image_type), 1, fp);
    CHECK_BYTES(1, f, status);
    
    return(status);
}

// Save a complete image to a file
//
SPStatus * im_save (SPImage *a, const char *fname, SPStatus *status)
{
    FILE *fp;
    int64_t f;
    
    CHECK_STATUS(status);
    CHECK_PTR(a,status);
    CHECK_IMCONT(a,status);
    CHECK_IMINTEG(a,status);
    
    fp = fopen(fname, "w");
    
    CHECK_FP(fp,status);
    
    CHECK_STATUS(status);
    
    im_save_header(a, fp, status);
    CHECK_STATUS(status);
    
    f = fwrite(a->data.v,  im_getsizeoftype(a->image_type), a->nx*a->ny, fp);
    CHECK_BYTES(a->nx*a->ny, f, status);
    
    fclose(fp);
    
    im_save_envi_hdr(a, fname, status);
    
    return(status);
}

// Save a complete image to a file but also specifying the
// corner coordinates of the image. Corner coordinates are specified
// as ECEF vectors of teh corners of the image starting with the top left
// and going clockwise around the image (tl,tr,br,bl)
//
SPStatus * im_save_with_geo_ecf(SPImage *a, SPVector tl, SPVector tr, SPVector br, SPVector bl, const char *fname, SPStatus *status)
{
    FILE *fp;
    int64_t f;
    SPVector tl_ll;
    SPVector tr_ll;
    SPVector br_ll;
    SPVector bl_ll;
    
    CHECK_STATUS(status);
    CHECK_PTR(a,status);
    CHECK_IMCONT(a,status);
    CHECK_IMINTEG(a,status);
    
    fp = fopen(fname, "w");
    
    CHECK_FP(fp,status);
    
    CHECK_STATUS(status);
    
    im_save_header(a, fp, status);
    CHECK_STATUS(status);
    
    f = fwrite(a->data.v,  im_getsizeoftype(a->image_type), a->nx*a->ny, fp);
    CHECK_BYTES(a->nx*a->ny, f, status);
    
    fclose(fp);
    
    ecef2latlon(&tl, &tl_ll, status);
    ecef2latlon(&tr, &tr_ll, status);
    ecef2latlon(&br, &br_ll, status);
    ecef2latlon(&bl, &bl_ll, status);
    
    im_save_envi_hdr_with_geo_ll(a, tl_ll, tr_ll, br_ll, bl_ll, fname, status);
    
    return(status);
}

// Initialises the line save info structure so that a file can be
// written a few lines at a time
//
SPStatus* im_save_line_init(SPImage * im, const char * fname, SPImageLineSaveInfo * info, SPStatus * status)
{
    int64_t temp;
    
    CHECK_STATUS(status);
    CHECK_PTR(im,status);
    CHECK_IMCONT(im,status);
    CHECK_IMINTEG(im,status);
    CHECK_PTR(fname,status);
    CHECK_PTR(info,status);
    
    info->fname = strdup(fname);
    info->fp = fopen(fname, "w+");
    CHECK_FP(info->fp,status);
    
    info->nx = im->nx;
    temp = im->ny;
    im->ny = 0;
    im_save_header(im, info->fp, status);
    im->ny = temp;
    
    CHECK_STATUS(status);
    
    return(status);
}

// Saves an image to a file that is already partially written. This is useful for writing
// an image to a file a piece at a time. The output file must be already initialised with
// im_save_line_init().
//
SPStatus* im_save_line(SPImage * im, SPImageLineSaveInfo * info, SPStatus * status)
{
    int64_t ny;
    int64_t f;
    
    CHECK_STATUS(status);
    CHECK_PTR(im,status);
    CHECK_IMCONT(im,status);
    CHECK_IMINTEG(im,status);
    CHECK_PTR(info,status);
    
    if (im->nx != info->nx)
    {
        fprintf(stderr, "im_save_line: Trying to write out line of different width! (%ld != %ld)\n", (long)im->nx, (long)info->nx);
        status->status = INPUT_NX_MISMATCHED;
        return(status);
    }
    
    fseek(info->fp, sizeof(int)+sizeof(im->nx), SEEK_SET);
    f = fread(&ny, sizeof(ny), 1, info->fp);
    CHECK_BYTES(1, f, status);
    
    ny += im->ny;
    
    fseek(info->fp, sizeof(int)+sizeof(im->nx), SEEK_SET);
    f = fwrite(&ny, sizeof(ny), 1, info->fp);
    CHECK_BYTES(1, f, status);
    
    fseek(info->fp, 0, SEEK_END);
    f = fwrite(im->data.v, im_getsizeoftype(im->image_type), im->nx * im->ny, info->fp);
    CHECK_BYTES(im->nx * im->ny, f, status);
    
    return(status);
}

// Closes down the line save info structure and file
//
SPStatus* im_save_line_close_with_geo_ecf(SPImageLineSaveInfo * info, SPVector p1, SPVector p2, SPVector p3, SPVector p4, SPStatus * status)
{
    SPVector tl_ll;
    SPVector tr_ll;
    SPVector br_ll;
    SPVector bl_ll;
    SPImage tmp;
    CHECK_STATUS(status);
    CHECK_PTR(info,status);
    CHECK_PTR(info->fp, status);
    
    ecef2latlon(&p1, &tl_ll, status);
    ecef2latlon(&p2, &tr_ll, status);
    ecef2latlon(&p3, &br_ll, status);
    ecef2latlon(&p4, &bl_ll, status);
    
    im_load_metadata(&tmp, info->fname, status);
    im_save_envi_hdr_with_geo_ll(&tmp, tl_ll, tr_ll, br_ll, bl_ll, info->fname, status);
    
    free(info->fname);
    fclose(info->fp);
    return(status);
}

// Closes down the line save info structure and file
//
SPStatus* im_save_line_close(SPImageLineSaveInfo * info, SPStatus * status)
{
    SPImage tmp;
    CHECK_STATUS(status);
    CHECK_PTR(info,status);
    CHECK_PTR(info->fp, status);
    
    im_load_metadata(&tmp, info->fname, status);
    im_save_envi_hdr(&tmp, info->fname, status);
    
    free(info->fname);
    fclose(info->fp);
    return(status);
}

// This saves out a text file format description file that ENVI can read
//
static
SPStatus * im_save_envi_hdr(SPImage * a, const char * fname, SPStatus * status)
{
    SPVector tl_ll;
    SPVector tr_ll;
    SPVector br_ll;
    SPVector bl_ll;
    
    tl_ll.x = 0.0; tl_ll.y = 0.0; tl_ll.z = 0.0;
    tr_ll.x = 0.0; tr_ll.y = 0.0; tr_ll.z = 0.0;
    br_ll.x = 0.0; br_ll.y = 0.0; br_ll.z = 0.0;
    bl_ll.x = 0.0; bl_ll.y = 0.0; bl_ll.z = 0.0;
    
    im_save_envi_hdr_with_geo_ll(a, tl_ll, tr_ll, br_ll, bl_ll, fname, status);
    
    return(status);
}

static
SPStatus * im_save_envi_hdr_with_geo_ll(SPImage * a, SPVector tl_ll, SPVector tr_ll, SPVector br_ll, SPVector bl_ll, const char * fname, SPStatus * status)
{
    FILE * fp;
    char * hdr_name;
    EnviImageType type;
    
    CHECK_STATUS(status);
    
    hdr_name = malloc(strlen(fname)+10);
    CHECK_PTR(hdr_name, status);
    
    sprintf(hdr_name, "%s.hdr", fname);
    
    fp = fopen(hdr_name, "w");
    
    CHECK_FP(fp, status);
    
    fprintf(fp, "ENVI\n");
    fprintf(fp, "description = {\n  File Imported into ENVI from SPImage.}\n");
    fprintf(fp, "samples = %ld\n", (long) a->nx);
    fprintf(fp, "lines   = %ld\n", (long) a->ny);
    fprintf(fp, "header offset = 40\n");
    fprintf(fp, "file type = ENVI Standard\n");
    
    if (a->image_type != ITYPE_VECTOR)
    {
        fprintf(fp, "bands   = 1\n");
        type = im_conv_to_envi_type(a->image_type);
        
        fprintf(fp, "data type = %d\n", type);
        fprintf(fp, "interleave = bsq\n");
        fprintf(fp, "sensor type = Unknown\n");
        if (im_machine_type() == IM_LITTLE_ENDIAN)
        {
            fprintf(fp, "byte order = 0\n");
        }
        else
        {
            fprintf(fp, "byte order = 1\n");
        }
        fprintf(fp, "wavelength units = Unknown\n");
        
        if (!(tl_ll.x == 0.0 && tl_ll.y == 0.0 && tl_ll.z == 0.0 &&
              tr_ll.x == 0.0 && tr_ll.y == 0.0 && tr_ll.z == 0.0 &&
              br_ll.x == 0.0 && br_ll.y == 0.0 && br_ll.z == 0.0 &&
              bl_ll.x == 0.0 && bl_ll.y == 0.0 && bl_ll.z == 0.0)) {
            fprintf(fp, "geo points = {\n");
            fprintf(fp, " 1.5000, 1.500, %15.12g, %15.12g,\n",                tl_ll.x,                    tl_ll.y);
            fprintf(fp, " %14.10g, 1.500, %15.12g, %15.12g,\n",   (double)a->nx+0.5,             tr_ll.x, tr_ll.y);
            fprintf(fp, " %14.10g, %14.10g, %15.12g, %15.12g,\n", (double)a->nx+0.5, (double)a->ny+0.5, br_ll.x, br_ll.y);
            fprintf(fp, " 1.5000, %14.10g, %15.12g, %15.12g}\n",  (double)a->ny+0.5,             bl_ll.x, bl_ll.y);
        }
        
        if (a->image_type == ITYPE_CMPL_FLOAT) {
            fprintf(fp, "complex function = Magnitude\n");
        }
        
    }
    else
    {
        fprintf(fp, "bands   = 3\n");
        
        fprintf(fp, "data type = 5\n");
        fprintf(fp, "interleave = bip\n");
        fprintf(fp, "sensor type = Unknown\n");
        if (im_machine_type() == IM_LITTLE_ENDIAN)
        {
            fprintf(fp, "byte order = 0\n");
        }
        else
        {
            fprintf(fp, "byte order = 1\n");  /* ???? */
        }
        fprintf(fp, "wavelength units = Unknown\n");
        if (a->image_type == ITYPE_CMPL_FLOAT) {
            fprintf(fp, "complex function = Magnitude\n");
        }
    }
    
    fclose(fp);
    free(hdr_name);
    
    return(status);
}

static SPStatus* load_meta(SPImage * a, FILE * fp, int * do_swap, SPStatus * status)
{
    int64_t f;
    int magic;
    
    *do_swap = 0;
    f = fread_byte_swap(&magic, sizeof(magic), 1, fp, *do_swap, status);
    CHECK_BYTES(1, f, status);
    
    if (magic != IMAGE_MAGIC)
    {
        *do_swap = 1;
        rewind(fp);
        f = fread_byte_swap(&magic, sizeof(magic), 1, fp, *do_swap, status);
        CHECK_BYTES(1, f, status);
        
        if (magic != IMAGE_MAGIC)
        {
            fprintf(stderr, "Could not read in magic number from start of file!\n");
            status->status = BAD_FILE;
            return (status);
        }
    }
    
    f = fread_byte_swap(&a->nx, sizeof(a->nx), 1, fp, *do_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&a->ny, sizeof(a->ny), 1, fp, *do_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&a->xspc, sizeof(a->xspc), 1, fp, *do_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&a->yspc, sizeof(a->yspc), 1, fp, *do_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&a->image_type, sizeof(a->image_type), 1, fp, *do_swap, status);
    CHECK_BYTES(1, f, status);
    
    return(status);
}

// Loads up just the metadata from file
//
SPStatus* im_load_metadata (SPImage *a, const char *fname, SPStatus *status)
{
    FILE *fp;
    int do_swap = 0;
    CHECK_STATUS(status);
    CHECK_PTR(a,status);
    im_init(a, status);
    
    fp = fopen(fname, "r");
    CHECK_FP(fp,status);
    CHECK_STATUS(status);
    
    load_meta(a, fp, &do_swap, status);
    
    fclose(fp);
    
    return(status);
}

// Loads up a complete image from file
//
SPStatus* im_load (SPImage *a, const char *fname, SPStatus *status)
{
    FILE *fp;
    int64_t f;
    int do_swap = 0;
    CHECK_STATUS(status);
    CHECK_PTR(a,status);
    
    im_init(a, status);
    
    fp = fopen(fname, "r");
    
    CHECK_FP(fp,status);
    CHECK_STATUS(status);
    
    load_meta(a, fp, &do_swap, status);
    
    im_create (a, a->image_type, a->nx, a->ny, a->xspc, a->yspc, status);
    
    CHECK_IMCONT(a,status);
    
    f = fread(a->data.v, im_getsizeoftype(a->image_type), a->nx*a->ny, fp);
    CHECK_BYTES(a->nx*a->ny, f, status);
    
    swap_data(a, do_swap);
    
    CHECK_IMINTEG(a, status);
    
    fclose(fp);
    
    return(status);
}

// Loads up a subset of an image from a file
//
SPStatus * im_load_subset (SPImage *a, const char *fname, int64_t ox, int64_t oy, SPStatus *status)
{
    FILE *fp;
    int64_t f;
    int64_t y;
    int do_swap = 0;
    SPImage orig;
    
    CHECK_STATUS(status);
    CHECK_PTR(a,status);
    
    fp = fopen(fname, "r");
    CHECK_FP(fp,status);
    CHECK_STATUS(status);
    
    load_meta(&orig, fp, &do_swap, status);
    
    CHECK_IMCONT(a,status);
    if (a->image_type != orig.image_type) {
        fclose(fp);
        fprintf(stderr, "Original image and subimage are different types!");
        status->status = TYPES_NOT_SAME;
        return status;
    }
    
    fseeko(fp, (ox + oy * orig.nx) * im_getsizeoftype(a->image_type), SEEK_CUR);
    
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

// Perform byte swapping on the data within an image
//
void swap_data(SPImage * a, int do_swap)
{
    int64_t nx;
    size_t sz = 0;
    
    if (do_swap)
    {
        nx = a->nx;
        switch(a->image_type)
        {
            case ITYPE_POLAR:
            case ITYPE_CMPL_FLOAT:
                nx *= 2;
                sz = sizeof(float);
                break;
                
            case ITYPE_DOUBLE:
            case ITYPE_INT64:
            case ITYPE_UINT64:
                sz = sizeof(double);
                break;
                
            case ITYPE_FLOAT:
            case ITYPE_INT32:
            case ITYPE_UINT32:
                sz = sizeof(float);
                break;
                
            case ITYPE_INT16:
            case ITYPE_UINT16:
                sz = sizeof(int16_t);
                break;
                
            case ITYPE_INT8:
            case ITYPE_UINT8:
            case ITYPE_CMPL_NIBBLE:
                sz = sizeof(int8_t);
                break;
                
            case ITYPE_VECTOR:
                nx *= 3;
                sz = sizeof(double);
                break;
                
            case ITYPE_CMPL_INT64:
                nx *= 2;
                sz = sizeof(int64_t);
                break;
                
            case ITYPE_CMPL_INT32:
                nx *= 2;
                sz = sizeof(int32_t);
                break;
                
            case ITYPE_CMPL_INT16:
                nx *= 2;
                sz = sizeof(int32_t);
                break;
                
            case ITYPE_CMPL_INT8:
                nx *= 2;
                sz = sizeof(int32_t);
                break;
                
            case ITYPE_UNKNOWN:
                nx = 0;
                sz = 0;
                break;
        }
        
        if (sz > 1) {
            im_swap_bytes(a->data.v, sz, a->ny * nx);
        }
    }
}


