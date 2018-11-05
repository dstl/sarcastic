/***************************************************************************
 * 
 *           Module :  dataio_hdf4.c
 *          Program :  sarclib
 *       Created by :  Matt Nottingham on 02/07/2009
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      This file contains functions for reading HDF4 byte files.
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

#include "dataio.h"
#include "error_defn.h"
#include "dataio_hdf4.h"

static SPStatus * im_read_header_hdf4(HDF4header * hdr, SPStatus * status);

SPStatus *
im_open_hdf4(HDF4header * hdr, const char * filename, SPStatus * status)
{
    int d[4];
    int i;
    
    CHECK_STATUS(status);
    CHECK_PTR(hdr, status);
    
    hdr->nx = -1;
    hdr->ny = -1;
    
    hdr->num_images = 0;
    hdr->imageOffset = NULL;
    
    hdr->fp = fopen(filename, "r");
    CHECK_FP(hdr->fp, status);
    
    for (i = 0; i < 4; i++) {
        d[i] = fgetc(hdr->fp);
    }
    
    if (d[0] != 0x0E || d[1] != 0x03 || d[2] != 0x13 || d[3] != 0x01) {
        fprintf(stderr, "File \"%s\" does not appear to be a valid HDF4 file", filename);
        status->status = BAD_FILE;
        return (status);
    }
    
    hdr->imageOffset = calloc(256, sizeof(int));
    
    CHECK_PTR(hdr->imageOffset, status);
    
    if (im_machine_type() == IM_LITTLE_ENDIAN) {
        hdr->do_swap = 1;
    } else {
        hdr->do_swap = 0;
    }
    
    im_read_header_hdf4(hdr, status);
    
    return(status);
}

static SPStatus *
im_read_header_hdf4(HDF4header * hdr, SPStatus * status)
{
    uint16_t ntags;
    uint32_t next_block;
    uint16_t tagnum;
    uint16_t refnum;
    uint32_t dataposition;
    uint32_t datasize;
    int i;
    off_t oldpos;
    
    do {
        fread_byte_swap(&ntags, sizeof(uint16_t), 1, hdr->fp, hdr->do_swap, status);
        fread_byte_swap(&next_block, sizeof(uint32_t), 1, hdr->fp, hdr->do_swap, status);
        
        for(i = 0; i < ntags; i++) {
            fread_byte_swap(&tagnum, sizeof(uint16_t), 1, hdr->fp, hdr->do_swap, status);
            fread_byte_swap(&refnum, sizeof(uint16_t), 1, hdr->fp, hdr->do_swap, status);
            fread_byte_swap(&dataposition, sizeof(uint32_t), 1, hdr->fp, hdr->do_swap, status);
            fread_byte_swap(&datasize, sizeof(uint32_t), 1, hdr->fp, hdr->do_swap, status);
            
            if (tagnum == 302) {
                hdr->imageOffset[hdr->num_images] = dataposition;
                hdr->num_images++;
            }
            
            if (tagnum == 300 && hdr->nx == -1 && hdr->ny == -1) {
                oldpos = ftello(hdr->fp);
                fseeko(hdr->fp, dataposition, SEEK_SET);
                fread_byte_swap(&hdr->nx, sizeof(uint32_t), 1, hdr->fp, hdr->do_swap, status);
                fread_byte_swap(&hdr->ny, sizeof(uint32_t), 1, hdr->fp, hdr->do_swap, status);
                fseeko(hdr->fp, oldpos, SEEK_SET);
            }
            
            if (tagnum == 7261 || tagnum == 101) {
                // this is metadata
            }
            
            if (tagnum == 14) {
                // This is version info
            }
            
        }
        
        if (next_block > 0) {
            fseeko(hdr->fp, next_block, SEEK_SET);
        }
    } while (next_block > 0);
    
    return (status);
}

// Read in a complete HDF4 file
//
SPStatus *
im_load_hdf4(SPImage * im, HDF4header * hdr, int frame, SPStatus * status) {
    
    CHECK_STATUS(status);
    CHECK_PTR(hdr, status);
    CHECK_PTR(im, status);
    
    if (frame > hdr->num_images || frame < 0) {
        fprintf(stderr, "im_load_hdf4: Frame requested (%d) is not between 0 and %d inclusive\n",
                frame, hdr->num_images - 1);
        
        status->status = OUT_OF_BOUND;
        
        return(status);
    }
    
    im_create(im, ITYPE_UINT8, hdr->nx, hdr->ny, 1.0, 1.0, status);
    CHECK_STATUS(status);
    
    im_load_hdf4_subset(im, hdr, frame, 0, 0, status);
    
    CHECK_STATUS(status);
    
    return(status);
}

// reads in a im->nx, im->ny subset of image from offset ox, oy
//
SPStatus*
im_load_hdf4_subset(SPImage *im, HDF4header * hdr, int frame, int ox, int oy, SPStatus * status) 
{
    long y;
    
    CHECK_STATUS(status);
    CHECK_PTR(hdr, status);
    
    if (frame > hdr->num_images || frame < 0) {
        fprintf(stderr, "im_load_hdf4: Frame requested (%d) is not between 0 and %d inclusive\n",
                frame, hdr->num_images-1);
        
        status->status = OUT_OF_BOUND;
        
        return(status);
    }
    
    CHECK_PTR(im, status);
    
    CHECK_STATUS(status);
    
    if (im->nx <= 0 || ox < 0 || (im->nx + ox) > hdr->nx)
    {
        fprintf(stderr, "im_load_hdf4_subset: Input image nx (%ld) + ox (%d) is either < 0 or > %d\n",
                (long) im->nx, ox, hdr->nx);
        status->status = INPUT_NX_MISMATCHED;
        return(status);
    }
    
    if (im->ny <= 0 || oy < 0 || (im->ny + oy) > hdr->ny)
    {
        fprintf(stderr, "im_load_hdf4_subset: Input image ny (%ld) + oy (%d) is either < 0 or > %d\n",
                (long) im->ny, oy, hdr->ny);
        status->status = INPUT_NY_MISMATCHED;
        return(status);
    }
    
    fseeko(hdr->fp, (oy * hdr->nx) * sizeof(uint8_t) + hdr->imageOffset[frame], SEEK_SET);
    
    for(y = 0; y < im->ny; y++)
    {
        fseeko(hdr->fp, ox * sizeof(uint8_t), SEEK_CUR);
        fread(&im->data.ui8[y * im->nx], sizeof(uint8_t), im->nx, hdr->fp);
        
        fseeko(hdr->fp, (hdr->nx - ox - im->nx) * sizeof(uint8_t), SEEK_CUR);
    }
    
    return(status);
}


SPStatus *
im_close_hdf4(HDF4header * hdr, SPStatus * status)
{
    CHECK_STATUS(status);
    CHECK_PTR(hdr, status);
    CHECK_FP(hdr->fp, status);
    
    fclose(hdr->fp);
    
    hdr->fp = NULL;
    
    return status;
}

SPStatus *
im_destroy_hdf4(HDF4header * hdr, SPStatus * status)
{
    CHECK_STATUS(status);
    CHECK_PTR(hdr, status);
    
    if (hdr->fp) {
        im_close_hdf4(hdr, status);
        CHECK_STATUS(status);
    }
    
    if (hdr->imageOffset) {
        free(hdr->imageOffset);
    }
    
    return status;
}
