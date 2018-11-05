/***************************************************************************
 * 
 *           Module :  dataio_sio.c
 *          Program :  sarclib
 *       Created by :  Matt Nottingham on 17/07/2006
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      Description: Files for reading/writing SIO format files.
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

SPStatus *
im_open_sio(SIOheader * hdr, const char * filename, SPStatus * status)
{
    size_t f;
    CHECK_STATUS(status);
    CHECK_PTR(hdr, status);
    printf("Filename %s\n", filename);
    hdr->fp = fopen(filename, "r");
    CHECK_FP(hdr->fp, status);
    
    hdr->do_byte_swap = 0;
    
    f = fread_byte_swap(&hdr->magic, sizeof(int), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    if (hdr->magic != SIO_MAGIC_BASIC && hdr->magic != SIO_MAGIC_EXTENDED)
    {
        fseeko(hdr->fp, -4, SEEK_CUR);   // Go back 4 bytes
        
        // Try reading it with byte swapping turned on
        //
        hdr->do_byte_swap = 1;
        
        f = fread_byte_swap(&hdr->magic, sizeof(int), 1, hdr->fp, hdr->do_byte_swap, status);
        CHECK_BYTES(1, f, status);
        
        if (hdr->magic != SIO_MAGIC_BASIC && hdr->magic != SIO_MAGIC_EXTENDED)
        {
            // We still don't have the magic number, so its not a valid SIO file.
            //
            fprintf(stderr, "Invalid SIO file - magic number = %x\n", hdr->magic);
            status->status = BAD_FILE;
            return (status);
        }
    }
    return status;
}

SPStatus *
im_read_header_sio(SIOheader * hdr, SPStatus * status)
{
    size_t f;
    int label;
    int label_length;
    int tmp;
    
    CHECK_STATUS(status);
    CHECK_PTR(hdr, status);
    
    f = fread_byte_swap(&hdr->nl, sizeof(int), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->ne, sizeof(int), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->d_type, sizeof(int), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->d_size, sizeof(int), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    hdr->field_labels = NULL;
    
    if (hdr->magic == SIO_MAGIC_EXTENDED)
    {
        hdr->field_labels = calloc(hdr->ne, sizeof(char *));
        CHECK_PTR(hdr->field_labels, status);
        
        for(label = 0; label < hdr->ne; label++)
        {
            f = fread_byte_swap(&tmp, sizeof(int), 1, hdr->fp, hdr->do_byte_swap, status);
            CHECK_BYTES(1, f, status);
            
            f = fread_byte_swap(&label_length, sizeof(int), 1, hdr->fp, hdr->do_byte_swap, status);
            CHECK_BYTES(1, f, status);
            
            hdr->field_labels[label] = calloc(label_length, sizeof(char));
            CHECK_PTR(hdr->field_labels[label], status);
            
            
            f = fread(hdr->field_labels[label], sizeof(char), label_length, hdr->fp);
            CHECK_BYTES(label_length, f, status);
            if (status->debug >= 10 ) {
                printf("Read in label %s\n", hdr->field_labels[label]);
            }
            
            f = fread_byte_swap(&tmp, sizeof(int), 1, hdr->fp, hdr->do_byte_swap, status);
            CHECK_BYTES(1, f, status);
            
        }
        
        f = fread_byte_swap(&tmp, sizeof(int), 1, hdr->fp, hdr->do_byte_swap, status);
        CHECK_BYTES(1, f, status);
        
    }
    
    hdr->data_start = ftello(hdr->fp);
    
    return status;
}

SPStatus *
im_close_sio(SIOheader * hdr, SPStatus * status)
{
    CHECK_STATUS(status);
    CHECK_PTR(hdr, status);
    CHECK_FP(hdr->fp, status);
    
    fclose(hdr->fp);
    
    hdr->fp = NULL;
    
    return status;
    
}

SPStatus *
im_destroy_sio(SIOheader * hdr, SPStatus * status)
{
    int label;
    
    CHECK_STATUS(status);
    CHECK_PTR(hdr, status);
    
    if (hdr->fp)
    {
        im_close_sio(hdr, status);
        CHECK_STATUS(status);
    }
    
    if (hdr->field_labels != NULL)
    {
        for (label = 0; label < hdr->ne; label++)
        {
            free(hdr->field_labels[label]);
        }
        free(hdr->field_labels);
    }
    
    hdr->field_labels = NULL;
    
    hdr->magic = 0;
    hdr->nl = 0;
    hdr->ne = 0;
    hdr->d_type = 0;
    hdr->d_size = 0;
    hdr->do_byte_swap = 0;
    hdr->data_start = 0;
    
    return status;
}

SPStatus *
im_load_sio_subset(SPImage * data, SIOheader * hdr, int startx, int starty, SPStatus * status)
{
    SPImageType itype;
    off_t first_sample;
    size_t ele_size;
    int y;
    CHECK_STATUS(status);
    CHECK_PTR(hdr, status);
    CHECK_PTR(data, status);
    
    im_sio_dtype_to_itype(hdr, &itype, status);
    CHECK_STATUS(status);
    
    if (data->nx + startx > hdr->ne)
    {
        fprintf(stderr, "Trying to read from outside of bounds of image (in X)! (%d + %d > %d)", (int) data->nx,
                (int) startx, (int) hdr->ne);
        status->status = OUTPUT_NX_MISMATCHED;
        return (status);
    }
    
    if (data->ny + starty > hdr->nl)
    {
        fprintf(stderr, "Trying to read from outside of bounds of image (in Y)! (%d + %d > %d)", (int) data->ny,
                (int) starty, (int) hdr->nl);
        status->status = OUTPUT_NY_MISMATCHED;
        return (status);
    }
    
    ele_size = im_getsizeoftype(itype);
    first_sample = hdr->data_start + (starty * hdr->ne + startx) * ele_size;
    fseeko(hdr->fp, first_sample, SEEK_SET);
    
    switch (itype)
    {
        case ITYPE_DOUBLE:
        {
            for(y = 0; y < data->ny; y++)
            {
                fread_byte_swap(&data->data.d[y * data->nx], sizeof(double), data->nx, hdr->fp, hdr->do_byte_swap, status);
                fseeko(hdr->fp, (startx + (hdr->ne - startx - data->nx)) * sizeof(double), SEEK_CUR);
            }
        }
            break;
        case ITYPE_CMPL_FLOAT:
        {
            for(y = 0; y < data->ny; y++)
            {
                fread_byte_swap(&data->data.f[y * data->nx * 2], sizeof(float), data->nx * 2, hdr->fp, hdr->do_byte_swap, status);
                fseeko(hdr->fp, (startx + (hdr->ne - startx - data->nx)) * sizeof(float) * 2, SEEK_CUR);
            }
        }
            break;
            
        case ITYPE_UINT8:
        {
            for(y = 0; y < data->ny; y++)
            {
                fread_byte_swap(&data->data.ui8[y * data->nx], sizeof(uint8_t), data->nx, hdr->fp, hdr->do_byte_swap, status);
                fseeko(hdr->fp, (startx + (hdr->ne - startx - data->nx)) * sizeof(uint8_t), SEEK_CUR);
            }
        }
            break;
            
        default:
            break;
    }
    
    return status;
}

SPStatus *
im_load_sio(SPImage * data, SIOheader * hdr, SPStatus * status)
{
    SPImageType itype;
    
    CHECK_STATUS(status);
    CHECK_PTR(hdr, status);
    CHECK_PTR(data, status);
    
    im_sio_dtype_to_itype(hdr, &itype, status);
    CHECK_STATUS(status);
    
    if (data->nx == 0 && data->ny == 0 && data->data.v == NULL)
    {
        im_create(data, itype, hdr->ne, hdr->nl, 1.0, 1.0, status);
        CHECK_STATUS(status);
    }
    
    if (data->nx != hdr->ne || data->ny != hdr->nl || data->image_type != itype)
    {
        printf("Input data image is wrong size (%ld, %ld) != (%d, %d) or wrong type (%d) != (%d)\n",
               (long)data->nx, (long)data->ny,
               hdr->ne, hdr->nl,
               data->image_type, itype);
        status->status = BAD_IMAGE;
        return status;
    }
    
    im_load_sio_subset(data, hdr, 0, 0, status);
    
    CHECK_STATUS(status);
    
    return status;
}

SPStatus *
im_sio_dtype_to_itype(SIOheader * hdr, SPImageType * itype, SPStatus * status)
{
    CHECK_STATUS(status);
    CHECK_PTR(itype, status);
    
    *itype = ITYPE_UNKNOWN;
    
    switch(hdr->d_type)
    {
        case 1:
            if (hdr->d_size == 1) {
                *itype = ITYPE_UINT8;
            }
            break;
        case 3:
            if (hdr->d_size == 8)
            {
                *itype = ITYPE_DOUBLE;
            }
            else if (hdr->d_size == 4)
            {
                *itype = ITYPE_FLOAT;
            }
            break;
            
        case 0x00d:
            if (hdr->d_size == 8)
            {
                *itype = ITYPE_CMPL_FLOAT;
            }
            break;
            
        default:
            *itype = ITYPE_UNKNOWN;
            break;
    }
    
    if (*itype == ITYPE_UNKNOWN)
    {
        status->status = INVALID_TYPE;
        printf("Unknown type (%d) and size (%d) in SIO file\n", hdr->d_type, hdr->d_size);
    }
    
    return status;
}

