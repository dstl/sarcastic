/***************************************************************************
 * 
 *           Module :  dataio_srf.c
 *          Program :  sarclib
 *       Created by :  Matt Nottingham on 18/08/2008
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      Description: Files for reading/writing SRF format files.
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

static SPStatus * im_read_header_srf(SRFheader * hdr, SPStatus * status);

SPStatus *
im_open_srf(SRFheader * hdr, const char * filename, SPStatus * status)
{
    char tmp_buffer[SRF_RECORD_LENGTH + 1];
    
    CHECK_STATUS(status);
    CHECK_PTR(hdr, status);
    
    hdr->fp = fopen(filename, "r");
    CHECK_FP(hdr->fp, status);
    
    if (fread(tmp_buffer, sizeof(char), SRF_RECORD_LENGTH, hdr->fp) != SRF_RECORD_LENGTH) {
        perror("Short read on opening SRF file!");
        status->status = BAD_FILE;
        return (status);
    }
    
    tmp_buffer[SRF_RECORD_LENGTH] = 0;
    
    if ((strcasecmp(tmp_buffer, "srf version 1.2") != 0) && (strcasecmp(tmp_buffer, "srf version v1.2") != 0)) {
        perror("Bad header reading in SRF file");
        status->status = BAD_FILE;
        return (status);
    }
    
    im_read_header_srf(hdr, status);
    
    return status;
}

static SPStatus *
im_read_header_srf(SRFheader * hdr, SPStatus * status)
{
    char tmp_buffer[SRF_RECORD_LENGTH + 1];
    size_t f;
    char dtype[8];
    
    CHECK_STATUS(status);
    CHECK_PTR(hdr, status);
    
    if (fread(tmp_buffer, sizeof(char), SRF_RECORD_LENGTH, hdr->fp) != SRF_RECORD_LENGTH) {
        perror("Short read with security record in SRF file!");
        status->status = BAD_FILE;
        return (status);
    }
    
    if (fread(hdr->datetime, sizeof(char), SRF_RECORD_LENGTH, hdr->fp) != SRF_RECORD_LENGTH) {
        perror("Short read with date/time record in SRF file!");
        status->status = BAD_FILE;
        return (status);
    }
    
    if (fread(tmp_buffer, sizeof(char), SRF_RECORD_LENGTH, hdr->fp) != SRF_RECORD_LENGTH) {
        perror("Short read with creation info record in SRF file!");
        status->status = BAD_FILE;
        return (status);
    }
    
    if (fread(tmp_buffer, sizeof(char), SRF_RECORD_LENGTH, hdr->fp) != SRF_RECORD_LENGTH) {
        perror("Short read with creation(2) info record in SRF file!");
        status->status = BAD_FILE;
        return (status);
    }
    
    if (fread(tmp_buffer, sizeof(char), SRF_RECORD_LENGTH, hdr->fp) != SRF_RECORD_LENGTH) {
        perror("Short read with mode record in SRF file!");
        status->status = BAD_FILE;
        return (status);
    }
    
    
    if (fread(tmp_buffer, sizeof(char), SRF_RECORD_LENGTH, hdr->fp) != SRF_RECORD_LENGTH) {
        perror("Short read with date record in SRF file!");
        status->status = BAD_FILE;
        return (status);
    }
    
    if (fread(dtype, sizeof(char), 8, hdr->fp) != 8) {
        perror("Short read with date record in SRF file!");
        status->status = BAD_FILE;
        return (status);
    }
    
    hdr->data_type = ITYPE_UNKNOWN;
    
    if (strcmp(dtype, "BYTE") == 0) {
        hdr->data_type = ITYPE_UINT8;
    }
    
    if (strcmp(dtype, "C8") == 0) {
        hdr->data_type = ITYPE_CMPL_FLOAT;
    }
    
    if (strcmp(dtype, "R4") == 0) {
        hdr->data_type = ITYPE_FLOAT;
    }
    
    if (hdr->data_type == ITYPE_UNKNOWN) {
        printf("SRF files doesn't contain BYTE, C8 or R4 data! (%s)\n\n", dtype);
        status->status = BAD_FILE;
        return (status);
    }
    
    if (fread(tmp_buffer, sizeof(char), 7, hdr->fp) != 7) {
        perror("Short read with endian field in SRF file!");
        status->status = BAD_FILE;
        return (status);
    }
    
    if (strcmp(tmp_buffer, "LITTLE") == 0) {
        if (im_machine_type() == IM_LITTLE_ENDIAN) {
            hdr->do_byte_swap = 0;
        } else {
            hdr->do_byte_swap = 1;
        }
    } else {
        if (im_machine_type() == IM_LITTLE_ENDIAN) {
            hdr->do_byte_swap = 1;
        } else {
            hdr->do_byte_swap = 0;
        }
    }
        
    fseeko(hdr->fp, 952, SEEK_SET);
    
    f = fread_byte_swap(&hdr->ncols, sizeof(int), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&hdr->nrows, sizeof(int), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&hdr->sx, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&hdr->sy, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&hdr->resX, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&hdr->resY, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    fseeko(hdr->fp, 1008, SEEK_SET);
    f = fread_byte_swap(&hdr->grp.x, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->grp.y, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->grp.z, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&hdr->img_grp, sizeof(int), 2, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(2, f, status);
    
    f = fread_byte_swap(&hdr->ipn.x, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->ipn.y, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->ipn.z, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&hdr->ip_rng.x, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->ip_rng.y, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->ip_rng.z, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&hdr->ip_cr.x, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->ip_cr.y, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->ip_cr.z, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&hdr->fpn.x, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->fpn.y, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->fpn.z, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    if (ftell(hdr->fp) != 1136L) {
        printf("We screwed up\n");
    }
    
    fseeko(hdr->fp, 1136, SEEK_SET);
    f = fread_byte_swap(&hdr->apc.x, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->apc.y, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->apc.z, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->apc_tx.x, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->apc_tx.y, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->apc_tx.z, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->apc_rx.x, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->apc_rx.y, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->apc_rx.z, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&hdr->syn.x, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->syn.y, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->syn.z, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&hdr->vel.x, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->vel.y, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->vel.z, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    printf("Vel = %.12g %.12g %.12g\n", hdr->vel.x, hdr->vel.y, hdr->vel.z);
    
    f = fread_byte_swap(&hdr->dof[0], sizeof(double), 2, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(2, f, status);
    
    f = fread_byte_swap(&hdr->ll1.x, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->ll1.y, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&hdr->ll2.x, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->ll2.y, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&hdr->ll3.x, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->ll3.y, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&hdr->ll4.x, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->ll4.y, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&hdr->datum, sizeof(char), 16, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(16, f, status);
    
    f = fread_byte_swap(&hdr->lambda[0], sizeof(double), 2, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(2, f, status);
    
    f = fread_byte_swap(&hdr->num_pulse, sizeof(int), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->pulse_offset, sizeof(int), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->num_samp, sizeof(int), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    f = fread_byte_swap(&hdr->samp_offset, sizeof(int), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&hdr->fft_dim, sizeof(int), 2, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(2, f, status);
    
    f = fread_byte_swap(&hdr->sl_ratio, sizeof(double), 2, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(2, f, status);
    
    f = fread_byte_swap(&hdr->rect_dim, sizeof(int), 2, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(2, f, status);
    
    f = fread_byte_swap(&hdr->rect0, sizeof(double), 2, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(2, f, status);
    
    f = fread_byte_swap(&hdr->rect_sp, sizeof(double), 2, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(2, f, status);
    
    printf("rect0 %.12g, %.12g   rect_sp %.12g,%.12g\n", hdr->rect0[0], hdr->rect0[1], hdr->rect_sp[0], hdr->rect_sp[1]);
    
    f = fread_byte_swap(&hdr->pulse_duration, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&hdr->tx_bw, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    f = fread_byte_swap(&hdr->rx_bw, sizeof(double), 1, hdr->fp, hdr->do_byte_swap, status);
    CHECK_BYTES(1, f, status);
    
    
    
    hdr->data_start = 2312L;
    
    return status;
}

SPStatus *
im_close_srf(SRFheader * hdr, SPStatus * status)
{
    CHECK_STATUS(status);
    CHECK_PTR(hdr, status);
    CHECK_FP(hdr->fp, status);
    
    fclose(hdr->fp);
    
    hdr->fp = NULL;
    
    return status;
    
}

SPStatus *
im_destroy_srf(SRFheader * hdr, SPStatus * status)
{
    CHECK_STATUS(status);
    CHECK_PTR(hdr, status);
    
    if (hdr->fp)
    {
        im_close_srf(hdr, status);
        CHECK_STATUS(status);
    }
    
    hdr->nrows = 0;
    hdr->ncols = 0;
    hdr->data_start = 0;
    
    return status;
}

SPStatus *
im_load_srf_subset(SPImage * data, SRFheader * hdr, int startx, int starty, SPStatus * status)
{
    off_t first_sample;
    size_t ele_size;
    int y;
    CHECK_STATUS(status);
    CHECK_PTR(hdr, status);
    CHECK_PTR(data, status);
    
    if (data->image_type != hdr->data_type) {
        fprintf(stderr, "Image and SRF file are not of same type! (%d != %d)\n", data->image_type, hdr->data_type);
        status->status = TYPES_NOT_SAME;
        return(status);
    }
    
    if (data->nx + startx > hdr->ncols)
    {
        fprintf(stderr, "Trying to read from outside of bounds of image (in X)! (%d + %d > %d)", (int) data->nx,
                (int) startx, (int) hdr->ncols);
        status->status = OUTPUT_NX_MISMATCHED;
        return (status);
    }
    
    if (data->ny + starty > hdr->nrows)
    {
        fprintf(stderr, "Trying to read from outside of bounds of image (in Y)! (%d + %d > %d)", (int) data->ny,
                (int) starty, (int) hdr->nrows);
        status->status = OUTPUT_NY_MISMATCHED;
        return (status);
    }
    
    ele_size = im_getsizeoftype(data->image_type);
    first_sample = hdr->data_start + (starty * hdr->ncols + startx) * ele_size;
    fseeko(hdr->fp, first_sample, SEEK_SET);
    
    if (hdr->data_type == ITYPE_CMPL_FLOAT) {
        for(y = 0; y < data->ny; y++)
        {
            fread_byte_swap(&data->data.f[y * data->nx * 2], sizeof(float), data->nx * 2, hdr->fp, hdr->do_byte_swap, status);
            fseeko(hdr->fp, (startx + (hdr->ncols - startx - data->nx)) * sizeof(float) * 2, SEEK_CUR);
        }
    } else {
        for(y = 0; y < data->ny; y++)
        {
            fread_byte_swap(&data->data.f[y * data->nx], sizeof(float), data->nx, hdr->fp, hdr->do_byte_swap, status);
            fseeko(hdr->fp, (startx + (hdr->ncols - startx - data->nx)) * sizeof(float), SEEK_CUR);
        }
    }
    
    return status;
}

SPStatus *
im_load_srf(SPImage * data, SRFheader * hdr, SPStatus * status)
{
    CHECK_STATUS(status);
    CHECK_PTR(hdr, status);
    CHECK_PTR(data, status);
    
    CHECK_STATUS(status);
    
    if (data->nx == 0 && data->ny == 0 && data->data.v == NULL)
    {
        im_create(data, hdr->data_type, hdr->ncols, hdr->nrows, hdr->sx, hdr->sy, status);
        CHECK_STATUS(status);
    }
    
    if (data->nx != hdr->ncols || data->ny != hdr->nrows || data->image_type != hdr->data_type)
    {
        printf("Input data image is wrong size (%ld, %ld) != (%d, %d) or wrong type (%d) != (%d)\n",
               (long)data->nx, (long)data->ny,
               hdr->ncols, hdr->nrows,
               data->image_type, hdr->data_type);
        status->status = BAD_IMAGE;
        return status;
    }
    
    im_load_srf_subset(data, hdr, 0, 0, status);
    
    CHECK_STATUS(status);
    
    return status;
}
