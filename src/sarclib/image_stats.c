/***************************************************************************
 *
 *       Module:    image_stats.c
 *      Program:    sarclib
 *   Created by:    Matt Notingham on 14/12/2005.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      Functions to compute image statistics
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

#include "sarclib.h"

// This loop is used to calculated image stats, it works out
// the min and max of the data, and their positions in the array
// it also works out the mean and the sum of the squares (to
// work out the standard deviation later) 
//
#define LOOP(type) {                        \
    info->min = a->data. type [0];          \
    info->max = info->min;                  \
    info->max_pos_x = 0;                    \
    info->max_pos_y = 0;                    \
    info->min_pos_x = 0;                    \
    info->min_pos_y = 0;	                \
                                            \
    for(i = 0; i < s; i++) {                \
        tmp = a->data. type [i];            \
        info->mean += tmp;                  \
        ssum += tmp * tmp;                  \
                                            \
        if (tmp > info->max) {              \
            info->max = tmp;                \
            info->max_pos_x = i % a->nx;    \
            info->max_pos_y = i / a->nx;    \
        }                                   \
                                            \
        if (tmp < info->min) {              \
            info->min = tmp;                \
            info->min_pos_x = i % a->nx;    \
            info->min_pos_y = i / a->nx;    \
        }                                   \
    }                                       \
}

// This macro sorts out which histogram bin a particular value belongs to then
// adds one to that bin
//
#define H_LINEARLOOP(type) {                                                    \
    for(i = 0; i < s; i++) {                                                    \
        tmp = (double) a->data. type [i];                                       \
        bin = (tmp - min_bin) * (double)(num_bins - 1)/ (max_bin - min_bin);    \
        if (bin >= 0 && bin < num_bins) {                                       \
                hist->data.i32[bin]++;                                          \
        }                                                                       \
    }                                                                           \
}

// This function calculates the min and max values of an image and the location
// of the min and max.  It also calculates the standard deviation and the mean
//
SPStatus*
im_stats(SPImage * a, SPStatsInfo * info, SPStatus * status)
{
    int64_t s;
    int64_t i;
    
    double ssum;
    double tmp;
    
    CHECK_STATUS(status);
    CHECK_IMINTEG(a, status);
    CHECK_PTR(info, status);
    
    s = a->nx * a->ny;
    
    ssum = 0.0;
    info->mean = 0.0;
    info->max_pos_x = 0;
    info->max_pos_y = 0;
    info->min_pos_x = 0;
    info->min_pos_y = 0;
    
    switch (a->image_type)
    {
        case ITYPE_DOUBLE:
            LOOP(d);
            break;
        case ITYPE_FLOAT:
            LOOP(f);
            break;
        case ITYPE_INT64:
            LOOP(i64);
            break;
        case ITYPE_INT32:
            LOOP(i32);
            break;
        case ITYPE_INT16:
            LOOP(i16);
            break;
        case ITYPE_INT8:
            LOOP(i8);
            break;
        case ITYPE_UINT64:
            LOOP(ui64);
            break;
        case ITYPE_UINT32:
            LOOP(ui32);
            break;
        case ITYPE_UINT16:
            LOOP(ui16);
            break;
        case ITYPE_UINT8:
            LOOP(ui8);
            break;
        default:
            printf("im_stats: Sorry, type %d (%s) is not supported\n", a->image_type, itype2string(a->image_type));
            status->status = INVALID_TYPE;
            return (status);
    }
    
    
    info->mean /= (double)s;
    ssum /= (double)s;
    
    // This is now the SD of the sample rather than the population - in practice it don't make a lot of difference
    //
    info->sd = sqrt(s) * sqrt(ssum - info->mean * info->mean) / sqrt(s - 1);
    
    return(status);
}

// This function makes a histogram of the data, the histogram is output as another image
//
SPStatus*
im_histo(SPImage * a, double min_bin, double max_bin, HistoMode mode, int64_t num_bins, SPImage * hist, SPStatus * status)
{
    int64_t s;
    int64_t i;
    double tmp;
    int64_t bin;
    
    CHECK_STATUS(status);
    CHECK_IMINTEG(a, status);
    CHECK_PTR(hist, status);
    
    s = a->nx * a->ny;
    
    if (hist->nx != 0 && hist->ny != 0) {
        if (hist->nx != num_bins || hist->ny != 1) {
            fprintf(stderr, "im_histo: histogram image is not the correct size (%ld, %ld) != (%ld, 1)\n",
                    (long)hist->nx, (long)hist->ny, (long)num_bins);
            status->status = OUTPUT_SIZES_MISMATCHED;
            return(status);
        }
        
        if (hist->image_type != ITYPE_INT32) {
            fprintf(stderr, "im_histo: histogram image is not correct type (%d (%s)), should be INT32\n", hist->image_type, itype2string(a->image_type));
            status->status = INVALID_TYPE;
            return (status);
        }
    }
    
    if (hist->nx == 0 && hist->ny == 0) {
        im_create(hist, ITYPE_INT32, num_bins, 1, 1.0, 1.0, status);
        CHECK_STATUS(status);
    }
    switch (mode)
    {
        case HISTO_LINEAR_BINS:
        {
            switch (a->image_type)
            {
                case ITYPE_DOUBLE:
                    H_LINEARLOOP(d);
                    break;
                case ITYPE_FLOAT:
                    H_LINEARLOOP(f);
                    break;
                case ITYPE_INT64:
                    H_LINEARLOOP(i64);
                    break;
                case ITYPE_INT32:
                    H_LINEARLOOP(i32);
                    break;
                case ITYPE_INT16:
                    H_LINEARLOOP(i16);
                    break;
                case ITYPE_INT8:
                    H_LINEARLOOP(i8);
                    break;
                case ITYPE_UINT64:
                    H_LINEARLOOP(ui64);
                    break;
                case ITYPE_UINT32:
                    H_LINEARLOOP(ui32);
                    break;
                case ITYPE_UINT16:
                    H_LINEARLOOP(ui16);
                    break;
                case ITYPE_UINT8:
                    H_LINEARLOOP(ui8);
                    break;
                default:
                    fprintf(stderr, "im_histo: Sorry, type %d (%s) is not supported\n", a->image_type, itype2string(a->image_type));
                    status->status = INVALID_TYPE;
                    return (status);
            }
        }
            break;
            
        default:
            fprintf(stderr, "im_histo: Sorry, that histogram mode (%d) is not supported\n", mode);
            status->status = INVALID_TYPE;
            break;
    }
    
    return (status);
}
