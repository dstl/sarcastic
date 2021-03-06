/***************************************************************************
 * 
 *           Module :  dataio_gdal_stubs.c
 *          Program :  sarclib
 *       Created by :  Matt Nottingham on 01/01/2011
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      This file contains the stubs for the GDAL functions if GDAL
 *      isn't available.
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

#include "dataio_gdal.h"

SPStatus * im_open_gdal(SPImage *a, const char * fname, SPStatus * status)
{
    fprintf(stderr, "\n***GDAL support was not compiled into sarclib***\n\n");
    status->status = BAD_FILE;
    return status;
}

SPStatus * im_load_gdal(SPImage *a, const char * fname, int band, SPStatus * status)
{
    fprintf(stderr, "\n***GDAL support was not compiled into sarclib***\n\n");
    status->status = BAD_FILE;
    return status;
}

SPStatus * im_load_gdal_subset(SPImage *a, const char * fname, int64_t ox, int64_t oy, int band, SPStatus * status)
{
    fprintf(stderr, "\n***GDAL support was not compiled into sarclib***\n\n");
    status->status = BAD_FILE;
    
    return status;
}

SPStatus * im_save_gdal(SPImage *a, const char * fname, const char * gdal_driver, SPStatus * status)
{
    fprintf(stderr, "\n***GDAL support was not compiled into sarclib***\n\n");
    status->status = BAD_FILE;
    
    return status;
}

SPStatus * im_save_gdal_with_geo_ecf(SPImage *a, SPVector p1, SPVector p2, SPVector p3, SPVector p4, const char * fname, const char * gdal_driver, SPStatus * status)
{
    fprintf(stderr, "\n***GDAL support was not compiled into sarclib***\n\n");
    status->status = BAD_FILE;
    
    return status;
}

SPStatus * im_save_gdal_with_geo_ll(SPImage *a, SPVector ll1, SPVector ll2, SPVector ll3, SPVector ll4, const char * fname, const char * gdal_driver, SPStatus * status)
{
    fprintf(stderr, "\n***GDAL support was not compiled into sarclib***\n\n");
    status->status = BAD_FILE;
    
    return status;
}

SPImageType
im_conv_from_gdal_type(GDALDataType itype)
{
    fprintf(stderr, "\n***GDAL support was not compiled into sarclib***\n\n");
    
    return (ITYPE_UNKNOWN);
}

GDALDataType
im_conv_to_gdal_type(SPImageType itype)
{
    fprintf(stderr, "\n***GDAL support was not compiled into sarclib***\n\n");
    
    return (0);
}
