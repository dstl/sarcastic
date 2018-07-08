/***************************************************************************
 *
 *       Module:    dataio_gdal.c
 *      Program:    sarclib
 *   Created by:    Matt Nottingham on 01/01/2011.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains functions for reading and writing to and from
 *      GDAL files.
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

#include <cpl_conv.h>
#include <gdal.h>
#include <ogr_api.h>
#include <ogr_srs_api.h>
#include "dataio_gdal.h"
#include "sarclib.h"

SPStatus * im_open_gdal(SPImage *a, const char * fname, SPStatus * status)
{
    GDALDatasetH hDataset;
    GDALRasterBandH hBand;
    GDALDataType gdal_data_type_on_disk;
    
    CHECK_STATUS(status);
    GDALAllRegister();
    hDataset = GDALOpen(fname, GA_ReadOnly);
    if (hDataset == NULL) {
        fprintf(stderr, "Failed to Gdal open  %s  (%s:%d)\n", fname, __FILE__, __LINE__);
        status->status = BAD_FILE;
        return status;
    }
    
    printf("Size is %dx%dx%d\n",
           GDALGetRasterXSize(hDataset),
           GDALGetRasterYSize(hDataset),
           GDALGetRasterCount(hDataset));
    
    hBand = GDALGetRasterBand(hDataset, 1);
    if (hBand == NULL) {
        fprintf(stderr, "Failed to get raster band (1) for %s  (%s:%d)\n", fname, __FILE__, __LINE__);
        status->status = BAD_FILE;
        GDALClose(hDataset);
        return status;
    }
    
    gdal_data_type_on_disk = GDALGetRasterDataType(hBand);
    im_init(a, status);
    a->nx = GDALGetRasterXSize(hDataset);
    a->ny = GDALGetRasterYSize(hDataset);
    a->image_type = im_conv_from_gdal_type(gdal_data_type_on_disk);
    if (a->image_type == ITYPE_UNKNOWN) {
        fprintf(stderr, "Unknown data type (%d) in %s  (%s:%d)\n", gdal_data_type_on_disk, fname, __FILE__, __LINE__);
        status->status = BAD_FILE;
    }
    
    // Once we're done, close properly the dataset
    //
    GDALClose(hDataset);
    
    return status;
}

SPStatus * im_load_gdal(SPImage *a, const char * fname, int band, SPStatus * status)
{
    CHECK_STATUS(status);
    
    im_open_gdal(a, fname, status);
    
    CHECK_STATUS(status);
    im_create(a, a->image_type, a->nx, a->ny, 1.0, 1.0, status);
    
    CHECK_STATUS(status);
    status = im_load_gdal_subset(a, fname, 0, 0, band, status);
    CHECK_STATUS(status);
    
    return status;
}

SPStatus * im_load_gdal_subset(SPImage *a, const char * fname, int64_t ox, int64_t oy, int band, SPStatus * status)
{
    GDALDatasetH hDataset;
    GDALRasterBandH hBand;
    GDALDataType gdal_data_type_on_disk;
    
    CHECK_STATUS(status);
    GDALAllRegister();
    hDataset = GDALOpen(fname, GA_ReadOnly);
    if (hDataset == NULL) {
        fprintf(stderr, "Failed to Gdal open  %s  (%s:%d)\n", fname, __FILE__, __LINE__);
        status->status = BAD_FILE;
        return status;
    }
    
    if (band > GDALGetRasterCount(hDataset)) {
    }
    
    hBand = GDALGetRasterBand(hDataset, band);
    
    gdal_data_type_on_disk = GDALGetRasterDataType(hBand);
    
    (void)GDALRasterIO(hBand, GF_Read, (int)ox, (int)oy, (int)a->nx, (int)a->ny, a->data.v, (int)a->nx, (int)a->ny, im_conv_to_gdal_type(a->image_type), 0, 0);
    GDALClose(hDataset);
    
    return status;
}

SPStatus * im_save_gdal(SPImage *a, const char * fname, const char * gdal_driver, SPStatus * status)
{
    GDALRasterBandH hBand;
    char **papszOptions = NULL;
    GDALDatasetH hDstDS;
    GDALDriverH hDriver;
    GDALDataType gdal_data_type;
    CPLErr err;
    
    CHECK_STATUS(status);
    gdal_data_type = im_conv_to_gdal_type(a->image_type);
    if (gdal_data_type == GDT_Unknown) {
        fprintf(stderr, "Unsupported data type (%d) for %s  (%s:%d)\n", a->image_type, fname, __FILE__, __LINE__);
        status->status = INVALID_TYPE;
        return status;
    }
    
    GDALAllRegister();
    hDriver = GDALGetDriverByName( gdal_driver );
    hDstDS  = GDALCreate(hDriver, fname, (int)a->nx, (int)a->ny, 1, gdal_data_type, papszOptions);
    
    if (hDstDS == NULL) {
        fprintf(stderr, "GDAL create failed for %s  (%s:%d)\n", fname, __FILE__, __LINE__);
        status->status = NULL_POINTER;
        return status;
    }
    hBand = GDALGetRasterBand(hDstDS, 1);
    err = GDALRasterIO(hBand, GF_Write, 0, 0, (int)a->nx, (int)a->ny,
                       a->data.v, (int)a->nx, (int)a->ny, gdal_data_type, 0, 0);
    if (err == CE_Failure || err == CE_Fatal) {
        fprintf(stderr, "GDALRastioIO failed for %s  (%s:%d)\n", fname, __FILE__, __LINE__);
        status->status = BAD_FILE;
    }
    
    // Once we're done, close properly the dataset
    //
    GDALClose(hDstDS);
    
    return status;
}

SPStatus * im_save_gdal_with_geo_ecf(SPImage *a, SPVector p1, SPVector p2, SPVector p3, SPVector p4, const char * fname, const char * gdal_driver, SPStatus * status)
{
    SPVector ll1;
    SPVector ll2;
    SPVector ll3;
    SPVector ll4;
    
    CHECK_STATUS(status);
    
    ecef2latlon(&p1, &ll1, status);
    ecef2latlon(&p2, &ll2, status);
    ecef2latlon(&p3, &ll3, status);
    ecef2latlon(&p4, &ll4, status);
    
    CHECK_STATUS(status);
    im_save_gdal_with_geo_ll(a, ll1, ll2, ll3, ll4, fname, gdal_driver, status);
    
    CHECK_STATUS(status);
    return status;
}

SPStatus * im_save_gdal_with_geo_ll(SPImage *a, SPVector ll1, SPVector ll2, SPVector ll3, SPVector ll4, const char * fname, const char * gdal_driver, SPStatus * status)
{
    double adfGeoTransform[6];
    OGRSpatialReferenceH hSRS;
    char *pszSRS_WKT = NULL;
    GDALRasterBandH hBand;
    char **papszOptions = NULL;
    GDALDatasetH hDstDS;
    GDALDriverH hDriver;
    GDALDataType gdal_data_type;
    CPLErr err;
    
    CHECK_STATUS(status);
    gdal_data_type = im_conv_to_gdal_type(a->image_type);
    if (gdal_data_type == GDT_Unknown) {
        fprintf(stderr, "Unsupported data type (%d) for %s   (%s:%d)\n", a->image_type, fname, __FILE__, __LINE__);
        status->status = INVALID_TYPE;
        return status;
    }
    
    adfGeoTransform[0] = ll1.y;
    adfGeoTransform[1] = (ll2.y - ll1.y) / (double) a->nx;
    adfGeoTransform[2] = (ll4.y - ll1.y) / (double) a->ny;
    adfGeoTransform[3] = ll1.x;
    adfGeoTransform[4] = (ll2.x - ll1.x) / (double) a->nx;
    adfGeoTransform[5] = (ll4.x - ll1.x) / (double) a->ny;
    
    GDALAllRegister();
    hDriver = GDALGetDriverByName(gdal_driver);
    hDstDS  = GDALCreate(hDriver, fname, (int)a->nx, (int)a->ny, 1, gdal_data_type, papszOptions);
    
    if (hDstDS == NULL) {
        fprintf(stderr, "GDAL create failed for %s  (%s:%d)\n", fname, __FILE__, __LINE__);
        status->status = NULL_POINTER;
        return status;
    }
    err = GDALSetGeoTransform(hDstDS, adfGeoTransform);
    if (err == CE_Failure || err == CE_Fatal) {
        fprintf(stderr, "GDALSetGeoTrans failed for %s  (%s:%d)\n", fname, __FILE__, __LINE__);
        status->status = BAD_FILE;
        return status;
    }
    
    hSRS = OSRNewSpatialReference(NULL);
    OSRSetWellKnownGeogCS(hSRS, "WGS84");
    OSRExportToWkt(hSRS, &pszSRS_WKT);
    OSRDestroySpatialReference(hSRS);
    
    err = GDALSetProjection(hDstDS, pszSRS_WKT);
    if (err == CE_Failure || err == CE_Fatal) {
        fprintf(stderr, "GDALSetProjection failed for %s  (%s:%d)\n", fname, __FILE__, __LINE__);
        status->status = BAD_FILE;
        return status;
    }
    
    CPLFree(pszSRS_WKT);
    
    hBand = GDALGetRasterBand(hDstDS, 1);
    if (hBand == NULL) {
    }
    
    err = GDALRasterIO(hBand, GF_Write, 0, 0, (int)a->nx, (int)a->ny,
                       a->data.v, (int)a->nx, (int)a->ny, gdal_data_type, 0, 0);
    
    if (err == CE_Failure || err == CE_Fatal) {
        fprintf(stderr, "GDALRastioIO failed for %s  (%s:%d)\n", fname, __FILE__, __LINE__);
        status->status = BAD_FILE;
    }
    // Once we're done, close properly the dataset
    //
    GDALClose(hDstDS);
    
    return status;
}

SPImageType
im_conv_from_gdal_type(GDALDataType itype)
{
    SPImageType image_type= ITYPE_UNKNOWN;
    
    switch(itype) {
        case GDT_Byte:
            image_type = ITYPE_UINT8;
            break;
            
        case GDT_Int16:
            image_type = ITYPE_INT16;
            break;
            
        case GDT_Int32:
            image_type = ITYPE_INT32;
            break;
            
        case GDT_Float32:
            image_type = ITYPE_FLOAT;
            break;
            
        case GDT_Float64:
            image_type = ITYPE_DOUBLE;
            break;
            
        case GDT_CFloat32:
            image_type = ITYPE_CMPL_FLOAT;
            break;
            
        case GDT_CFloat64:
            image_type = ITYPE_UNKNOWN;
            break;
            
        case GDT_UInt16:
            image_type = ITYPE_UINT16;
            break;
            
        case GDT_UInt32:
            image_type = ITYPE_UINT32;
            break;
            
        case GDT_CInt16:
            image_type = ITYPE_CMPL_INT16;
            break;
            
        case GDT_CInt32:
            image_type = ITYPE_CMPL_INT32;
            break;
            
        default:
            image_type = ITYPE_UNKNOWN;
            break;
    }
    
    return (image_type);
}

GDALDataType
im_conv_to_gdal_type(SPImageType itype)
{
    GDALDataType image_type = GDT_Unknown;
    
    switch (itype) {
        case ITYPE_POLAR:
            image_type = GDT_Unknown;
            break;
            
        case ITYPE_CMPL_FLOAT:
            image_type = GDT_CFloat32;
            break;
            
        case ITYPE_FLOAT:
            image_type = GDT_Float32;
            break;
            
        case ITYPE_DOUBLE:
            image_type = GDT_Float64;
            break;
            
        case ITYPE_INT64:
            image_type = GDT_Unknown;
            break;
            
        case ITYPE_INT32:
            image_type = GDT_Int32;
            break;
            
        case ITYPE_INT16:
            image_type = GDT_Int16;
            break;
            
        case ITYPE_INT8:
            image_type = GDT_Unknown;
            break;
            
        case ITYPE_UINT64:
            image_type = GDT_Unknown;
            break;
            
        case ITYPE_UINT32:
            image_type = GDT_UInt32;
            break;
            
        case ITYPE_UINT16:
            image_type = GDT_UInt16;
            break;
            
        case ITYPE_UINT8:
            image_type = GDT_Byte;
            break;
            
        case ITYPE_VECTOR:
            image_type = GDT_Unknown;
            break;
            
        case ITYPE_CMPL_INT64:
            image_type = GDT_Unknown;
            break;
            
        case ITYPE_CMPL_INT32:
            image_type = GDT_CInt32;
            break;
            
        case ITYPE_CMPL_INT16:
            image_type = GDT_CInt16;
            break;
            
        case ITYPE_CMPL_INT8:
            image_type = GDT_Unknown;
            break;
            
        default:
            image_type = GDT_Unknown;
            break;
    }
    
    return (image_type);
}
