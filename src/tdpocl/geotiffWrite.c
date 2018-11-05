/***************************************************************************
 * 
 *           Module :  geotiffWrite.c
 *          Program :  tdpocl
 *       Created by :  Darren Muff on 01/06/2013
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *          Function to write image out as a Geotiff
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
#include <stdio.h>
#include <mtlib/mtlib.h>
#include <gdal.h>

void geotiffWrite( SPImage *im, char *filename, SPVector corners[4], SPStatus *status ){
    
	GDALAllRegister();
	
    const char *pszFormat = "GTiff";
    GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
    
    GDALDatasetH hDstDS;
    char **papszOptions = NULL;
    
//    hDstDS = GDALCreate(hDriver, filename, (int)im->nx, (int)im->ny, 1, GDT_Float32, papszOptions);
    hDstDS = GDALCreate(hDriver, filename, (int)im->nx, (int)im->ny, 1, GDT_Int16, papszOptions);

    
//    adfGeoTransform[0] /* top left x */
//    adfGeoTransform[1] /* w-e pixel resolution */
//    adfGeoTransform[2] /* rotation, 0 if image is "north up" */
//    adfGeoTransform[3] /* top left y */
//    adfGeoTransform[4] /* rotation, 0 if image is "north up" */
//    adfGeoTransform[5] /* n-s pixel resolution */
//      or, in human form....
//    padfTransform[0] = corner lon
//    padfTransform[1] = cos(alpha)*(scaling)
//    padfTransform[2] = -sin(alpha)*(scaling)
//    padfTransform[3] = corner lat
//    padfTransform[4] = sin(alpha)*(scaling)
//    padfTransform[5] = cos(alpha)*(scaling)
    
    SPVector cornersLL[4] ;
    for(int i=0; i<4; i++){
        ecef2latlon(&(corners[i]), &(cornersLL[i]), status) ;
    }
    
    double adfGeoTransform[6] ;
    double theta ; // rotation angle of image
    SPVector Z, E, B; // Z is Earth axis north up; E is local east; B is bearing vector
    VECT_CREATE(0, 0, 1, Z) ;
    VECT_CROSS(Z, corners[0], E);
    VECT_UNIT(E, E) ;
    VECT_SUB(corners[1], corners[0], B) ;
    theta = acos(VECT_DOT(E, B));
    
    adfGeoTransform[0] = cornersLL[0].y ;
    adfGeoTransform[1] =  cos(theta)*im->xspc ;
    adfGeoTransform[2] = -sin(theta)*im->xspc ;
    adfGeoTransform[3] =  cornersLL[0].x ;
    adfGeoTransform[4] =  sin(theta)*im->yspc ;
    adfGeoTransform[5] =  cos(theta)*im->yspc ;
    
    GDALSetGeoTransform( hDstDS, adfGeoTransform );
    
    GDALRasterBandH hBand;
//    float *abyRaster;
    short *abyRaster ;
//    abyRaster = (float *)calloc(sizeof(float), im->nx*im->ny);
    abyRaster = (short *)calloc(sizeof(short), im->nx*im->ny);

    if(abyRaster == NULL){
//        printf("ERROR: Failed to calloc %lld bytes for tiff output\n",sizeof(float)*im->nx*im->ny);
        printf("ERROR: Failed to calloc %lld bytes for tiff output\n",sizeof(short)*im->nx*im->ny);

        exit(1);
    }
    
    // load data into gdal output array
    //
    SPCmplx t;
    for(int y=0; y<im->ny; y++){
        for(int x=0; x<im->nx; x++){
            t = im->data.cmpl_f[y*im->nx+x] ;
//            abyRaster[y*im->nx+x] = CMPLX_MAG(t) ;
            abyRaster[y*im->nx+x] = y;

        }
    }
    
    hBand = GDALGetRasterBand( hDstDS, 1 );
//    GDALRasterIO( hBand, GF_Write, 0, 0, (int)im->nx, (int)im->ny,abyRaster, (int)im->nx, (int)im->ny, GDT_Float32, 0, 0 );
    GDALRasterIO( hBand, GF_Write, 0, 0, (int)im->nx, (int)im->ny,abyRaster, (int)im->nx, (int)im->ny, GDT_Int16, 0, 0 );
    
    /* Once we're done, close properly the dataset */
    GDALClose( hDstDS );
    
    return ; 
}

