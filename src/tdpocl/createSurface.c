/***************************************************************************
 *
 *       Module:    createSurface.c
 *      Program:    tdp2
 *   Created by:    Darren Muff on 15/03/2013.
 *                  Copyright (c) 2013 [Dstl]. All rights reserved.
 *
 *   Description:
 *   <ENTER DESCRIPTION HERE>
 Creates a surface file and saves it to disk
 *
 *
 *   CLASSIFICATION        :  <PENDING>
 *   Date of CLASSN        :  15/03/2013
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

#include <stdio.h>
#include "tdp.h"

int create_surface(SPVector sc_pos, SPVector sensor_pos, char * surface_fname, int nx, int ny, double sx, double sy, SPStatus * status)
{
    int x;
    int y;
    double azimuth;
    SPImage surface;
    
    SPVector R;
    SPVector R_;
    SPVector KDP_;
    SPVector JDP_;
    SPVector IDP_;
    
    SPVector rng_dir;
    SPVector azi_dir;
    SPImageLineSaveInfo info;
    
    printf("Sensor position       : %f %f %f\n", sensor_pos.x, sensor_pos.y, sensor_pos.z);
    printf("Scene centre position : %f %f %f\n", sc_pos.x, sc_pos.y, sc_pos.z);
    im_create(&surface, ITYPE_VECTOR, nx, 1, sx, sy, status);
    im_save_line_init(&surface, surface_fname, &info, status);
    
    VECT_SUB(sensor_pos, sc_pos, R);
    VECT_UNIT(R, R_);
    
    /* IDP_, KDP_ & IDP_ are unit vectors for the aim point, with JDP in the range direction and IDP in the azimuth direction */
    VECT_UNIT(sc_pos, KDP_);         // == zbar
    VECT_PROJ(R_, KDP_, JDP_);
    VECT_CROSS(JDP_, KDP_, IDP_);
    
    azimuth = -atan2(IDP_.z, JDP_.z);
    if (azimuth < 0) {
        azimuth += 2.0 * M_PI;
    }
    
    printf("Azimuth angle %g, Grazing angle %g\n", azimuth * 180.0 / M_PI, asin ( VECT_DOT(R_, KDP_))* 180.0 / M_PI); /* Print it out just to check I got it OK */
    
    VECT_SCMULT(JDP_, -sy, rng_dir);
    VECT_SCMULT(IDP_,  sx, azi_dir);
    
    SPVector itl, itr, ibl, ibr, otl, otr, obl, obr ;
    itl.x = surface.data.vect[0].x = sc_pos.x + (0 - surface.nx/2) * azi_dir.x + (0 - ny/2) * rng_dir.x;
    itl.y = surface.data.vect[0].y = sc_pos.y + (0 - surface.nx/2) * azi_dir.y + (0 - ny/2) * rng_dir.y;
    itl.z = surface.data.vect[0].z = sc_pos.z + (0 - surface.nx/2) * azi_dir.z + (0 - ny/2) * rng_dir.z;
    itr.x = surface.data.vect[nx-1].x = sc_pos.x + (nx-1 - surface.nx/2) * azi_dir.x + (0 - ny/2) * rng_dir.x;
    itr.y = surface.data.vect[nx-1].y = sc_pos.y + (nx-1 - surface.nx/2) * azi_dir.y + (0 - ny/2) * rng_dir.y;
    itr.z = surface.data.vect[nx-1].z = sc_pos.z + (nx-1 - surface.nx/2) * azi_dir.z + (0 - ny/2) * rng_dir.z;
    ibl.x = surface.data.vect[0].x = sc_pos.x + (0 - surface.nx/2) * azi_dir.x + (ny-1 - ny/2) * rng_dir.x;
    ibl.y = surface.data.vect[0].y = sc_pos.y + (0 - surface.nx/2) * azi_dir.y + (ny-1 - ny/2) * rng_dir.y;
    ibl.z = surface.data.vect[0].z = sc_pos.z + (0 - surface.nx/2) * azi_dir.z + (ny-1 - ny/2) * rng_dir.z;
    ibr.x = surface.data.vect[nx-1].x = sc_pos.x + (nx-1 - surface.nx/2) * azi_dir.x + (ny-1 - ny/2) * rng_dir.x;
    ibr.y = surface.data.vect[nx-1].y = sc_pos.y + (nx-1 - surface.nx/2) * azi_dir.y + (ny-1 - ny/2) * rng_dir.y;
    ibr.z = surface.data.vect[nx-1].z = sc_pos.z + (nx-1 - surface.nx/2) * azi_dir.z + (ny-1 - ny/2) * rng_dir.z;
    ecef2latlon(&itl, &otl, status);
    ecef2latlon(&itr, &otr, status);
    ecef2latlon(&ibl, &obl, status);
    ecef2latlon(&ibr, &obr, status);
    printf("Image Definition (lat deg, lon deg, alt m): \n");
    printf("   %2.5f,%2.5f,%2.1f  --------  %2.5f,%2.5f,%2.1f\n",otl.lat,otl.lon,otl.alt,otr.lat,otr.lon,otr.alt);
    printf("           |                            |\n");
    printf("           |                            |\n");
    printf("   %2.5f,%2.5f,%2.1f  --------  %2.5f,%2.5f,%2.1f\n",obl.lat,obl.lon,obl.alt,obr.lat,obr.lon,obr.alt);
    
    for(y = 0; y < ny; y++) {
        for(x = 0; x < nx; x++) {
            surface.data.vect[x].x = sc_pos.x + (x - surface.nx/2) * azi_dir.x + (y - ny/2) * rng_dir.x;
            surface.data.vect[x].y = sc_pos.y + (x - surface.nx/2) * azi_dir.y + (y - ny/2) * rng_dir.y;
            surface.data.vect[x].z = sc_pos.z + (x - surface.nx/2) * azi_dir.z + (y - ny/2) * rng_dir.z;
        }
        im_save_line(&surface, &info, status);
        if ((y % 64) == 0) {
            printf("   Done %.3g%% of the surface         \r", (double) y * 100.0 / (double) ny);
            fflush(stdout);
        }
    }
    
    im_save_line_close(&info, status);
    printf("\n");
    
    im_destroy(&surface, status);
    return(0);
}
