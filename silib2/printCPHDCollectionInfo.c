/***************************************************************************
 *
 *       Module:    printCPHDCollectionInfo.c
 *      Program:    SILib2
 *   Created by:    Darren on 26/07/2013.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      Print useful collection info about a CPHD file
 *
 *
 *   CLASSIFICATION        :  <PENDING>
 *   Date of CLASSN        :  14/03/2013
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

int printCPHDCollectionInfo(CPHDHeader *hdr, SPStatus *status){
    
    CHECK_STATUS(status) ;
    if ( hdr == NULL ){
        if (status->debug > 10 ) {
            printf("ERROR : CPHD File not specified in argument 1 of printCPHDCollectionInfo\n") ;
        }
        status->status = CPHD_FILE_INVALID ;
        return (status->status) ;
    }
    
    SPVector nullgrp;
    VECT_CREATE(0, 0, 0, nullgrp) ;
    collectionGeom geom1 ;
    int err ;
    
    err = collectionGeometry(hdr, -1, nullgrp, &geom1, status);
    CHECK_STATUS(status) ;
    SPVector srp_lla ;
    ecef2latlon(&geom1.SRP, &srp_lla, status);
    CHECK_STATUS(status) ;
    
    printf("\n");
    printf("                   Collection Summary                    \n");
    printf("---------------------------------------------------------\n");
    printf("Image ID                    : %15s \n",hdr->dataSetID);
    printf("Time of Collect             : %s   \n",hdr->dateTime);
    printf("Dataset Classification      : %s   \n",hdr->classification) ;
    printf("Collection mode             : %s   \n",hdr->mode) ;
    printf("Collection Geometry         : %s   \n",hdr->geometry) ;
    if (hdr->fixedSRP == 1) {
        printf("Is SRP fixed                : Yes   \n") ;
    }else if (hdr->fixedSRP == 0){
        printf("Is SRP fixed                : No   \n") ;
    }else{
        printf("Is SRP fixed                : UNKNOWN   \n") ;
    }
    printf("Sensor                      : %s \n",hdr->sensor) ;
    char data_type_str[8];
    switch (hdr->data_type){
        case ITYPE_CMPL_FLOAT:
            strcpy(data_type_str,"cmplxf");
            break;
        case ITYPE_CMPL_INT16:
            strcpy(data_type_str,"cmplxs");
            break;
        case ITYPE_CMPL_INT8:
            strcpy(data_type_str,"cmplxb");
            break;
        default:
            strcpy(data_type_str,"cmplxn");
    }
    printf("CPHD file data type         : %s \n",data_type_str) ;
    printf("Number of channels          : %d \n",hdr->nchan) ;
    if(!(hdr->polarisation == NULL)){
        for(int i=0; i<hdr->nchan; ++i){
            printf("Channel %d polarisation     : %s\n",i,hdr->polarisation[i]) ;
            printf("  Azimuth 3dB beam width    : %f\n", hdr->antenna_width_az) ;
            printf("  Elevation 3dB beam Width  : %f\n", hdr->antenna_width_el) ;
        }
    }
    if (hdr->deskewed) {
        printf("Deskew applied?             : Yes\n");
    }else{
        printf("Deskew applied?             : No\n");
    }
    printf("Tasked SRP                  : %f %f %f (deg, deg, m)\n", srp_lla.lat, srp_lla.lon, srp_lla.alt );
    printf("Centre Azimuth is           : %5.1f deg \n",geom1.azimuthRad*GeoConsts_RADTODEG);
    printf("Centre Grazing angle is     : %5.1f deg \n",geom1.grazingRad*GeoConsts_RADTODEG);
    printf("Centre Squint angle is      : %f deg \n",geom1.squintRad*GeoConsts_RADTODEG);
    printf("Centre range is             : %5.1f metres \n",geom1.range);
    printf("Number of pulses            : %5d \n",hdr->num_azi);
    printf("ADC Samples per pulse       : %5d \n",hdr->nsamp);
    printf("Chirp rate                  : %3.3e Hz/sec\n",hdr->chirp_gamma);
    printf("ADRate                      : %3.3e Samps/sec\n",hdr->clock_speed);
    printf("TX Pulse duration           : %5.6f milliseconds  \n",1000*hdr->pulse_length);
    if (geom1.look == LeftLook ){
        printf("Look direction is           : 'Left'\n");
    }else if(geom1.look == RightLook ){
        printf("Look direction is           : 'Right'\n");
    }else {
        printf("Look direction is           : 'not calculated'\n");
    }
    // We have to be a bit crafty calculating the SA duration as first and
    // last pulses are usually junk
    //
    double saStart =  hdr->pulses[0].sat_rx_time ;
    double saEnd   = hdr->pulses[hdr->num_azi-1].sat_rx_time  ;
    double meanPRF = 1/((saEnd-saStart) / (hdr->num_azi)) ;
    printf("Mean PRF                    : %7.1f Hz\n",meanPRF ) ;
    double Tsa = hdr->pulses[hdr->num_azi-1].sat_tx_time - hdr->pulses[0].sat_tx_time ;
    printf("Collection duration         : %7.1f secs\n",Tsa) ;
    double bandwidth = hdr->chirp_gamma * hdr->pulse_length ;
    printf("Transmit bandwidth          : %7.3f MHz\n",bandwidth/1e6);
    double slantResolution = SIPC_c/(2 * bandwidth);
    printf("Slant Range resolution      : %7.1f m\n",slantResolution) ;
    double groundResolution = slantResolution / cos(geom1.grazingRad) ;
    printf("Ground Range resolution     : %7.1f m\n",groundResolution) ;
    double fcent = hdr->freq_centre;
    printf("Transmit Centre Frequency   : %7.3f GHz\n", fcent/1e9) ;
    double fmin = fcent - bandwidth/2 ;
    double fmax = fcent + bandwidth/2 ;
    printf("Transmit Chirp start/stop   : %7.3f - %7.3f MHz\n", fmin/1e6, fmax/1e6) ;
    double wavelength = SIPC_c / fcent ;
    printf("Centre wavelength           : %7.3f m \n", wavelength) ;
    SPVector centPls = hdr->pulses[(hdr->num_azi/2)+1].sat_ps_tx ;
//    double centPlsT  = hdr->pulses[(hdr->num_azi/2)+1].sat_tx_time ;
    SPVector centMin = hdr->pulses[(hdr->num_azi/2)-1].sat_ps_tx ;
//    double centMinT  = hdr->pulses[(hdr->num_azi/2)-1].sat_tx_time ;
    SPVector diff; VECT_SUB(centPls, centMin, diff) ;
//    double sdiff = VECT_MAG(diff) ;
    printf("Synthetic aperture centre   : %f,%f,%f metres\n",
           hdr->pulses[hdr->num_azi/2].sat_ps_tx.x, hdr->pulses[hdr->num_azi/2].sat_ps_tx.y, hdr->pulses[hdr->num_azi/2].sat_ps_tx.z) ;
//    double VelCentre = sdiff / (centPlsT - centMinT) ;
//    printf("Velocity at aperture Centre : %7.1f m/s\n", VelCentre) ;
    printf("Velocity at aperture Centre : %7.1f m/s\n", VECT_MAG( geom1.vel )) ;
    VECT_SUB(hdr->pulses[hdr->num_azi-1].sat_ps_tx, hdr->pulses[0].sat_ps_tx, diff);
    double SADistance = VECT_MAG(diff) ;
    printf("Synthetic aperture distance : %7.1f m\n", SADistance) ;
    double azResolution = fabs(geom1.range * wavelength / (2.0 * SADistance * cos(geom1.squintRad)));
    printf("Finest azimuth resolution   : %7.3f m\n",azResolution) ;
    double rSwath = hdr->clock_speed * SIPC_c / (2 * hdr->chirp_gamma);
    printf("Maximum range scene size    : %7.3f m \n",rSwath) ;
    if(hdr->antenna_width_az >1e-10){
        printf("Illuminated area (az x ra)  : %7.1f x %7.1fm\n", geom1.range * hdr->antenna_width_az, (geom1.range * hdr->antenna_width_el) / sin(geom1.grazingRad)) ;
    }
    
    return (status->status) ;
}
