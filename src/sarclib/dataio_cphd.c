/***************************************************************************
 *
 *       Module:    dataio_cphd.c
 *      Program:    sarclib
 *   Created by:    Matt Nottingham on 25/08/2009.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains functions for reading Compensated Phase History
 *      DATA (CPHD) formatted data.
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

// private function prototypes
//
int read_cphd3_header(CPHDHeader * hdr);
int read_cphd3_pulses(CPHDHeader * hdr, int start_pulse, int stop_pulse);
int load_cphd3(SPImage * data, CPHDHeader * hdr, int start_pulse, int stop_pulse, SPStatus * status);
int read_cphdx_header(CPHDHeader *  hdr);
int read_cphdx_pulses(CPHDHeader * hdr, int start_pulse, int stop_pulse);
int load_cphdx(SPImage * data, CPHDHeader * hdr, int start_pulse, int stop_pulse, SPStatus * status);

// Read in the header of a CPHD3 or CPHDX file. This also reads the narrowband
// data for each of the pulses and loads the information into the CPHDHeader structure
//
int read_cphd_header(const char * fname, CPHDHeader *  hdr) {
    SPStatus status ;
    im_init_status(status, 0) ;
    return readCPHDHeader(fname, hdr, &status) ;
}

int readCPHDHeader(const char *fname, CPHDHeader * hdr, SPStatus * status){
    CHECK_STATUS(status) ;
    char line[CPHD_MAX_LINE_LENGTH];
    
    pthread_mutex_init(&hdr->file_lock, NULL);
    pthread_mutex_lock(&hdr->file_lock);
    
    hdr->fp = fopen(fname, "r");
    
    if (hdr->fp == NULL) {
        if(status->debug>=10){
            fprintf(stderr, "Failed to open file %s\n", fname);
            perror("Error opening file");
        }
        status->status = BAD_FILE ;
        return (status->status) ;
    }
    
    hdr->sensor = NULL;
    hdr->dateTime = NULL;
    hdr->dataSetID = NULL;
    hdr->mode = NULL;
    hdr->geometry = NULL;
    hdr->classification = NULL ;
    
    fgets(line, CPHD_MAX_LINE_LENGTH, hdr->fp);
    rewind(hdr->fp);
    if (strncmp(line, "CPHD/0.3", 8) == 0) {
        readCPHDXHeaderV0p3(hdr, status) ;
    }else if(strncmp(line, "CPHD/1.0", 8) == 0) {
        printf("CPHDX 1.0 not supported\n");
        status->status = CPHD_READ_ERROR ;
        return status->status;
    }else {
        readCPHD3Header(hdr, status) ;
    }
    
    pthread_mutex_unlock(&hdr->file_lock);
    
    return(0);
}

int read_cphd_pulses(CPHDHeader * hdr, int start_pulse, int stop_pulse){
    SPStatus status ;
    im_init_status(status, 0);
    return (readCPHDPulses(hdr, start_pulse, stop_pulse, &status)) ;
}
int readCPHDPulses(CPHDHeader * hdr, int start_pulse, int stop_pulse, SPStatus *status){
    int ret;
    CHECK_STATUS(status) ;
    
    if (hdr->version == 'x') {
        ret = readCPHDXPulsesV0p3(hdr, start_pulse, stop_pulse, status) ;
    } else {
        ret = readCPHD3Pulses(hdr, start_pulse, stop_pulse, status) ;
    }
    
    return (status->status);
}

int load_cphd(SPImage * data, CPHDHeader * hdr, int start_pulse, int stop_pulse, SPStatus * status){
    int ret;
    
    if (hdr->version == 'x') {
        ret = load_cphdx(data, hdr, start_pulse, stop_pulse, status);
    } else {
        ret = load_cphd3(data, hdr, start_pulse, stop_pulse, status);
    }
    
    return ret;
}

int destroy_cphd_header(CPHDHeader * hdr, SPStatus * status){
    CHECK_STATUS(status);
    
    if (hdr->pulses != NULL) {
        destroy_cphd_pulses(hdr, 0, hdr->num_azi, status);
        pthread_mutex_lock(&hdr->file_lock);
        free(hdr->pulses);
        hdr->pulses = NULL;
        pthread_mutex_unlock(&hdr->file_lock);
        pthread_mutex_destroy(&hdr->file_lock);
    }
    
    if (hdr->sensor) {
        free(hdr->sensor);
    }
    
    free(hdr->dateTime);
    free(hdr->dataSetID);
    free(hdr->mode);
    free(hdr->geometry);
    
    fclose(hdr->fp);
    
    return (0);
}

int destroy_cphd_pulses(CPHDHeader * hdr, int start_pulse, int last_pulse, SPStatus * status){
    int i;
    
    CHECK_STATUS(status);
    
    pthread_mutex_lock(&hdr->file_lock);
    
    if (hdr->pulses != NULL) {
        for(i = start_pulse; i < last_pulse; i++) {
            if (hdr->pulses[i].data.v != NULL) {
                free(hdr->pulses[i].data.v);
                hdr->pulses[i].data.v = NULL;
            }
        }
    }
    
    pthread_mutex_unlock(&hdr->file_lock);
    
    return (0);
}

int writeCPHDFile(CPHDHeader *hdr, SPStatus * status){
    
    CHECK_STATUS(status) ;

    pthread_mutex_init(&hdr->file_lock, NULL);
    pthread_mutex_lock(&hdr->file_lock);

    if (hdr->version == 'x' ) {
        writeCPHDXFile(hdr, status) ;
        CHECK_STATUS(status) ;

    }else if(hdr->version == '3') {
        writeCPHD3File(hdr, status) ;
        CHECK_STATUS(status) ;
    }else {
        printf("Error: Version flag incorrectly set in CPHD data structure\n");
        status->status = CPHD_WRITE_ERROR ;
        return status->status ;
    }
    
    pthread_mutex_unlock(&hdr->file_lock);
    
    return status->status ;
}
