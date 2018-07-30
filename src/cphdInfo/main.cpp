/***************************************************************************
 *
 *       Module:    main.cpp
 *      Program:    readCPHD
 *   Created by:    Darren Muff on 8th Dec 12016
 *                  Copyright (c) 12013 [dstl]. All rights reserved.
 *
 *   Description:
 *      Main routine for programme to read a CPHD file and print out useful
 *      information about it
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  8/12/12016
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

#include <iostream>
#include <sarclib/sarclib.h>
#include "colourCodes.h"
#include "readCPHD_Version.h"

void banner (){
    printf(" \n");
    printf(DARK GREEN "              cphdInfo - A Programme to read CPHD Files\n" NORMAL);
    printf(BLUE "                         Version :" RED" %s \n", SHORT_VERSION);
    printf(BLUE "                          Commit :" RED" %s \n", COMMIT);
    printf(BLUE "                Revision: " RED"%s, %s \n",REVISION, VERSION_DATE);
    printf(BLUE "           Copyright (c) 2016 " WHITE"[" BLUE"Dstl" WHITE"]" BLUE". All rights reserved.\n" RESETCOLOR);
    printf(" \n");
    
    return ;
}

using namespace std;

int main(int argc, const char * argv[]) {
    
    banner() ;
    
    int start=0, num=1;
    
    SPStatus status ;
    im_init_status(status, 0) ;
    im_init_lib(&status, (char *)"cphdInfo", argc, (char **)argv);
    
    char *fname, *fname_cstr;
    fname_cstr = (char *)sp_malloc(sizeof(char) *255) ;
    sprintf(fname_cstr,"cphdFile.cph") ;
    fname = input_string("Enter cphd filename", "Fname",
                         "A CPHD filename to examine", fname_cstr) ;
    free(fname_cstr);
    
    bool opNB = false;
    opNB = input_yesno("Do you want to see the narrowband data in each pulse", "showNBData",
                       "If yes then the Narrowband data for each pulse in the CPHD file will be output", opNB) ;
    
    CPHDHeader hdr ;
    readCPHDHeader(fname, &hdr, &status);
    
    bool opWB = false;
    opWB = input_yesno("Do you want to save the wideband data to a file?", "showWBData",
                       "if 'yes' then you will be asked for a filename to store the wideband data to. \n \
                       The data is stored in a .dat file and covers only the pulese that you have specified", opWB);
    
    if(opNB || opWB){
        do {
            start=hdr.num_azi/2;
            start = input_int("Narrowband start pulse (default: NPulses/2)", "startPulse",
                              "The first pulse in the CPHD file to show the narrowband data for", start);
            num=hdr.num_azi-start ;
            num   = input_int("Narrowband number of pulses", "NPulses",
                              "The number of pulses to show the narrowband data for", num) ;
        } while ((start > hdr.num_azi-2) && (start+num > hdr.num_azi));
        
    }
    
    char *ofname=NULL, *prompt=NULL;
    if(opWB){
        FILE *fp;
        bool goodFile = false;
        prompt = (char *)sp_malloc(sizeof(char) * 255);
        sprintf(prompt,"/tmp/wbdata_%05d-%05d.dat", start, start+num);
        do{
            ofname = input_string("Name of file to save WB data to", "WBFilename",
                                  "Name of file to save WB data to. If it is not writable you will be asked again.", prompt);
            fp = fopen(ofname,"w");
            if (fp == NULL) {
                printf("File %s is not avaialble for writing\n",ofname);
                goodFile = false;
            }else{
                goodFile = true;
            }
            fclose(fp) ;
        }while (goodFile == false) ;
        free(prompt);
    }
    
    collectionGeom geom;
    SPVector nullgrp;
    VECT_CREATE(0, 0, 0, nullgrp) ;
    collectionGeometry(&hdr, hdr.num_azi/2, nullgrp, &geom, &status);
    printCPHDCollectionInfo(&hdr, &status) ;
    
    if(opNB){
        printf("\n");
        printf("                   Narrowband Data [%5d - %5d]           \n", start,start+num);
        printf("--------------------------------------------------------------------------------\n");
        for(int i=start; i<start+num ; ++i){
            printf("[%d]  Tx Time: %f, TxPos: %f,%f,%f\n",i,
                   hdr.pulses[i].sat_tx_time, hdr.pulses[i].sat_ps_tx.x, hdr.pulses[i].sat_ps_tx.y,hdr.pulses[i].sat_ps_tx.z);
            printf("[%d]  Rx Time: %f, RxPos: %f,%f,%f\n",i,
                   hdr.pulses[i].sat_rx_time, hdr.pulses[i].sat_ps_rx.x, hdr.pulses[i].sat_ps_rx.y,hdr.pulses[i].sat_ps_rx.z);
            printf("[%d]     fx0 : %f, fx_step: %f, Amp0: %f\n",i, hdr.pulses[i].fx0, hdr.pulses[i].fx_step_size, hdr.pulses[i].amp_sf0);
        }
    }
    
    if(opWB){
        printf("Saving Wideband data to %s ... \n",ofname);
        readCPHDPulses(&hdr, start, start+num, &status) ;
        SPImage pulse ;
        SPImage wbData ;
        im_create(&pulse,  ITYPE_CMPL_FLOAT, hdr.nsamp, 1, 1.0, 1.0, &status);
        im_create(&wbData, ITYPE_CMPL_FLOAT, hdr.nsamp, num, 1.0, 1.0, &status);
        
        for(int i=start; i<start+num; ++i){
            switch (hdr.data_type) {
                case ITYPE_CMPL_INT8:
                    for(int j=0; j<hdr.nsamp; ++j){
                        pulse.data.cmpl_i8[j].r = hdr.pulses[i].data.cmpl_i8[j].r ;
                        pulse.data.cmpl_i8[j].i = hdr.pulses[i].data.cmpl_i8[j].i ;
                    }
                    break;
                    
                case ITYPE_CMPL_INT16:
                    for(int j=0; j<hdr.nsamp; ++j){
                        pulse.data.cmpl_i16[j].r = hdr.pulses[i].data.cmpl_i16[j].r ;
                        pulse.data.cmpl_i16[j].i = hdr.pulses[i].data.cmpl_i16[j].i ;
                    }
                    break;
                    
                case ITYPE_CMPL_FLOAT:
                    for(int j=0; j<hdr.nsamp; ++j){
                        pulse.data.cmpl_f[j].r = hdr.pulses[i].data.cmpl_f[j].r ;
                        pulse.data.cmpl_f[j].i = hdr.pulses[i].data.cmpl_f[j].i ;
                    }
                    break;
                    
                default:
                    printf("Unknown data type : %d\n",hdr.data_type);
                    exit(0) ;
                    break;
            }
            im_insert(&pulse, 0, i-start, &wbData, &status);
        }
        im_fftw(&wbData, FFT_2D, &status) ;
        im_circshift(&wbData, wbData.nx/2 , wbData.ny/2, &status) ;
        im_save(&wbData, ofname, &status);
        im_destroy(&pulse, &status) ;
        im_destroy(&wbData, &status) ;
    }
    im_close_lib(&status);
    
    return 0;
}

