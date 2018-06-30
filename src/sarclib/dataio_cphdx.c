/***************************************************************************
 *
 *       Module:    dataio_cphdx.c
 *      Program:    sarclib
 *   Created by:    Matt Nottingham on 08/09/2009.
 *                  and Darren Muff 24/07/2013
 *                  Copyright (c) 2017 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains functions for reading Compensated Phase History
 *      DATA (CPHD) formatted data in cphdx format
 *
 *   CLASSIFICATION        :  OFFICIAL
 *   Date of CLASSN        :  29/08/2016
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

#include "alloc.h"
#include "sarclib.h"
#include "sarclib_xml.h"
#include "sensitive.h"
#include "dataio_cphdx.h"

static int readCPHDXNarrowbandV0p3(CPHDHeader * hdr, SPStatus * status);
static int read_cphdx_xml   (CPHDHeader * hdr, SPStatus * status);
static int convert_cphdx_xml(CPHDHeader * hdr, SPStatus * status);

int readCPHDXHeaderV0p3(CPHDHeader * hdr, SPStatus * status) {
    
    char line[CPHD_MAX_LINE_LENGTH];
    char * match;
    
    CHECK_STATUS(status) ;
    
    hdr->version = 'x';
    
    do {
        fgets(line, CPHD_MAX_LINE_LENGTH, hdr->fp);
        match = strstr(line, " := ");
        if (match) {
            if (strncmp(line, "XML_DATA_SIZE", 13) == 0) {
                hdr->xml_data_size = atol(&line[13 + 4]);
            }
            if (strncmp(line, "XML_BYTE_OFFSET", 15) == 0) {
                hdr->xml_byte_offset = atol(&line[15 + 4]);
            }
            if (strncmp(line, "VB_DATA_SIZE", 12) == 0) {
                hdr->vb_data_size = atol(&line[12 + 4]);
            }
            if (strncmp(line, "VB_BYTE_OFFSET", 14) == 0) {
                hdr->vb_byte_offset = atol(&line[14 + 4]);
            }
            if (strncmp(line, "CPHD_DATA_SIZE", 14) == 0) {
                hdr->cphd_data_size = atol(&line[14 + 4]);
            }
            if (strncmp(line, "CPHD_BYTE_OFFSET", 16) == 0) {
                hdr->cphd_byte_offset = atol(&line[16 + 4]);
                hdr->start_wb_data = hdr->cphd_byte_offset;
            }
            
        }
    } while (line[0] != '\f' && line[1] != '\n');
    
    if (im_machine_type() == IM_LITTLE_ENDIAN) {
        hdr->byte_swap = TRUE;
    } else {
        hdr->byte_swap = FALSE;
    }
    hdr->polarisation = NULL;
    
    // Now read in the XML
    //
    read_cphdx_xml(hdr, status);
    CHECK_STATUS(status) ;
    
    // Now parse the XML into normal CPHD header
    //
    convert_cphdx_xml(hdr, status);
    CHECK_STATUS(status) ;
    
    // Now read in the NB data
    //
    readCPHDXNarrowbandV0p3(hdr, status);
    CHECK_STATUS(status) ;
    
    hdr->pulse_size = im_getsizeoftype(hdr->data_type) * hdr->nsamp;
    
    return (status->status);
}

int readCPHDXPulsesV0p3(CPHDHeader *hdr, int start_pulse, int stop_pulse, SPStatus * status){
    
    FILE *fp = hdr->fp;
    int i;
    int s;
    int p;
    
    pthread_mutex_lock(&hdr->file_lock);
    
    if(hdr->pulses == NULL){
        if(status->debug >= 10)fprintf(stderr, "Error : NarrowBand data does not exist: %s, line %d\n",__FILE__,__LINE__);
        status->status = CPHD_READ_ERROR ;
    }
    
    CHECK_STATUS(status) ;
    
    fseeko(fp, hdr->cphd_byte_offset + hdr->pulse_size * start_pulse, SEEK_SET);
    
    for(i = start_pulse; i < stop_pulse; i++) {
        hdr->pulses[i].data.v =  calloc(hdr->nsamp, im_getsizeoftype(hdr->data_type));
        if (hdr->pulses[i].data.v == NULL)
        {
            if(status->debug>=10)fprintf(stderr, "Failed to calloc data for pulse in %s, line %d\n",__FILE__,__LINE__);
            status->status = OUT_OF_MEMORY ;
            return (status->status);
        }
        switch (hdr->data_type) {
            case ITYPE_CMPL_INT8:
                // This only reads in bytes, so we don't need to swap
                //
                if (fread(hdr->pulses[i].data.cmpl_i8, sizeof(signed char), 2 * hdr->nsamp, fp) != 2 * hdr->nsamp) {
                    return 1;
                }
                for(s = 0; s < (int)((hdr->pulses[i].fx1 - hdr->pulses[i].fx0) / hdr->pulses[i].fx_step_size); s++) {
                    hdr->pulses[i].data.cmpl_i8[s].r = 0;
                    hdr->pulses[i].data.cmpl_i8[s].i = 0;
                }
                for(s = (int)((hdr->pulses[i].fx2 - hdr->pulses[i].fx0) / hdr->pulses[i].fx_step_size) + 1; s < hdr->nsamp; s++) {
                    hdr->pulses[i].data.cmpl_i8[s].r = 0;
                    hdr->pulses[i].data.cmpl_i8[s].i = 0;
                }
                
                break;
                
            case ITYPE_CMPL_INT16:
                // To be safe, lets check if there is any padding
                //
                if (sizeof(SPCmplxInt16) == 4) {
                    if (fread_byte_swap(hdr->pulses[i].data.cmpl_i16, sizeof(int16_t), 2 * hdr->nsamp, fp, hdr->byte_swap, status) != 2 * hdr->nsamp) {
                        return 1;
                    }
                } else {
                    for(p = 0; p < hdr->nsamp; p++) {
                        if (fread_byte_swap(&hdr->pulses[i].data.cmpl_i16[p].r, sizeof(int16_t), 1, fp, hdr->byte_swap, status) != 1) {
                            return 1;
                        }
                        if (fread_byte_swap(&hdr->pulses[i].data.cmpl_i16[p].i, sizeof(int16_t), 1, fp, hdr->byte_swap, status) != 1) {
                            return 1;
                        }
                    }
                }
                for(s = 0; s < (int)((hdr->pulses[i].fx1 - hdr->pulses[i].fx0) / hdr->pulses[i].fx_step_size); s++) {
                    hdr->pulses[i].data.cmpl_i16[s].r = 0;
                    hdr->pulses[i].data.cmpl_i16[s].i = 0;
                }
                for(s = (int)((hdr->pulses[i].fx2 - hdr->pulses[i].fx0) / hdr->pulses[i].fx_step_size); s < hdr->nsamp; s++) {
                    hdr->pulses[i].data.cmpl_i16[s].r = 0;
                    hdr->pulses[i].data.cmpl_i16[s].i = 0;
                }
                break;
                
            case ITYPE_CMPL_FLOAT:
                if (fread_byte_swap(hdr->pulses[i].data.cmpl_f, sizeof(float), 2 * hdr->nsamp, fp, hdr->byte_swap, status) != 2 * hdr->nsamp) {
                    return 1;
                }
                for(s = 0; s < (int)((hdr->pulses[i].fx1 - hdr->pulses[i].fx0) / hdr->pulses[i].fx_step_size); s++) {
                    hdr->pulses[i].data.cmpl_f[s].r = 0;
                    hdr->pulses[i].data.cmpl_f[s].i = 0;
                }
                for(s = (int)((hdr->pulses[i].fx2 - hdr->pulses[i].fx0) / hdr->pulses[i].fx_step_size); s < hdr->nsamp; s++) {
                    hdr->pulses[i].data.cmpl_f[s].r = 0;
                    hdr->pulses[i].data.cmpl_f[s].i = 0;
                }
                break;
                
            default:
                if(status->debug >= 10)printf("Unhandled type in WB sample read (%d) (%s:%d)\n", hdr->data_type,__FILE__,__LINE__);
                status->status = INVALID_TYPE ;
                return ( status->status );
                break;
        }
    }
    
    pthread_mutex_unlock(&hdr->file_lock);
    
    return(0);
}

int load_cphdx(SPImage * data, CPHDHeader * hdr, int start, int stop, SPStatus * status) {
    int i, j;
    
    CHECK_STATUS(status);
    
    if(hdr->pulses == NULL){
        if(status->debug >= 10)fprintf(stderr, "Error : NarrowBand data does not exist: %s, line %d\n",__FILE__,__LINE__);
        status->status = CPHD_READ_ERROR ;
    }
    
    if (readCPHDXPulsesV0p3(hdr, start, stop, status)) {
        status->status = BAD_FILE;
    }
    
    CHECK_STATUS(status);
    if (data->nx == 0 && data->ny == 0) {
        im_create(data, ITYPE_CMPL_FLOAT, hdr->nsamp, stop-start, 1.0, 1.0, status);
    }
    CHECK_STATUS(status);
    
    if (data->nx != hdr->nsamp) {
        status->status = INPUT_NX_MISMATCHED;
    }
    
    if (data->ny != (stop-start)) {
        status->status = INPUT_NY_MISMATCHED;
    }
    
    CHECK_STATUS(status);
    
    for(i = 0; i < stop-start; i++)
    {
        if (hdr->data_type == ITYPE_CMPL_INT8) {
            for(j = 0; j < hdr->nsamp; j++)
            {
                data->data.cmpl_f[j + data->nx * i].r = (float)hdr->pulses[i + start].data.cmpl_i8[j].r;
                data->data.cmpl_f[j + data->nx * i].i = (float)hdr->pulses[i + start].data.cmpl_i8[j].i;
            }
        }
        if (hdr->data_type == ITYPE_CMPL_INT16) {
            for(j = 0; j < hdr->nsamp; j++)
            {
                data->data.cmpl_f[j + data->nx * i].r = (float)hdr->pulses[i + start].data.cmpl_i16[j].r;
                data->data.cmpl_f[j + data->nx * i].i = (float)hdr->pulses[i + start].data.cmpl_i16[j].i;
            }
        }
        if (hdr->data_type == ITYPE_CMPL_FLOAT) {
            for(j = 0; j < hdr->nsamp; j++)
            {
                data->data.cmpl_f[j + data->nx * i].r = (float)hdr->pulses[i + start].data.cmpl_f[j].r;
                data->data.cmpl_f[j + data->nx * i].i = (float)hdr->pulses[i + start].data.cmpl_f[j].i;
            }
        }
        
    }
    
    return(status->status);
}

int writeCPHDXFile(CPHDHeader *hdr, SPStatus * status){
    CHECK_STATUS(status) ;
    if (hdr->fp == NULL) {
        printf("Error: Output filename not set in CPHD Data structure\n");
        printf("Error: in writeCPHDXFile()\n");
        status->status = CPHD_WRITE_ERROR ;
        return status->status;
    }
    writeCPHDXHeader(hdr, hdr->fp, status);
    CHECK_STATUS(status) ;
    writeCPHDXNarrowband(hdr, 0, hdr->fp, status);
    CHECK_STATUS(status) ;

    fseeko(hdr->fp, hdr->cphd_byte_offset, SEEK_SET);
    
    for (int i = 0; i < hdr->num_azi; i++) {
        
        if ( fwrite_byte_swap(&(hdr->pulses[i].data.cmpl_f[0]), sizeof(float), 2 * hdr->nsamp, hdr->fp, hdr->byte_swap, status) != 2 * hdr->nsamp) {
            if(status->debug >= 10)fprintf(stderr, "Write file for IQ data\n");
            status->status = CPHD_WRITE_ERROR ;
            return (status->status);
        }
    }
    return status->status ;
}

int writeCPHDXHeader(CPHDHeader * hdr, FILE * fp, SPStatus * status) {
    CHECK_STATUS(status) ;
    
    // first set up the XML size/offsets in the header
    //
    
    int VBPBytes =
    8       // TxTime
    + 24    // TxPos x,y,z
    + 8     // RxTime
    + 24    // RxPos
    + 24    // SRPpos
    + 8     // Fx0
    + 8     // Fx_SS
    + 8     // Fx1
    + 8     // Fx2
    ;
    hdr->vb_data_size = VBPBytes * hdr->num_azi ;
    hdr->cphd_data_size = 2 * sizeof(float) * hdr->num_azi * hdr->nsamp ;
    
    // set up a memory buffer to write the XML into
    //
    int xmlmetabufsize = 20 * 1024 * 1024 ;    // 20MBytes
    char *xmlmetabuf = (char *)malloc(sizeof(char) * xmlmetabufsize);
    char *xmlhdr = (char *)malloc(sizeof(char) * BUFFSIZE);
    char *ptr;
    ptr=xmlmetabuf;
    
    
    // XML Metadata
    //
    writeToBufText(ptr, "<CPHD xmlns=\"urn:CPHD:0.3\">\n") ;
    writeToBufText(ptr, "  <CollectionInfo>\n") ;
    writeToBufParm(ptr, "      <CollectorName>%s</CollectorName>\n",             hdr->sensor) ;
    writeToBufParm(ptr, "      <CoreName>%s</CoreName>\n",                       hdr->dataSetID);
    writeToBufParm(ptr, "      <RadarMode><ModeType>%s</ModeType></RadarMode>\n",hdr->mode);
    writeToBufParm(ptr, "      <Classification>%s</Classification>\n",           hdr->classification);
    writeToBufText(ptr, "  </CollectionInfo>\n");
    writeToBufText(ptr, "  <Data>\n") ;
    writeToBufText(ptr, "    <SampleType>RE32F_IM32F</SampleType>\n") ;            // Complex Float
    writeToBufParm(ptr, "    <NumCPHDChannels>%d</NumCPHDChannels>\n",           hdr->nchan);
    writeToBufParm(ptr, "    <NumBytesVBP>%d</NumBytesVBP>\n",                   VBPBytes);
    writeToBufText(ptr, "    <ArraySize>\n") ;
    writeToBufParm(ptr, "        <NumVectors>%d</NumVectors>\n",                 hdr->num_azi);
    writeToBufParm(ptr, "        <NumSamples>%d</NumSamples>\n",                 hdr->nsamp);
    writeToBufText(ptr, "    </ArraySize>") ;
    writeToBufText(ptr, "  </Data>\n") ;
    writeToBufText(ptr, "  <Global>\n") ;
    writeToBufText(ptr, "    <DomainType>FX</DomainType>\n") ;
    writeToBufText(ptr, "    <PhaseSGN>-1</PhaseSGN>\n");
    writeToBufParm(ptr, "    <CollectStart>%s</CollectStart>\n",hdr->dateTime);
    writeToBufParm(ptr, "    <CollectDuration>%f</CollectDuration>\n",hdr->pulses[hdr->num_azi-1].sat_rx_time - hdr->pulses[0].sat_rx_time);
    writeToBufParm(ptr, "    <TxTime1>%f</TxTime1>\n",hdr->pulses[0].sat_tx_time);
    writeToBufParm(ptr, "    <TxTime2>%f</TxTime2>\n",hdr->pulses[hdr->num_azi-1].sat_tx_time);
    
    collectionGeom geom ;
    SPVector NULLVect; VECT_CREATE(0, 0, 0, NULLVect);
    collectionGeometry(hdr, -1, NULLVect, &geom, status) ;
    writeToBufText(ptr, "    <ImageArea>\n");
    writeToBufText(ptr, "      <Corners>\n");
    writeToBufText(ptr, "        <ACP index=\"1\">\n");
    writeToBufParm(ptr, "          <Lat>%3.12f</Lat>\n", geom.imgTL.lat);
    writeToBufParm(ptr, "          <Lon>%3.12f</Lon>\n", geom.imgTL.lon);
    writeToBufParm(ptr, "          <HAE>%3.12f</HAE>\n", geom.imgTL.alt);
    writeToBufText(ptr, "        </ACP>\n");
    writeToBufText(ptr, "        <ACP index=\"2\">\n");
    writeToBufParm(ptr, "          <Lat>%3.12f</Lat>\n", geom.imgTR.lat);
    writeToBufParm(ptr, "          <Lon>%3.12f</Lon>\n", geom.imgTR.lon);
    writeToBufParm(ptr, "          <HAE>%3.12f</HAE>\n", geom.imgTR.alt);
    writeToBufText(ptr, "        </ACP>\n");
    writeToBufText(ptr, "        <ACP index=\"3\">\n");
    writeToBufParm(ptr, "          <Lat>%3.12f</Lat>\n", geom.imgBR.lat);
    writeToBufParm(ptr, "          <Lon>%3.12f</Lon>\n", geom.imgBR.lon);
    writeToBufParm(ptr, "          <HAE>%3.12f</HAE>\n", geom.imgBR.alt);
    writeToBufText(ptr, "        </ACP>\n");
    writeToBufText(ptr, "        <ACP index=\"4\">\n");
    writeToBufParm(ptr, "          <Lat>%3.12f</Lat>\n", geom.imgBL.lat);
    writeToBufParm(ptr, "          <Lon>%3.12f</Lon>\n", geom.imgBL.lon);
    writeToBufParm(ptr, "          <HAE>%3.12f</HAE>\n", geom.imgBL.alt);
    writeToBufText(ptr, "        </ACP>\n");
    writeToBufText(ptr, "      </Corners>\n");
    writeToBufText(ptr, "    </ImageArea>\n");
    writeToBufText(ptr, "  </Global>\n") ;
    
    writeToBufText(ptr, "  <Channel>\n") ;
    for (int chan=0; chan < hdr->nchan; ++chan) {
        writeToBufParm(ptr, "    <Parameters index=\"%d\">\n",chan) ;
        writeToBufParm(ptr, "      <SRP_Index>%d</SRP_Index>\n",0) ;
        writeToBufText(ptr, "      <NomTOARateSF>1</NomTOARateSF>\n") ;
        writeToBufParm(ptr, "      <FxCtrNom>%f</FxCtrNom>\n",hdr->freq_centre ) ;
        writeToBufParm(ptr, "      <BWSavedNom>%f</BWSavedNom>\n",hdr->chirp_gamma * hdr->pulse_length) ;
        writeToBufParm(ptr, "      <TOASavedNom>%22.21f</TOASavedNom>\n",hdr->TOASaved) ;
        writeToBufText(ptr, "    </Parameters>\n") ;
    }
    writeToBufText(ptr, "  </Channel>\n") ;
    writeToBufText(ptr, "  <SRP>\n") ;
    writeToBufText(ptr, "    <SRPType>STEPPED</SRPType>\n") ;       // We will always write SRP into each vector so 'STEPPED'
    writeToBufText(ptr, "    <NumSRPs>0</NumSRPs>\n") ;             // STEPPED SRP type so always 0
    writeToBufText(ptr, "  </SRP>\n") ;
    writeToBufText(ptr, "  <VectorParameters>\n") ;
    writeToBufParm(ptr, "    <TxTime>%d</TxTime>\n",(int)sizeof(double)) ;
    writeToBufParm(ptr, "    <TxPos>%d</TxPos>\n",(int)sizeof(double)*3) ;
    writeToBufParm(ptr, "    <RcvTime>%d</RcvTime>\n",(int)sizeof(double)) ;
    writeToBufParm(ptr, "    <RcvPos>%d</RcvPos>\n",(int)sizeof(double)*3) ;
    writeToBufParm(ptr, "    <SRPPos>%d</SRPPos>\n",(int)sizeof(double)*3) ;
    writeToBufText(ptr, "    <FxParameters>\n") ;
    writeToBufParm(ptr, "       <Fx0>%d</Fx0>\n",(int)sizeof(double)) ;
    writeToBufParm(ptr, "       <Fx_SS>%d</Fx_SS>\n",(int)sizeof(double)) ;
    writeToBufParm(ptr, "       <Fx1>%d</Fx1>\n",(int)sizeof(double)) ;
    writeToBufParm(ptr, "       <Fx2>%d</Fx2>\n",(int)sizeof(double)) ;
    writeToBufText(ptr, "    </FxParameters>\n") ;
    writeToBufText(ptr, "  </VectorParameters>\n") ;
    writeToBufText(ptr, "</CPHD>\n") ;
    writeToBufText(ptr, "\f\n") ;
    
    hdr->xml_data_size  = strlen(xmlmetabuf)-2 ;
    
    // File Header
    //
    hdr->xml_byte_offset = 512 ;
    do {
        ptr = xmlhdr ;
        hdr->xml_byte_offset--;
        hdr->vb_byte_offset   = hdr->xml_byte_offset + hdr->xml_data_size ;
        hdr->cphd_byte_offset = hdr->vb_byte_offset  + hdr->vb_data_size  ;
        writeToBufText(ptr,"CPHD/0.3\n");
        writeToBufParm(ptr,"CLASSIFICATION := %s\n",        hdr->classification   ) ;
        writeToBufParm(ptr,"XML_DATA_SIZE := %lld\n",       hdr->xml_data_size    ) ;
        writeToBufParm(ptr,"XML_BYTE_OFFSET := %lld\n",     hdr->xml_byte_offset  ) ;
        writeToBufParm(ptr,"VB_DATA_SIZE := %lld\n",        hdr->vb_data_size     ) ;
        writeToBufParm(ptr,"VB_BYTE_OFFSET := %lld\n",      hdr->vb_byte_offset   ) ;
        writeToBufParm(ptr,"CPHD_DATA_SIZE := %lld\n",      hdr->cphd_data_size   ) ;
        writeToBufParm(ptr,"CPHD_BYTE_OFFSET := %lld\n",    hdr->cphd_byte_offset ) ;
        writeToBufText(ptr,"\f\n");
        
    } while (strlen(xmlhdr) != hdr->xml_byte_offset && hdr->xml_byte_offset > 0);
    
    if(hdr->xml_byte_offset <= 0){
        printf("Failed to find xml_byte_offset\n");
        status->status = CPHD_WRITE_ERROR;
        return status->status;
    }
    
    fprintf(fp, "%s", xmlhdr) ;
    fprintf(fp, "%s", xmlmetabuf);
    
    return 0;
}

int writeCPHDXNarrowband(CPHDHeader * hdr, int channel, FILE * out_fp, SPStatus * status){
    
    int idx;
    size_t f;
    
    CHECK_STATUS(status) ;
    
    if (im_machine_type() == IM_LITTLE_ENDIAN) {
        hdr->byte_swap = TRUE;
    } else {
        hdr->byte_swap = FALSE;
    }
    
    
    if(channel != 0){
        printf("Error: CPHDX writer does not support multi-channel datasets\n");
        status->status = CPHD_WRITE_ERROR ;
        return status->status ;
    }
    fseeko(out_fp, hdr->vb_byte_offset, SEEK_SET);
    
    for(idx = 0; idx < hdr->num_azi; idx++) {
        
        // Tx Time
        //
        f = fwrite_byte_swap(&hdr->pulses[idx].sat_tx_time, sizeof(double), 1, out_fp,  hdr->byte_swap, status) ;
        CHECK_BYTES(f, 1, status);
        
        // Tx Pos
        //
        f = fwrite_byte_swap(&hdr->pulses[idx].sat_ps_tx.x, sizeof(double), 1, out_fp,  hdr->byte_swap, status) ;
        CHECK_BYTES(f, 1, status);
        f = fwrite_byte_swap(&hdr->pulses[idx].sat_ps_tx.y, sizeof(double), 1, out_fp,  hdr->byte_swap, status) ;
        CHECK_BYTES(f, 1, status);
        f = fwrite_byte_swap(&hdr->pulses[idx].sat_ps_tx.z, sizeof(double), 1, out_fp,  hdr->byte_swap, status) ;
        CHECK_BYTES(f, 1, status);
        
        // Rx Time
        //
        f = fwrite_byte_swap(&hdr->pulses[idx].sat_rx_time, sizeof(double), 1, out_fp,  hdr->byte_swap, status) ;
        CHECK_BYTES(f, 1, status);
        
        // Rx Pos
        //
        f = fwrite_byte_swap(&hdr->pulses[idx].sat_ps_rx.x, sizeof(double), 1, out_fp,  hdr->byte_swap, status) ;
        CHECK_BYTES(f, 1, status);
        f = fwrite_byte_swap(&hdr->pulses[idx].sat_ps_rx.y, sizeof(double), 1, out_fp,  hdr->byte_swap, status) ;
        CHECK_BYTES(f, 1, status);
        f = fwrite_byte_swap(&hdr->pulses[idx].sat_ps_rx.z, sizeof(double), 1, out_fp,  hdr->byte_swap, status) ;
        CHECK_BYTES(f, 1, status);
        
        // SRP Pos
        //
        f = fwrite_byte_swap(&hdr->pulses[idx].srp.x, sizeof(double), 1, out_fp,  hdr->byte_swap, status) ;
        CHECK_BYTES(f, 1, status);
        f = fwrite_byte_swap(&hdr->pulses[idx].srp.y, sizeof(double), 1, out_fp,  hdr->byte_swap, status) ;
        CHECK_BYTES(f, 1, status);
        f = fwrite_byte_swap(&hdr->pulses[idx].srp.z, sizeof(double), 1, out_fp,  hdr->byte_swap, status) ;
        CHECK_BYTES(f, 1, status);
        
        // Fx0
        //
        f = fwrite_byte_swap(&hdr->pulses[idx].fx0, sizeof(double), 1, out_fp,  hdr->byte_swap, status) ;
        CHECK_BYTES(f, 1, status);
        
        // Fx Step Size
        //
        f = fwrite_byte_swap(&hdr->pulses[idx].fx_step_size, sizeof(double), 1, out_fp,  hdr->byte_swap, status) ;
        CHECK_BYTES(f, 1, status);
        
        // Fx1
        //
        f = fwrite_byte_swap(&hdr->pulses[idx].fx1, sizeof(double), 1, out_fp,  hdr->byte_swap, status) ;
        CHECK_BYTES(f, 1, status);
        
        // Fx2
        //
        f = fwrite_byte_swap(&hdr->pulses[idx].fx2, sizeof(double), 1, out_fp,  hdr->byte_swap, status) ;
        CHECK_BYTES(f, 1, status);
        
    }
    
    return (status->status);
    
}

int writeCPHDXWideband(CPHDHeader * hdr, SPImage * data, FILE * out_fp, SPStatus * status){
    
    int i;
    
    CHECK_STATUS(status) ;
    
    if (data->image_type != ITYPE_CMPL_FLOAT) {
        if(status->debug >= 10)printf("Data needs to be complex float write_cphd_wb_vectors (%s:%d)\n", __FILE__, __LINE__);
        status->status = INVALID_TYPE;
        return (status->status);
    }
    
    fseeko(out_fp, hdr->cphd_byte_offset, SEEK_SET);
    
    for (i = 0; i < data->ny; i++) {
        if (fwrite_byte_swap(&(data->data.cmpl_f[i * data->nx]), sizeof(float), 2 * data->nx, out_fp, hdr->byte_swap, status)  != 2 * data->nx) {
            if(status->debug >= 10)fprintf(stderr, "Write file for IQ data\n");
            status->status = CPHD_WRITE_ERROR ;
            return (status->status);
        }
    }
    
    return (status->status);
}


// Depracated Functions
//
int read_cphdx_pulses(CPHDHeader * hdr, int start_pulse, int stop_pulse) {
    SPStatus status ;
    im_init_status(status, 0) ;
    return (readCPHDXPulsesV0p3(hdr, start_pulse, stop_pulse, &status)) ;
}

// Local helper functions
//
static void XMLCALL start(void *data, const char *el, const char **attr) {
    strcpy(Names[Depth], el);
    Depth++;
}

static void XMLCALL end(void *data, const char *el) {
    Names[Depth][0] = 0;
    Depth--;
}

static void XMLCALL charhdl(void * data, const char * s, int len) {
    int i;
    int count;
    struct SPXMLListEle * newele;
    struct SPXMLListEle * j;
    char * cp;
    
    if (!(len == 1 && s[0] == '\n')) {
        
        newele = calloc(1, sizeof(struct SPXMLListEle));
        
        newele->value = (char *)calloc(len+1, sizeof(char));
        count = 0;
        for(i = 0; i < Depth; i++) {
            count += strlen(Names[i])+2;
        }
        newele->name = (char *)calloc(count, sizeof(char));
        
        count = 0;
        for(i = 0; i < Depth; i++) {
            cp = strchr(Names[i], ':');
            newele->name[count] = '/';
            if (cp) {
                memcpy(&newele->name[count+1], &cp[1], strlen(cp));             // Takes into account of the null term of the string
            } else {
                memcpy(&newele->name[count+1], Names[i], strlen(Names[i])+1);   // Takes into account of the null term of the string
            }
            count = (int)strlen(newele->name);
        }
        
        memcpy(newele->value, s, len);
        
        if (cphdXMLlist) {
            for (j = cphdXMLlist; j->next != NULL; j = j->next);
            j->next = newele;
        } else {
            cphdXMLlist = newele;
        }
    }
}

static XML_Parser create_myXMLParser(void) {
    XML_Parser p;
    p = XML_ParserCreate(NULL);
    if (! p) {
        fprintf(stderr, "Couldn't allocate memory for parser\n");
        exit(-1);
    }
    
    XML_SetElementHandler(p, start, end);
    XML_SetCharacterDataHandler(p, charhdl);
    
    return p;
}

static int convert_cphdx_xml(CPHDHeader * hdr, SPStatus * status) {
    struct SPXMLListEle * j;
    
    CHECK_STATUS(status) ;
    
    unsigned short requiredCheck[NUMCPHDXRequiredItems] ;
    for(int i=0; i<NUMCPHDXRequiredItems; ++i) requiredCheck[i] = 0 ;
    
    hdr->geometry = strdup("bistatic");
    hdr->polarisation = NULL;
    hdr->fix_up_sensor = FALSE;
    
    hdr->pos_nb[0] = CPHD_NB_ChannelNumber;
    hdr->pos_nb[1] = CPHD_NB_VectorNumber;
    hdr->num_nb_items = 2;
    
    for (j = hdr->XMLlist; j != NULL; j = j->next) {
        if(status->debug >= 10 )printf("%s : %s\n", j->name, j->value);
        if (strcmp(j->name, "/CPHD/CollectionInfo/CollectorName") == 0) {
            hdr->sensor = strdup(j->value);
            requiredCheck[CollectorName] = 1;
        }
        if (strcmp(j->name, "/CPHD/CollectionInfo/Classification") == 0) {
            hdr->classification = (char *)malloc(sizeof(char) * strlen(j->value));
            hdr->classification = strdup(j->value);
            requiredCheck[Classification] = 1 ;
        }
        if (strcmp(j->name, "/CPHD/CollectionInfo/CoreName") == 0) {
            hdr->dataSetID = strdup(j->value);
            requiredCheck[CoreName] = 1;
        }
        if (strcmp(j->name, "/CPHD/CollectionInfo/RadarMode/ModeType") == 0) {
            if (strcasecmp(j->value, "spotlight") == 0) {
                hdr->mode = strdup("SPOTLIGHT");
            } else if (strncasecmp(j->value, "dynamic", 6) == 0) {
                hdr->mode = strdup("DYNAMIC");
            } else {
                hdr->mode = strdup(j->value);
            }
            requiredCheck[ModeType] = 1 ;
        }
        
        if (strcmp(j->name, "/CPHD/Data/SampleType") == 0) {
            hdr->data_type = ITYPE_UNKNOWN;
            
            if (strcmp(j->value, "RE08I_IM08I") == 0) {
                hdr->data_type = ITYPE_CMPL_INT8;
            }
            if (strcmp(j->value, "RE16I_IM16I") == 0) {
                hdr->data_type = ITYPE_CMPL_INT16;
            }
            if (strcmp(j->value, "RE32F_IM32F") == 0) {
                hdr->data_type = ITYPE_CMPL_FLOAT;
            }
            requiredCheck[SampleType] = 1 ;
        }
        if (strcmp(j->name, "/CPHD/Data/NumCPHDChannels") == 0) {
            hdr->nchan = (int)atol(j->value);
            requiredCheck[NumCPHDChannels] = 1 ;
        }
        if (strcmp(j->name, "/CPHD/Data/NumBytesVBP") == 0) {
            hdr->vb_data_size = (int)atol(j->value);
            requiredCheck[NumBytesVBP] = 1 ;
        }
        if (strcmp(j->name, "/CPHD/Data/ArraySize/NumVectors") == 0) {
            hdr->num_azi = (int)atol(j->value);
            requiredCheck[NumVectors] = 1 ;
        }
        if (strcmp(j->name, "/CPHD/Data/ArraySize/NumSamples") == 0) {
            hdr->nsamp = (int)atol(j->value);
            requiredCheck[NumSamples] = 1;
        }
        
        if (strcmp(j->name, "/CPHD/Global/DomainType") == 0) {
            if (strcmp(j->value,"FX")) {
                printf("ERROR : CPHDX file has domain type %s. Only type \"FX\" is supported\n",j->value) ;
                status->status = CPHD_READ_ERROR ;
                return status->status ;
            }
            requiredCheck[DomainType] = 1 ;
        }
        if (strcmp(j->name, "/CPHD/Global/PhaseSGN") == 0) {
            hdr->phaseSgn = (int)atol(j->value);
            requiredCheck[PhaseSGN] = 1 ;
        }
        if (strcmp(j->name, "/CPHD/Global/CollectStart") == 0) {
            hdr->dateTime = strdup(j->value);
            requiredCheck[CollectStart] = 1 ;
        }
        if (strcmp(j->name, "/CPHD/Global/CollectDuration") == 0) {
            requiredCheck[CollectDuration] = 1;
        }
        if (strcmp(j->name, "/CPHD/Global/TxTime1") == 0) {
            requiredCheck[TxTime1] = 1;
        }
        if (strcmp(j->name, "/CPHD/Global/TxTime2") == 0) {
            requiredCheck[TxTime2] = 1;
        }
        if (strcmp(j->name, "/CPHD/Global/ImageArea/Corners/ACP/Lat") == 0) {
            requiredCheck[ACP] = 1;
        }
        
        if (strcmp(j->name, "/CPHD/Channel/Parameters/SRP_Index") == 0) {
            requiredCheck[SRP_Index] = 1 ;
        }
        if (strcmp(j->name, "/CPHD/Channel/Parameters/NomTOARateSF") == 0) {
            requiredCheck[NomTOARateSF] = 1 ;
        }
        if (strcmp(j->name, "/CPHD/Channel/Parameters/FxCtrNom") == 0) {
            hdr->freq_centre = atof(j->value);
            requiredCheck[FxCtrNom] = 1 ;
        }
        if (strcmp(j->name, "/CPHD/Channel/Parameters/BWSavedNom") == 0) {
            requiredCheck[BWSavedNom] = 1 ;
        }
        if (strcmp(j->name, "/CPHD/Channel/Parameters/TOASavedNom") == 0) {
            hdr->TOASaved = atof(j->value);
            requiredCheck[TOASavedNom] = 1 ;
        }
        
        if (strcmp(j->name, "/CPHD/SRP/SRPType") == 0) {
            if (strcmp(j->value, "FIXEDPT") == 0) {
                hdr->fixedSRP = 1 ;
            } else {
                hdr->fixedSRP = 0 ;
            }
            requiredCheck[SRPType] = 1 ;
        }
        if (strcmp(j->name, "/CPHD/SRP/NumSRPs") == 0) {
            requiredCheck[NumSRPs] = 1 ;
        }
        
        if (strcmp(j->name, "/CPHD/VectorParameters/TxTime") == 0) {
            hdr->pos_nb[hdr->num_nb_items] = CPHD_NB_TxTime;
            hdr->num_nb_items++;
            if((int)atol(j->value) == 8) requiredCheck[TxTime] = 1 ;
        }
        if (strcmp(j->name, "/CPHD/VectorParameters/TxPos") == 0) {
            hdr->pos_nb[hdr->num_nb_items] = CPHD_NB_TxPos;
            hdr->num_nb_items++;
            if((int)atol(j->value) == 24) requiredCheck[TxPos] = 1 ;
        }
        if (strcmp(j->name, "/CPHD/VectorParameters/RcvTime") == 0) {
            hdr->pos_nb[hdr->num_nb_items] = CPHD_NB_RcvTime;
            hdr->num_nb_items++;
            if((int)atol(j->value) == 8) requiredCheck[RcvTime] = 1 ;
        }
        if (strcmp(j->name, "/CPHD/VectorParameters/RcvPos") == 0) {
            hdr->pos_nb[hdr->num_nb_items] = CPHD_NB_RcvPos;
            hdr->num_nb_items++;
            if((int)atol(j->value) == 24) requiredCheck[RcvPos] = 1 ;
        }
        if (strcmp(j->name, "/CPHD/VectorParameters/SRPTime") == 0) {
            hdr->pos_nb[hdr->num_nb_items] = CPHD_NB_SRPTime;
            hdr->num_nb_items++;
        }
        if (strcmp(j->name, "/CPHD/VectorParameters/SRPPos") == 0) {
            hdr->pos_nb[hdr->num_nb_items] = CPHD_NB_SRP;
            hdr->num_nb_items++;
            if((int)atol(j->value) == 24) requiredCheck[SRPPos] = 1 ;
        }
        if (strcmp(j->name, "/CPHD/VectorParameters/AmpSF") == 0) {
            hdr->pos_nb[hdr->num_nb_items] = CPHD_NB_AmpSF0;
            hdr->num_nb_items++;
        }
        if (strcmp(j->name, "/CPHD/VectorParameters/FxParameters/Fx0") == 0) {
            hdr->pos_nb[hdr->num_nb_items] = CPHD_NB_Fx0;
            hdr->num_nb_items++;
            if((int)atol(j->value) == 8) requiredCheck[Fx0] = 1 ;
        }
        if (strcmp(j->name, "/CPHD/VectorParameters/FxParameters/Fx_SS") == 0) {
            hdr->pos_nb[hdr->num_nb_items] = CPHD_NB_FxStepSize;
            hdr->num_nb_items++;
            if((int)atol(j->value) == 8) requiredCheck[Fx_SS] = 1 ;
        }
        if (strcmp(j->name, "/CPHD/VectorParameters/FxParameters/Fx1") == 0) {
            hdr->pos_nb[hdr->num_nb_items] = CPHD_NB_Fx1;
            hdr->num_nb_items++;
            if((int)atol(j->value) == 8) requiredCheck[Fx1] = 1 ;
        }
        
        if (strcmp(j->name, "/CPHD/VectorParameters/FxParameters/Fx2") == 0) {
            hdr->pos_nb[hdr->num_nb_items] = CPHD_NB_Fx2;
            hdr->num_nb_items++;
            if((int)atol(j->value) == 8) requiredCheck[Fx2] = 1;
        }
        
    }
    
    int truCnt = 0;
    for(int i=0; i<NUMCPHDXRequiredItems; ++i){
        if (requiredCheck[i] == 1) {
            truCnt++;
        }
    }
    if(truCnt != NUMCPHDXRequiredItems){
        printf("Error : Failed Intigrity Check CPHDX File does not contain all required elements\n");
        printf("Missing elements are:\n");
        for (int i=0; i<NUMCPHDXRequiredItems; ++i) {
            if(requiredCheck[i] == 0){
                printf("    %s\n",CPHDXReqdElements[i]);
            }
        }
//        status->status = CPHD_READ_ERROR ;
//        return status->status ;
    }
    
    // Check for clasified sensor to correct for frequency offsets required by security
    //
    if (hdr->sensor != NULL) {
        if ((strncmp(hdr->sensor, SENSEPARAM04, 3) == 0) && (hdr->freq_centre < SENSEPARAM03)) {
            hdr->freq_centre += SENSEPARAM03;
        }
    }
    
    hdr->pulse_length = hdr->TOASaved * 0.95; // Why-o-why isn't this saved in a cphd-x file?
    hdr->clock_speed = ((double) hdr->nsamp) / hdr->TOASaved;
    
    hdr->deskewed = TRUE;   // CPHD-X files are assumed to be deskewed
    hdr->nchan = 1;
    
    hdr->num_wb_items = 3;
    
    return (status->status);
}

static int read_cphdx_xml(CPHDHeader * hdr, SPStatus * status) {
    struct SPXMLListEle * j;
    CHECK_STATUS(status) ;
    
    XML_Parser p_hdr = create_myXMLParser();
    
    cphdXMLlist = NULL;
    
    fseeko(hdr->fp, hdr->xml_byte_offset, SEEK_SET);
    
    for (;;) {
        int done;
        int len;
        
        len = (int)fread(Buff, 1, BUFFSIZE, hdr->fp);
        if (ferror(hdr->fp)) {
            if(status->debug >= 10 )fprintf(stderr, "Read error\n");
            status->status = XML_READ_ERROR ;
            return (status->status) ;
        }
        done = feof(hdr->fp);
        
        if (XML_Parse(p_hdr, Buff, len, done) == XML_STATUS_ERROR) {
            done = 1;
        }
        
        if (done)
            break;
    }
    
    if(status->debug >= 10 ){
        printf("XML=\n");
        for (j = cphdXMLlist; j != NULL; j = j->next) {
            printf("%s : %s\n", j->name, j->value);
        }
    }
    
    XML_ParserFree(p_hdr);
    
    hdr->XMLlist = cphdXMLlist;
    
    return (status->status);
}

static int readCPHDXNarrowbandV0p3(CPHDHeader * hdr, SPStatus * status){
    
    int p;
    int i;
    
    CHECK_STATUS(status);
    
    hdr->pulses = (CPHDPulse *)calloc(hdr->num_azi, sizeof(CPHDPulse));
    fseeko(hdr->fp, hdr->vb_byte_offset, SEEK_SET);
    
    for(p = 0; p < hdr->num_azi; p++) {
        hdr->pulses[p].pulse_number = p;
        hdr->pulses[p].amp_sf0 = 1.0;
        
        for(i = 0; i < hdr->num_nb_items; i++) {
            switch (hdr->pos_nb[i]) {
                case CPHD_NB_TxTime:
                    fread_byte_swap(&hdr->pulses[p].sat_tx_time, sizeof(double), 1, hdr->fp, hdr->byte_swap, status);
                    break;
                    
                case CPHD_NB_TxPos:
                    fread_byte_swap(&hdr->pulses[p].sat_ps_tx, sizeof(double), 3, hdr->fp, hdr->byte_swap, status);
                    break;
                    
                case CPHD_NB_RcvTime:
                    fread_byte_swap(&hdr->pulses[p].sat_rx_time, sizeof(double), 1, hdr->fp, hdr->byte_swap, status);
                    break;
                    
                case CPHD_NB_RcvPos:
                    fread_byte_swap(&hdr->pulses[p].sat_ps_rx, sizeof(double), 3, hdr->fp, hdr->byte_swap, status);
                    break;
                    
                case CPHD_NB_SRPTime:
                    fread_byte_swap(&hdr->pulses[p].srp_time, sizeof(double), 1, hdr->fp, hdr->byte_swap, status);
                    break;
                    
                case CPHD_NB_SRP:
                    fread_byte_swap(&hdr->pulses[p].srp, sizeof(double), 3, hdr->fp, hdr->byte_swap, status);
                    break;
                    
                case CPHD_NB_AmpSF0:
                    fread_byte_swap(&hdr->pulses[p].amp_sf0, sizeof(double), 1, hdr->fp, hdr->byte_swap, status);
                    break;
                    
                case CPHD_NB_Fx0:
                    fread_byte_swap(&hdr->pulses[p].fx0, sizeof(double), 1, hdr->fp, hdr->byte_swap, status);
                    break;
                    
                case CPHD_NB_FxStepSize:
                    fread_byte_swap(&hdr->pulses[p].fx_step_size, sizeof(double), 1, hdr->fp, hdr->byte_swap, status);
                    break;
                    
                case CPHD_NB_Fx1:
                    fread_byte_swap(&hdr->pulses[p].fx1, sizeof(double), 1, hdr->fp, hdr->byte_swap, status);
                    break;
                    
                case CPHD_NB_Fx2:
                    fread_byte_swap(&hdr->pulses[p].fx2, sizeof(double), 1, hdr->fp, hdr->byte_swap, status);
                    break;
                    
                default:
                    continue;
            }
        }
        
        // Add frequency offset for sensor. Can be classified
        // so its included in sensitive.h
        //
        hdr->pulses[p].fx0 += SENSEPARAM03;
        hdr->pulses[p].fx1 += SENSEPARAM03;
        hdr->pulses[p].fx2 += SENSEPARAM03;
        
    }
    
    hdr->chirp_gamma = hdr->pulses[0].fx_step_size * hdr->clock_speed;
    
    return ( status->status );
}
