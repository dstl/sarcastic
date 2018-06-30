/***************************************************************************
 *
 *       Module:    dataio_cphd3.c
 *      Program:    sarclib
 *   Created by:    Darren Muff, Matt Nottingham on 08/09/2009.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains functions for reading Compensated Phase History
 *      DATA (CPHD) formatted data in cphd3 format
 *
 *   CLASSIFICATION        :  UK OFFICIAL
 *   Date of CLASSN        :  30th May 2016
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
#include "sensitive.h"

static int readCPHD3Narrowband(CPHDHeader * hdr, SPStatus *status) ;
static void cphd3_get_nb_order(CPHDHeader * hdr);
static void cphd3_get_wb_order(CPHDHeader * hdr);

// These need to match CPHD_NB_items & CPHD_WB_items respectively
//
char * nb_names[] = {"ChannelNumber", "VectorNumber", "SRP", "SRPTime", "TxPos", "RcvPos", "TxTime", "RcvTime", "Fx0", "FxStepSize", "Fx1", "Fx2", "AmpSF0"};
char * wb_names[] = {"ChannelNumber", "VectorNumber", "SampleData"};

int read_cphd3_header(CPHDHeader *  hdr) {
    SPStatus status ;
    im_init_status(status, 0) ;
    return (readCPHD3Header(hdr, &status));
}

int readCPHD3Header(CPHDHeader *  hdr, SPStatus *status) {
    char line[CPHD_MAX_LINE_LENGTH];
    char name[CPHD_MAX_LINE_LENGTH];
    char value1[CPHD_MAX_LINE_LENGTH];
    char value2[CPHD_MAX_LINE_LENGTH];
    unsigned int magic;
    
    CHECK_STATUS(status) ;
    
    hdr->sensor = NULL;
    hdr->dateTime = NULL;
    hdr->dataSetID = NULL;
    hdr->mode = NULL;
    hdr->geometry = NULL;
    
    hdr->version = '3';
    
    do {
        fgets(line, CPHD_MAX_LINE_LENGTH, hdr->fp);
        value1[0] = 0;
        value2[0] = 0;
        sscanf(line, "%s%s%s", name, value1, value2);
        if (strcmp(name, "Nvectors") == 0) {
            hdr->num_azi = atoi(value1);
        }
        
        if (strcmp(name, "Nsamples") == 0) {
            hdr->nsamp = atoi(value1);
        }
        
        if (strcmp(name, "XmitPulseDuration") == 0) {
            hdr->pulse_length = atof(value1);
        }
        
        if (strcmp(name, "Sensor") == 0) {
            hdr->sensor = strdup(value1);
        }
        
        if (strcmp(name, "NominalADRate") == 0) {
            hdr->clock_speed = atof(value1);
        }
        
        if (strcmp(name, "NominalCenterFreq") == 0) {
            hdr->freq_centre = atof(value1);
        }
        
        if (strcmp(name, "NominalChirpRate") == 0) {
            hdr->chirp_gamma = atof(value1);
        }
        
        if (strcmp(name, "Nchannels") == 0) {
            hdr->nchan = atoi(value1);
            if (hdr->nchan != 1) {
                if(status->debug >= 10)printf("Can only handle 1 channel data, not %d\n", hdr->nchan);
                status->status = CPHD_CHANNEL_ERROR ;
                return (status->status) ;
            }
            hdr->polarisation = (char **)sp_malloc(sizeof(char **) * hdr->nchan) ;
            for(int i=0; i<hdr->nchan; ++i){
                hdr->polarisation[i] = (char *)sp_malloc(sizeof(char) * 3);
            }
        }
        if (strcmp(name, "ChanDesc") == 0) {
            int nc = atoi(value1);
            strcpy(hdr->polarisation[nc], value2) ;
        }
        if (strcmp(name, "DateTime") == 0) {
            hdr->dateTime = strdup(value1);
        }
        
        if (strcmp(name, "DataSetID") == 0) {
            hdr->dataSetID = strdup(value1);
        }
        
        if (strcasecmp(name, "mode") == 0) {
            hdr->mode = strdup(value1);
        }
        
        if (strcasecmp(name, "geometry") == 0) {
            hdr->geometry = strdup(value1);
        }
        
        if (strcasecmp(name, "FixedSRP") == 0) {
            hdr->fixedSRP = (strcmp(value1, "Yes") == 0);
        }
        
        if (strcmp(name, "PhaseSgn") == 0) {
            hdr->phaseSgn = atoi(value1);
        }
        
        if (strcmp(name, "TOASaved") == 0) {
            hdr->TOASaved = atof(value1);
        }
        
        if (strcmp(name, "DeskewApplied") == 0) {
            hdr->deskewed = (strcmp(value1, "Yes") == 0);
        }
        
        if (strcmp(name, "PHDataType") == 0) {
            hdr->data_type = ITYPE_UNKNOWN;
            
            if (strcmp(value1, "cmplxn") == 0) {
                hdr->data_type = ITYPE_CMPL_NIBBLE;
                if(status->debug >= 10)printf("Unsupported type NIBBLE in PHDataType\n");
                status->status = CPHD_DATATYPE_ERROR ;
                return (status->status) ;
            }
            if (strcmp(value1, "cmplxb") == 0) {
                hdr->data_type = ITYPE_CMPL_INT8;
            }
            if (strcmp(value1, "cmplxs") == 0) {
                hdr->data_type = ITYPE_CMPL_INT16;
            }
            if (strcmp(value1, "cmplxf") == 0) {
                hdr->data_type = ITYPE_CMPL_FLOAT;
            }
            if (hdr->data_type == ITYPE_UNKNOWN) {
                if(status->debug >= 10)printf("Unknown data type! (%s)\n", value1);
                status->status = CPHD_DATATYPE_ERROR ;
                return (status->status) ;
            }
        }
        
        if (strcmp(name, "Interleaved") == 0) {
            hdr->interleaved = (strcmp(value1, "Yes") == 0);
            if (hdr->interleaved) {
                if(status->debug >= 10)printf("Cannot currently handle interleaved data!\n");
                status->status = CPHD_DATATYPE_ERROR ;
                return (status->status) ;
            }
        }
        
        if (strcmp(name, "AntennaPatternWidth") == 0) {
            if (strcmp(value1, "HalfPowerWidthAz") == 0) {
                hdr->antenna_width_az = atof(value2);
            }
            if (strcmp(value1, "HalfPowerWidthEl") == 0) {
                hdr->antenna_width_el = atof(value2);
            }
        }
        
        if (strcmp(name, "BegNBVector") == 0) {
            cphd3_get_nb_order(hdr);
        }
        
        if (strcmp(name, "BegWBVector") == 0) {
            cphd3_get_wb_order(hdr);
        }
        
        if (strcmp(name, "Classification") == 0) {
            hdr->classification = strdup(value1);
        }
        
    } while (strcmp(name, "EndPreamble") != 0);
    
    hdr->fix_up_sensor = FALSE;
    
    if (hdr->sensor) {
        if (strcmp(hdr->sensor, SENSEPARAM01) == 0 && hdr->freq_centre < SENSEPARAM02) {
            hdr->freq_centre += SENSEPARAM03;
            hdr->phaseSgn *= -1;
            hdr->fix_up_sensor = TRUE;
        }
    }
    
    if (fread(&magic, sizeof(unsigned int), 1, hdr->fp) != 1) {
        return 1;
    }
    if (magic == CPHD_MAGIC_NUMBER) {
        hdr->byte_swap = 0;
    } else if (magic == CPHD_MAGIC_NUMBER_BYTE_SWAPPED) {
        hdr->byte_swap = 1;
    } else {
        if(status->debug >= 10)fprintf(stderr, "Failed to read in correct magic number (read in %x)!!\n", magic);
        status->status = CPHD_MAGICNUM_ERROR ;
        return (status->status) ;
    }
    
    hdr->start_nb_data = ftello(hdr->fp);
    hdr->pulse_size = im_getsizeoftype(hdr->data_type) * hdr->nsamp + 2 * sizeof(int);
    hdr->pulses = NULL;
    
    if (readCPHD3Narrowband(hdr, status)) {
        if (status->debug >= 10)fprintf(stderr, "ERROR: readCPHD3Narrowband returned bad status value %d\n",status->status);
        return (status->status) ;
    }
    
    pthread_mutex_unlock(&hdr->file_lock);
    
    return(status->status);
}

static void cphd3_get_nb_order(CPHDHeader * hdr) {
    char line[CPHD_MAX_LINE_LENGTH];
    char name[CPHD_MAX_LINE_LENGTH];
    
    int finished = FALSE;
    int pos = 0;
    int parsed;
    
    hdr->nb_vector_size = 0;
    
    do {
        fgets(line, CPHD_MAX_LINE_LENGTH, hdr->fp);
        sscanf(line, "%s", name);
        parsed = FALSE;
        
        if (strcmp(name, "ChannelNumber") == 0) {
            hdr->pos_nb[pos] = CPHD_NB_ChannelNumber;
            pos++;
            parsed = TRUE;
            hdr->nb_vector_size += sizeof(int);
        }
        
        if (strcmp(name, "VectorNumber") == 0) {
            hdr->pos_nb[pos] = CPHD_NB_VectorNumber;
            pos++;
            parsed = TRUE;
            hdr->nb_vector_size += sizeof(int);
        }
        
        if (strcmp(name, "SRP") == 0) {
            hdr->pos_nb[pos] = CPHD_NB_SRP;
            pos++;
            parsed = TRUE;
            hdr->nb_vector_size += sizeof(double) * 3;
        }
        
        if (strcmp(name, "TxPos") == 0) {
            hdr->pos_nb[pos] = CPHD_NB_TxPos;
            pos++;
            parsed = TRUE;
            hdr->nb_vector_size += sizeof(double) * 3;
        }
        
        if (strcmp(name, "RcvPos") == 0) {
            hdr->pos_nb[pos] = CPHD_NB_RcvPos;
            pos++;
            parsed = TRUE;
            hdr->nb_vector_size += sizeof(double) * 3;
        }
        
        if (strcmp(name, "TxTime") == 0) {
            hdr->pos_nb[pos] = CPHD_NB_TxTime;
            pos++;
            parsed = TRUE;
            hdr->nb_vector_size += sizeof(double);
        }
        
        if (strcmp(name, "RcvTime") == 0) {
            hdr->pos_nb[pos] = CPHD_NB_RcvTime;
            pos++;
            parsed = TRUE;
            hdr->nb_vector_size += sizeof(double);
        }
        
        if (strcmp(name, "Fx0") == 0) {
            hdr->pos_nb[pos] = CPHD_NB_Fx0;
            pos++;
            parsed = TRUE;
            hdr->nb_vector_size += sizeof(double);
        }
        
        if (strcmp(name, "FxStepSize") == 0) {
            hdr->pos_nb[pos] = CPHD_NB_FxStepSize;
            pos++;
            parsed = TRUE;
            hdr->nb_vector_size += sizeof(double);
        }
        
        if (strcmp(name, "AmpSF0") == 0) {
            hdr->pos_nb[pos] = CPHD_NB_AmpSF0;
            pos++;
            parsed = TRUE;
            hdr->nb_vector_size += sizeof(double);
        }
        
        if (strcmp(name, "EndNBVector") == 0) {
            finished = TRUE;
            parsed = TRUE;
            hdr->num_nb_items = pos;
        }
        
        if (!parsed) {
            printf("Unknown NB data item \"%s\"\n", name);
            pos++;
        }
        
        if (pos >= CPHD_MAX_NB_ITEMS) {
            printf("Too many items in the NBVector description!\n");
            exit(864548);
        }
        
    } while(!finished);
}

static void cphd3_get_wb_order(CPHDHeader * hdr) {
    char line[CPHD_MAX_LINE_LENGTH];
    char name[CPHD_MAX_LINE_LENGTH];
    
    int finished = FALSE;
    int pos = 0;
    int parsed;
    
    do {
        fgets(line, CPHD_MAX_LINE_LENGTH, hdr->fp);
        sscanf(line, "%s", name);
        parsed = FALSE;
        
        if (strcmp(name, "ChannelNumber") == 0) {
            hdr->pos_wb[pos] = CPHD_WB_ChannelNumber;
            pos++;
            parsed = TRUE;
        }
        
        if (strcmp(name, "VectorNumber") == 0) {
            hdr->pos_wb[pos] = CPHD_WB_VectorNumber;
            pos++;
            parsed = TRUE;
        }
        
        if (strcmp(name, "SampleData") == 0) {
            hdr->pos_wb[pos] = CPHD_WB_Samples;
            pos++;
            parsed = TRUE;
        }
        if (strcmp(name, "EndWBVector") == 0) {
            finished = TRUE;
            parsed = TRUE;
            hdr->num_wb_items = pos;
        }
        
        if (!parsed) {
            printf("Unknown WB data item \"%s\"\n", name);
            pos++;
        }
        
        if (pos >= CPHD_MAX_WB_ITEMS) {
            printf("Too many items in the WBVector description!\n");
            exit(864545);
        }
        
    } while(!finished);
}

static int readCPHD3Narrowband(CPHDHeader * hdr, SPStatus *status) {
    FILE * fp = hdr->fp;
    int i;
    int ch;
    int idx;
    
    CHECK_STATUS(status) ;
    
    hdr->pulses = (CPHDPulse *)calloc(hdr->num_azi, sizeof(CPHDPulse)) ;
    if (hdr->pulses == NULL) {
        if(status->debug >=10 )fprintf(stderr, "Failed to allocate enough mem for %d narrowband items\n", hdr->num_azi);
        status->status = OUT_OF_MEMORY ;
        return (status->status);
    }
    
    for(idx = 0; idx < hdr->num_azi; idx++) {
        
        // Each pulse contains narrowband data pertaining to the pulse and an array of wideband data
        // Set the wideband data to null for now until it is read in through a call to
        // readWBData()
        //
        hdr->pulses[idx].data.v = NULL ;
        
        // Amp_sf0 is an optional entry, so lets set it to 1 to make everyone who uses the data have
        // a sane value if it isn't set
        //
        hdr->pulses[idx].amp_sf0 = 1.0;
        
        for(i = 0; i < hdr->num_nb_items; i++) {
            switch (hdr->pos_nb[i]) {
                case   CPHD_NB_ChannelNumber:
                    if (fread_byte_swap(&ch, sizeof(int), 1, fp, hdr->byte_swap, status) != 1) {
                        return 1;
                    }
                    break;
                    
                case   CPHD_NB_VectorNumber:
                    if (fread_byte_swap(&hdr->pulses[idx].pulse_number, sizeof(int), 1, fp, hdr->byte_swap, status) != 1) {
                        return 1;
                    }
                    break;
                    
                case   CPHD_NB_SRP:
                    if (fread_byte_swap(&hdr->pulses[idx].srp.x, sizeof(double), 1, fp, hdr->byte_swap, status) != 1) {
                        return 1;
                    }
                    if (fread_byte_swap(&hdr->pulses[idx].srp.y, sizeof(double), 1, fp, hdr->byte_swap, status) != 1) {
                        return 1;
                    }
                    if (fread_byte_swap(&hdr->pulses[idx].srp.z, sizeof(double), 1, fp, hdr->byte_swap, status) != 1) {
                        return 1;
                    }
                    break;
                    
                case   CPHD_NB_TxPos:
                    if (fread_byte_swap(&hdr->pulses[idx].sat_ps_tx.x, sizeof(double), 1, fp, hdr->byte_swap, status) != 1) {
                        return 1;
                    }
                    if (fread_byte_swap(&hdr->pulses[idx].sat_ps_tx.y, sizeof(double), 1, fp, hdr->byte_swap, status) != 1) {
                        return 1;
                    }
                    if (fread_byte_swap(&hdr->pulses[idx].sat_ps_tx.z, sizeof(double), 1, fp, hdr->byte_swap, status) != 1) {
                        return 1;
                    }
                    break;
                    
                case   CPHD_NB_RcvPos:
                    if (fread_byte_swap(&hdr->pulses[idx].sat_ps_rx.x, sizeof(double), 1, fp, hdr->byte_swap, status) != 1) {
                        return 1;
                    }
                    if (fread_byte_swap(&hdr->pulses[idx].sat_ps_rx.y, sizeof(double), 1, fp, hdr->byte_swap, status) != 1) {
                        return 1;
                    }
                    if (fread_byte_swap(&hdr->pulses[idx].sat_ps_rx.z, sizeof(double), 1, fp, hdr->byte_swap, status) != 1) {
                        return 1;
                    }
                    break;
                    
                case   CPHD_NB_TxTime:
                    if (fread_byte_swap(&hdr->pulses[idx].sat_tx_time, sizeof(double), 1, fp, hdr->byte_swap, status) != 1) {
                        return 1;
                    }
                    break;
                    
                case   CPHD_NB_RcvTime:
                    if (fread_byte_swap(&hdr->pulses[idx].sat_rx_time, sizeof(double), 1, fp, hdr->byte_swap, status) != 1) {
                        return 1;
                    }
                    break;
                    
                case   CPHD_NB_Fx0:
                    if (fread_byte_swap(&hdr->pulses[idx].fx0, sizeof(double), 1, fp, hdr->byte_swap, status) != 1) {
                        return 1;
                    }
                    if (hdr->fix_up_sensor) {
                        hdr->pulses[idx].fx0 += SENSEPARAM03;
                    }
                    break;
                    
                case   CPHD_NB_FxStepSize:
                    if (fread_byte_swap(&hdr->pulses[idx].fx_step_size, sizeof(double), 1, fp, hdr->byte_swap, status) != 1) {
                        return 1;
                    }
                    break;
                    
                case   CPHD_NB_AmpSF0:
                    if (fread_byte_swap(&hdr->pulses[idx].amp_sf0, sizeof(double), 1, fp, hdr->byte_swap, status) != 1) {
                        return 1;
                    }
                    break;
            }
        }
    }
    
    hdr->grp = hdr->pulses[0].srp;
    
    hdr->start_wb_data = ftello(fp);
    
    pthread_mutex_unlock(&hdr->file_lock);
    
    return(status->status);
}

int read_cphd3_pulses(CPHDHeader * hdr, int start_pulse, int stop_pulse) {
    SPStatus status ;
    im_init_status(status, 0) ;
    return (readCPHD3Pulses(hdr, start_pulse, stop_pulse, &status)) ;
}

int readCPHD3Pulses(CPHDHeader * hdr, int start_pulse, int stop_pulse, SPStatus *status) {
    FILE *fp = hdr->fp;
    int ch, vect;
    int i;
    int n;
    int p;
    
    pthread_mutex_lock(&hdr->file_lock);
    
    CHECK_STATUS(status) ;
    
    if(hdr->pulses == NULL){
        if(status->debug >= 10)fprintf(stderr, "Error : NarrowBand data does not exist: %s, line %d\n",__FILE__,__LINE__);
        status->status = CPHD_READ_ERROR ;
        return (status->status) ;
    }
    
    fseeko(fp, hdr->start_wb_data + hdr->pulse_size * start_pulse, SEEK_SET);
    
    for(i = start_pulse; i < stop_pulse; i++) {
        for(n = 0; n < hdr->num_wb_items; n++) {
            switch (hdr->pos_wb[n]) {
                case CPHD_WB_ChannelNumber:
                    if (fread_byte_swap(&ch, sizeof(int), 1, fp, hdr->byte_swap, status) != 1) {
                        return 1;
                    }
                    break;
                    
                case CPHD_WB_VectorNumber:
                    if (fread_byte_swap(&vect, sizeof(int), 1, fp, hdr->byte_swap, status) != 1) {
                        return 1;
                    }
                    break;
                    
                case CPHD_WB_Samples:
                {
                    hdr->pulses[i].data.v =  calloc(hdr->nsamp, im_getsizeoftype(hdr->data_type));
                    if (hdr->pulses[i].data.v == NULL)
                    {
                        if(status->debug >= 10 )fprintf(stderr, "Failed to calloc data for pulse\n");
                        status->status = OUT_OF_MEMORY ;
                        return (status->status) ;
                    }
                    switch (hdr->data_type) {
                        case ITYPE_CMPL_INT8:
                            // This only reads in bytes, so we don't need to swap
                            //
                            if (fread(hdr->pulses[i].data.cmpl_i8, sizeof(signed char), 2 * hdr->nsamp, fp) != 2 * hdr->nsamp) {
                                if(status->debug >= 10)fprintf(stderr, "Error : reading complex bytes in %s, line %d\n",__FILE__,__LINE__);
                                status->status = CPHD_READ_ERROR ;
                                return (status->status) ;
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
                                        if(status->debug >= 10)fprintf(stderr, "Error : reading complex int16 in %s, line %d\n",__FILE__,__LINE__);
                                        status->status = CPHD_READ_ERROR ;
                                        return (status->status) ;
                                    }
                                    if (fread_byte_swap(&hdr->pulses[i].data.cmpl_i16[p].i, sizeof(int16_t), 1, fp, hdr->byte_swap, status) != 1) {
                                        if(status->debug >= 10)fprintf(stderr, "Error : reading complex int16 in %s, line %d\n",__FILE__,__LINE__);
                                        status->status = CPHD_READ_ERROR ;
                                        return (status->status) ;
                                    }
                                }
                            }
                            break;
                            
                        case ITYPE_CMPL_FLOAT:
                            if (fread_byte_swap(hdr->pulses[i].data.cmpl_f, sizeof(float), 2 * hdr->nsamp, fp, hdr->byte_swap, status) != 2 * hdr->nsamp) {
                                if(status->debug >= 10)fprintf(stderr, "Error : reading complex floats in %s, line %d\n",__FILE__,__LINE__);
                                status->status = CPHD_READ_ERROR ;
                                return (status->status) ;
                            }
                            break;
                            
                        default:
                            if(status->debug >= 10)printf("Unhandled type in WB sample read (%d)\n", hdr->data_type);
                            status->status = CPHD_READ_ERROR ;
                            return (status->status) ;
                            break;
                    }
                }
                    break;
                    
            }
        }
    }
    
    pthread_mutex_unlock(&hdr->file_lock);
    
    return(status->status);
}

int load_cphd3(SPImage * data, CPHDHeader * hdr, int start, int stop, SPStatus * status) {
    int i, j;
    
    CHECK_STATUS(status);
    
    if(hdr->pulses == NULL){
        if(status->debug >= 10)fprintf(stderr, "Error : NarrowBand data does not exist: %s, line %d\n",__FILE__,__LINE__);
        status->status = CPHD_READ_ERROR ;
    }
    
    if (readCPHD3Pulses(hdr, start, stop, status)) {
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

int write_cphd3_wb_vectors(CPHDHeader * hdr, SPImage * data, int pulse, FILE * out_fp, SPStatus * status) {
    
    CHECK_STATUS(status) ;

    int i;
    int chan_num = 0;
    int pulse_num;
    
    if (data->image_type != ITYPE_CMPL_FLOAT) {
        if(status->debug >= 10)printf("Data needs to be complex float write_cphd_wb_vectors (%s:%d)\n", __FILE__, __LINE__);
        status->status = INVALID_TYPE;
        return (status->status);
    }
    
    for (i = 0; i < data->ny; i++) {
        
        if (fwrite(&chan_num, sizeof(int), 1, out_fp) != 1) {
            if(status->debug >= 10)fprintf(stderr, "Write file for chan num data\n");
            status->status = CPHD_WRITE_ERROR ;
            return (status->status);
        }
        
        pulse_num = pulse + i;
        if (fwrite(&pulse_num, sizeof(int), 1, out_fp) != 1) {
            if(status->debug >= 10)fprintf(stderr, "Write file for pulse number\n");
            status->status = CPHD_WRITE_ERROR ;
            return (status->status);
        }
        
        if (fwrite(&data->data.cmpl_f[i * data->nx], sizeof(SPCmplx), data->nx, out_fp) != data->nx) {
            if(status->debug >= 10)fprintf(stderr, "Write file for IQ data\n");
            status->status = CPHD_WRITE_ERROR ;
            return (status->status);
        }
    }
    
    return (status->status);
}

int write_cphd3_nb_vectors(CPHDHeader * hdr, int channel, FILE * out_fp, SPStatus * status) {
    int idx;
    int i;
    size_t f;
    
    for(idx = 0; idx < hdr->num_azi; idx++) {
        for(i = 0; i < hdr->num_nb_items; i++) {
            switch (hdr->pos_nb[i]) {
                case   CPHD_NB_ChannelNumber:
                    f = fwrite(&channel, sizeof(int), 1, out_fp);
                    CHECK_BYTES(f, 1, status);
                    break;
                    
                case   CPHD_NB_VectorNumber:
                    f = fwrite(&hdr->pulses[idx].pulse_number, sizeof(int), 1, out_fp);
                    CHECK_BYTES(f, 1, status);
                    break;
                    
                case   CPHD_NB_SRP:
                    f = fwrite(&hdr->pulses[idx].srp.x, sizeof(double), 1, out_fp);
                    CHECK_BYTES(f, 1, status);
                    f = fwrite(&hdr->pulses[idx].srp.y, sizeof(double), 1, out_fp);
                    CHECK_BYTES(f, 1, status);
                    f = fwrite(&hdr->pulses[idx].srp.z, sizeof(double), 1, out_fp);
                    CHECK_BYTES(f, 1, status);
                    break;
                    
                case   CPHD_NB_TxPos:
                    f = fwrite(&hdr->pulses[idx].sat_ps_tx.x, sizeof(double), 1, out_fp);
                    CHECK_BYTES(f, 1, status);
                    f = fwrite(&hdr->pulses[idx].sat_ps_tx.y, sizeof(double), 1, out_fp);
                    CHECK_BYTES(f, 1, status);
                    f = fwrite(&hdr->pulses[idx].sat_ps_tx.z, sizeof(double), 1, out_fp);
                    CHECK_BYTES(f, 1, status);
                    break;
                    
                case   CPHD_NB_RcvPos:	
                    f = fwrite(&hdr->pulses[idx].sat_ps_rx.x, sizeof(double), 1, out_fp);
                    CHECK_BYTES(f, 1, status);
                    f = fwrite(&hdr->pulses[idx].sat_ps_rx.y, sizeof(double), 1, out_fp);
                    CHECK_BYTES(f, 1, status);
                    f = fwrite(&hdr->pulses[idx].sat_ps_rx.z, sizeof(double), 1, out_fp);
                    CHECK_BYTES(f, 1, status);
                    break;
                    
                case   CPHD_NB_TxTime:
                    f = fwrite(&hdr->pulses[idx].sat_tx_time, sizeof(double), 1, out_fp);
                    CHECK_BYTES(f, 1, status);
                    break;
                    
                case   CPHD_NB_RcvTime:
                    f = fwrite(&hdr->pulses[idx].sat_rx_time, sizeof(double), 1, out_fp);
                    CHECK_BYTES(f, 1, status);
                    break;
                    
                case   CPHD_NB_Fx0:
                    f = fwrite(&hdr->pulses[idx].fx0, sizeof(double), 1, out_fp);
                    CHECK_BYTES(f, 1, status);
                    break;
                    
                case   CPHD_NB_FxStepSize:
                    f = fwrite(&hdr->pulses[idx].fx_step_size, sizeof(double), 1, out_fp);
                    CHECK_BYTES(f, 1, status);
                    break;
                    
                case   CPHD_NB_AmpSF0:
                    f = fwrite(&hdr->pulses[idx].amp_sf0, sizeof(double), 1, out_fp);
                    CHECK_BYTES(f, 1, status);
                    break;
                    
                default:
                    continue;   // These items exist only in CPHD-X spec
            }
        }
    }
    
    return (status->status);
}

int write_cphd3_header(CPHDHeader * hdr,   FILE * fp){
    SPStatus status;
    im_init_status(status, 0) ;
    return (writeCPHD3Header(hdr, fp, &status) ) ;
}

int writeCPHD3Header(CPHDHeader *hdr, FILE *fp, SPStatus *status) {
    int i;
    char * data_type_str;
    unsigned int magic;
    double tmp_az, tmp_el;
    
    CHECK_STATUS(status) ;
    fprintf(fp, "BegPreamble\n");
    fprintf(fp, "Version                  3.0\n");
    fprintf(fp, "Classification           %s\n", hdr->classification);
    fprintf(fp, "DateTime                 %s\n", hdr->dateTime);
    fprintf(fp, "DataSetID                %s\n", hdr->dataSetID);
    fprintf(fp, "Mode                     %s\n", hdr->mode);
    fprintf(fp, "Geometry                 %s\n", hdr->geometry);
    
    fprintf(fp, "FixedSRP                 ");
    if (hdr->fixedSRP) {
        fprintf(fp, "Yes\n");
    } else {
        fprintf(fp, "No\n");
    }
    fprintf(fp, "Datum                    WGS84\n");
    switch(hdr->data_type){
        case ITYPE_CMPL_FLOAT:
            data_type_str = "cmplxf";
            break;
        case ITYPE_CMPL_INT16:
            data_type_str = "cmplxs";
            break;
        case ITYPE_CMPL_INT8:
            data_type_str = "cmplxb";
            break;
        default:
            data_type_str = "cmplxn";
    }
    fprintf(fp, "PHDataType               %s\n", data_type_str);
    fprintf(fp, "Interleaved              No\n");
    fprintf(fp, "Nchannels                %d\n", hdr->nchan);
    fprintf(fp, "Nvectors                 %d\n", hdr->num_azi);
    fprintf(fp, "Nsamples                 %d\n", hdr->nsamp);
    fprintf(fp, "TOASaved                 %.14g\n", hdr->TOASaved);
    fprintf(fp, "PhaseSgn                 %d\n", hdr->phaseSgn);
    fprintf(fp, "Sensor                   %s\n", hdr->sensor);
    fprintf(fp, "FreqReferenceIndex       0\n");
    
    fprintf(fp, "Grid                     Polar\n");
    
    fprintf(fp, "DeskewApplied            ");
    if (hdr->deskewed) {
        fprintf(fp, "Yes\n");
    } else {
        fprintf(fp, "No\n");
    }
    fprintf(fp, "NominalCenterFreq        %.14g\n", hdr->freq_centre);
    fprintf(fp, "NominalChirpRate         %.14g\n", hdr->chirp_gamma);
    fprintf(fp, "NominalADRate            %.14g\n", hdr->clock_speed);
    fprintf(fp, "XmitPulseDuration        %.14g\n", hdr->pulse_length);
    if (hdr->polarisation) {
        for (i = 0; i < hdr->nchan; i++) {
            fprintf(fp, "ChanDesc                 %d %s\n", i, hdr->polarisation[i]);
        }
    }
    
    for (i = 0; i < hdr->nchan; i++) {
        fprintf(fp, "AntennaPatternWidth  channel           %d\n", i);
        if (hdr->antenna_width_az <= 0.0) {
            tmp_az = 0.003660;
        } else {
            tmp_az = hdr->antenna_width_az;
        }
        
        if (hdr->antenna_width_el <= 0.0) { 
            tmp_el = 0.003660;
        } else {
            tmp_el = hdr->antenna_width_el;
        }
        fprintf(fp, "AntennaPatternWidth  HalfPowerWidthAz %.8f\n", tmp_az);
        fprintf(fp, "AntennaPatternWidth  HalfPowerWidthEl %.8f\n", tmp_el);
    }
    fprintf(fp, "BegNBVector\n");
    for(i = 0; i < hdr->num_nb_items; i++) {
        switch (hdr->pos_nb[i]) {
            case CPHD_NB_ChannelNumber:
            case CPHD_NB_VectorNumber:
            case CPHD_NB_SRP:
            case CPHD_NB_TxPos:
            case CPHD_NB_RcvPos:
            case CPHD_NB_TxTime:
            case CPHD_NB_RcvTime:
            case CPHD_NB_Fx0:
            case CPHD_NB_FxStepSize:
            case CPHD_NB_AmpSF0:
                
                fprintf(fp, "           %s\n", nb_names[hdr->pos_nb[i]]);
                break;
                
            default:
                break;
        }
    }
    fprintf(fp, "EndNBVector\n");
    
    fprintf(fp, "BegWBVector\n");
    for(i = 0; i < hdr->num_wb_items; i++) {
        fprintf(fp, "           %s\n", wb_names[i]);
    }
    fprintf(fp, "EndWBVector\n");
    
    fprintf(fp, "EndPreamble\n");
    
    magic = CPHD_MAGIC_NUMBER;
    
    fwrite(&magic, sizeof(unsigned int), 1, fp);
    
    return status->status;
}

int writeCPHD3File(CPHDHeader *hdr, SPStatus * status){
    
    int chan_num = 0;
    int pulse_num ;
    
    CHECK_STATUS(status) ;
    
    if (hdr->fp == NULL) {
        printf("Error: File pointer not set in cphd data structure\n");
        status->status = CPHD_WRITE_ERROR ;
        return status->status ;
    }
    
    writeCPHD3Header(hdr, hdr->fp, status);
    CHECK_STATUS(status) ;
    write_cphd3_nb_vectors(hdr, 1, hdr->fp, status) ;
    CHECK_STATUS(status) ;
    
    for (int i = 0; i < hdr->num_azi; i++) {
        
        if (fwrite(&chan_num, sizeof(int), 1, hdr->fp) != 1) {
            if(status->debug >= 10)fprintf(stderr, "Write file for chan num data\n");
            status->status = CPHD_WRITE_ERROR ;
            return (status->status);
        }
        
        pulse_num = i;
        if (fwrite(&pulse_num, sizeof(int), 1, hdr->fp) != 1) {
            if(status->debug >= 10)fprintf(stderr, "Write file for pulse number\n");
            status->status = CPHD_WRITE_ERROR ;
            return (status->status);
        }
        
        if (fwrite(&(hdr->pulses[i].data.cmpl_f[0]), sizeof(SPCmplx), hdr->nsamp, hdr->fp) != hdr->nsamp) {
            if(status->debug >= 10)fprintf(stderr, "Write file for IQ data\n");
            status->status = CPHD_WRITE_ERROR ;
            return (status->status);
        }
    }
    
    return status->status ;
}

