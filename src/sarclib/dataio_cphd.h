/***************************************************************************
 * 
 *           Module :  dataio_cphd.h
 *          Program :  sarclib
 *       Created by :  Matt Nottingham on 25/08/2009
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      This file contains functions for reading Compensated Phase History
 *      DATA (CPHD) formatted data.
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

#ifndef sarclib_DATAIO_CPHD_H__
#define sarclib_DATAIO_CPHD_H__

#include "sarclib.h"

#define CPHD_MAX_NB_ITEMS  (64)
#define CPHD_MAX_WB_ITEMS  (16)
#define CPHD_MAGIC_NUMBER (0x4d6f6a6f)
#define CPHD_MAGIC_NUMBER_BYTE_SWAPPED (0x6f6a6f4d)

typedef enum {
    CPHD_NB_ChannelNumber,
    CPHD_NB_VectorNumber,
    CPHD_NB_SRP,
    CPHD_NB_SRPTime,
    CPHD_NB_TxPos,
    CPHD_NB_RcvPos,
    CPHD_NB_TxTime,
    CPHD_NB_RcvTime,
    CPHD_NB_Fx0,
    CPHD_NB_FxStepSize,
    CPHD_NB_Fx1,
    CPHD_NB_Fx2,
    CPHD_NB_AmpSF0
} CPHD_NB_items;

typedef enum {
    CPHD_WB_ChannelNumber,
    CPHD_WB_VectorNumber,
    CPHD_WB_Samples
} CPHD_WB_items;

typedef struct
{
    int pulse_number;               ///< Index of this pulse
    double sat_tx_time;             ///< Time of transmit relative to collection time in DateTime keyword (sec)
    double sat_rx_time;             ///< Time of receive relative to DateTime keyword (sec)
    SPVector sat_ps_tx;             ///< Transmit aperture position for this vector.
    SPVector sat_ps_rx;             ///< Receive aperture position for this vector.
    SPVector srp;                   ///< Position of Stabilization Reference Point for this vector.
    double srp_time;
    double fx0;                     ///< frequency associated with sample m=0 of this vector (Hz)
    double fx1;
    double fx2;
    double fx_step_size;            ///< Frequency sample spacing of this vector (Hz)
    double amp_sf0;                 ///< Optional amplitude scale factor for per-vector scaling of wideband data.
    
    union {
        SPCmplxInt8  * cmpl_i8;
        SPCmplxInt16 * cmpl_i16;
        SPCmplx      * cmpl_f;
        void         * v;
    } data;
    
} CPHDPulse;

/// Structure define the CPHDHeader
///
///
typedef struct
{
    struct SPXMLListEle * XMLlist;
    
    FILE * fp;
    pthread_mutex_t file_lock;
    char * sensor;                  ///< Radar sensor or system identifier. A name chosen to uniquely ID the data source
    char * dateTime;                ///< Precise only to seconds in UTC - more precise times in NBVector block
                                    ///< Date and time of start of collection (UTC). Format is YYYYMMDDHHMMSS.
    char * dataSetID;               ///< Data set identifier. The format varies and is generally unique to a platform
    char * classification;          ///< classification of the CPHDFile or system.
    char * mode;
    char * geometry;
    char ** polarisation;
    int fix_up_sensor;
    SPCmplx dc_bias;
    CPHDPulse *pulses ;
    int pos_nb[CPHD_MAX_NB_ITEMS];
    int pos_wb[CPHD_MAX_WB_ITEMS];
    SPVector grp;
    size_t pulse_size;
    double chirp_gamma;             ///< FM rate of the transmit waveform (Hz/Sec) for Linear FM waveforms.
    double clock_speed;             ///< Nominal AtoD Rate in Samples per Second
    double freq_centre;             ///< Centre RF frequency of the WB data in the file (Hz).
    ///< This may be exact, or relative (see FreqReferenceIndex keyword).
    double pulse_length;            ///< Transmit pulse length (sec).
    double TOASaved;                ///< Span in Time of Arrival (TOA) saved (seconds) for WB data in file
    double antenna_width_az;        ///< Antenna azimuth half power width in radians
    double antenna_width_el;        ///< Antenna elevation half power width in radians
    SPImageType data_type;          ///< Data type for dataset
    int deskewed;                   ///< Has deskew been applied to the data to correct for residual video phase (RVP)
    int nchan;                      ///< Number of data channels contained in file.
    ///< Channels are indexed by c=0, 1, 2, ...., Nchannels-1.
    int nsamp;                      ///< Number of complex data samples per vector. Samples are
    ///< indexed by m = 0, 1, 2, ...., Nsamples-1.
    int isamp;
    int num_azi;                    ///< Number of NB and WB vectors per channel. Vectors are
    ///< indexed by n=0, 1, 2, ......, Nvectors-1.
    int byte_swap;
    int num_nb_items;
    int num_wb_items;
    int nb_vector_size;
    int fixedSRP;
    int phaseSgn;                   ///< Defines the sign of the phase vs frequency for targets
    ///< off of the SRP. Allowed values are -1 or +1.
    int interleaved;                ///< ONLY in Version 3.0, Defines format for Once per Vector parameters.
    ///< Must precede once per Vector marker tags. Allowed values - Yes, No
    ///< Yes - Once per vector NB data is interleaved with pulses.
    ///< No - NB data is in a single NBVector block.
    off_t start_nb_data;
    off_t start_wb_data;
    
    char version;                   ///< currently either 'x' or '3'
    off_t xml_data_size;
    off_t xml_byte_offset;
    off_t vb_data_size;
    off_t vb_byte_offset;
    off_t cphd_data_size;
    off_t cphd_byte_offset;
    
} CPHDHeader;

#define CPHD_MAX_LINE_LENGTH (256)

// Generic cphd data handlers
//

/// Read in just the header information from the cphd file specified in filename
/// This includes reading all the narrowband information for the CPHD pulses
///
int readCPHDHeader(const char * filename, CPHDHeader * hdr, SPStatus *status);

/// Read wideband (pulse) data from a file associated with a CPHD Header file
/// and place it in hdr pulse data structure
///
int readCPHDPulses(CPHDHeader * hdr, int start_pulse, int stop_pulse, SPStatus *status);

/// Read wideband (pulse) data from a file associated with a CPHD header file
/// (using read_cphd_pulses) and then use the data to create an SPImage dataset
///
int load_cphd(SPImage * data, CPHDHeader * hdr, int start_pulse, int stop_pulse, SPStatus * status);

/// Destroy a cphd header structure and reclaim allocated memory
///
int destroy_cphd_header(CPHDHeader * hdr, SPStatus * status);

/// Destroy pulse data associated with a cphd header and reclaim memory
///
int destroy_cphd_pulses(CPHDHeader * hdr, int start_pulse, int last_pulse, SPStatus * status);

int writeCPHDFile(CPHDHeader *hdr, SPStatus * status) ;


// CPHD3 version data handling functions
//

/// Write the narrowband information to a cphd3 file
///
int write_cphd3_nb_vectors(CPHDHeader * hdr, int channel, FILE * out_fp, SPStatus * status);

/// write the cphd3 header information to a file
///
int writeCPHD3Header(CPHDHeader * hdr, FILE * fp, SPStatus *status);

/// Write the cphd3 wideband data to a file
///
int write_cphd3_wb_vectors(CPHDHeader * hdr, SPImage *data, int pulse, FILE * out_fp, SPStatus * status);

/// Read the header and narrowband data of a CPHD3 file. Used of the file is known to be CPHD3.
/// If the file type is unknown then tuse the function readCPHDHeader()
///
int readCPHD3Header(CPHDHeader * hdr, SPStatus *status);

/// Read in the wideband pulse data from a CPHD3 file. If the file type is unknwown then use
/// the function readCPHDPulses()
///
int readCPHD3Pulses(CPHDHeader * hdr, int start_pulse, int stop_pulse, SPStatus *status);

/// Load in the wideband data form a CPHD3 file and store it as an SPImage. If teh data type is unknown
/// then use the function load_cphd()
///
int load_cphd3(SPImage * data, CPHDHeader * hdr, int start_pulse, int stop_pulse, SPStatus * status);
int readCPHDXHeaderV0p3(CPHDHeader *  hdr, SPStatus *status);
int readCPHDXPulsesV0p3(CPHDHeader * hdr, int start_pulse, int stop_pulse, SPStatus *status);
int load_cphdx(SPImage * data, CPHDHeader * hdr, int start_pulse, int stop_pulse, SPStatus * status);

/// write the cphdX header information to a file
///
int writeCPHDXHeader(CPHDHeader * hdr, FILE * fp, SPStatus * status) ;

/// Write the narrowband information to a cphdX file
///
int writeCPHDXNarrowband(CPHDHeader * hdr, int channel, FILE * out_fp, SPStatus * status) ;

/// Write wideband data stored in the SPImage 'data' to a cphdX file given by 'out_fp'
///
int writeCPHDXWideband(CPHDHeader * hdr, SPImage * data, FILE * out_fp, SPStatus * status) ;

/// Version specific function to write a cphd data structure as a cphdx (v0.3) formatted file
/// The header must be correctly populated with narrowband and wideband data and
/// the file pointer correctly set.
///
int writeCPHDXFile(CPHDHeader *hdr, SPStatus * status) ;

/// Version specific function to write a cphd data structure as a cphd3 formatted file
/// The header must be correctly populated with narrowband and wideband data and
/// the file pointer correctly set.
///
int writeCPHD3File(CPHDHeader *hdr, SPStatus * status) ;

    
// Deprecated CPHD Function calls
//
/// Read in just the header information from the cphd file specified in filename
/// This includes reading all the narrowband information for the CPHD pulses
/// *** This function is deprecated - please use readCPHDHeader ***
///
int read_cphd_header(const char * filename, CPHDHeader * hdr) __attribute__((deprecated("Function deprecated use 'readCPHDHeader' instead !!!"))) ;
/// Read wideband (pulse) data from a file associated with a CPHD Header file
/// and place it in hdr pulse data structure
/// *** This function is deprecated - please use readCPHDPulses ***
///
int read_cphd_pulses(CPHDHeader * hdr, int start_pulse, int stop_pulse) __attribute__((deprecated("Function deprecated use 'readCPHDPulses' instead !!!"))) ;
/// write the cphd3 header information to a file
/// *** This function is deprecated - please use writeCPHD3Header ***
///
int write_cphd3_header(CPHDHeader * hdr, FILE * fp) __attribute__((deprecated("Function deprecated use 'writeCPHD3Header' instead !!!"))) ;
int read_cphd3_header(CPHDHeader * hdr) __attribute__((deprecated("Function deprecated use 'readCPHD3Header' instead !!!"))) ;
int read_cphd3_pulses(CPHDHeader * hdr, int start_pulse, int stop_pulse) __attribute__((deprecated("Function deprecated use 'readCPHD3Pulses' instead !!!"))) ;
int read_cphdx_header(CPHDHeader *  hdr) __attribute__((deprecated("Function deprecated use 'readCPHDXHeader' instead !!!"))) ;
int read_cphdx_pulses(CPHDHeader * hdr, int start_pulse, int stop_pulse) __attribute__((deprecated("Function deprecated use 'readCPHDXPulses' instead !!!"))) ;

#endif

