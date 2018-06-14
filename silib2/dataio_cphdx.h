/***************************************************************************
 *
 *       Module:    dataio_cphdx.h
 *      Program:    SIlib2
 *   Created by:    Matt Nottingham on 08/09/2009.
 *                  and Darren Muff 24/07/2013
 *                  Copyright (c) 2017 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains functions for reading Compensated Phase History
 *      DATA (CPHD) formatted data in cphdx format
 *
 *   CLASSIFICATION        :  OFFICIAL
 *   Date of CLASSN        :  24/07/2017
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

#ifndef dataio_cphdx_h
#define dataio_cphdx_h


// Enumeration of Required types for CPHDv0.3
// The specification allows conditionals that are not yet implemented
// This implementation supports the following conditionals:
//  DomainType : FX
//  SRPType : STEPPED
//  NUMSRPS : 1 (Prolly a bug in the spec as 'STEPPED' allows SRPs to change per pulse/vector)
//
#define NUMCPHDXRequiredItems 32
typedef enum {
    CollectorName,      // Radar platform identifier.
    CoreName,           // Collection & imaging data set identifier
    ModeType,           // Radar imaging mode.
    Classification,     // Text based field containing human- readable banner
    SampleType,         // Indicates the PHD sample format of the PHD array(s).
                        // All arrays have the sample type. Real and imaginary
                        // components stored in adjacent bytes, real component stored first.
                        // All binary data parameters are stored with Most Significant Byte
                        // first (i.e. "Big Endian").
                        // Allowed Values: "RE08I_IM08I", "RE16I_IM16I", "RE32F_IM32F"
    NumCPHDChannels,    // Number of CPHD channels stored in the product file.
    NumBytesVBP,        // Number of bytes per set of Vector Based Parameters. One set of
                        // VBPs for each CPHD vector.
    NumVectors,         // Number of slow time vectors in the PHD array for channel n
    NumSamples,         // Number of samples per vector in the PHD array for channel n
    DomainType,         // Indicates the domain represented by the sample dimension of the
                        // PHD arrays. All CPHD channels are the same domain type.
    PhaseSGN,           // Phase SGN associated target phase and positive deltaTOA
    CollectStart,       // Collection start date and time (UTC) on UTC the Collector
                        // platform. Format is YYYY-MM-DDThh:mm:ss.s+
    CollectDuration,    // Collection duration (sec) for the Collector platform The CPHD array(s)
                        // in the file may or may not span the whole collection.
    TxTime1,            // Earliest TxTime for any CPHD vector in the file. Time relative to collection start.
    TxTime2,            // Latest TxTime for any CPHD vector in the file. Time relative to collection start.
    ACP,                // Corner point parameters. Corners indexed x = 1, 2, 3, 4 clockwise.
                        // -90.0 < Lat < 90.0, -180.0 < Lon < 180.0
    SRP_Index,          // Index to identify the SRP position - function used for the channel.
                        // 'Required' but not used as SRPType is STEPPED
    NomTOARateSF,       // Scale factor to indicate the fraction of
                        // the Doppler spectrum that is clear.
    FxCtrNom,           // Nominal center transmit frequency associated with the channel (Hz).
    BWSavedNom,         // Nominal transmit bandwidth associated with the channel (Hz).
    TOASavedNom,        // Nominal span in TOA saved for the channel. TOASavedNom is the bandwidth
                        // saved for all vectors.
    SRPType,            // STEPPED: The SRP positions are not constant and are specified in each vector
    NumSRPs,            // Number of unique SRP position functions used. Required but must be '1'
    TxTime,             // Number of bytes for TxTime in per-vector Narrowband data (8)
    TxPos,              // Number of bytes for TxPos in per-vector Narrowband data (3x8 = 24)
    RcvTime,            // Number of bytes for RxTime in per-vector Narrowband data (8)
    RcvPos,             // Number of bytes for RxPos in per-vector Narrowband data (3x8 = 24)
    SRPPos,             // Number of bytes for SRPPos in per-vector Narrowband data (3x8 = 24)
    Fx0,                // Number of bytes for Fx0 in per-vector Narrowband data (8)
    Fx_SS,              // Number of bytes for Fx_ss in per-vector Narrowband data (8)
    Fx1,                // Number of bytes for Fx1 in per-vector Narrowband data (8)
    Fx2                 // Number of bytes for Fx2 in per-vector Narrowband data (8)
} CPHDXRequired;


static const char *CPHDXReqdElements[] = {
    "CollectorName","CoreName","ModeType","Classification","SampleType","NumCPHDChannels","NumBytesVBP","NumVectors","NumSamples",
    "DomainType","PhaseSGN","CollectStart","CollectDuration","TxTime1","TxTime2","ACP","SRP_Index","NomTOARateSF","FxCtrNom",
    "BWSavedNom","TOASavedNom","SRPType","NumSRPs","TxTime","TxPos","RcvTime","RcvPos","SRPPos","Fx0","Fx_SS","Fx1","Fx2"
};


#define BUFFSIZE        8192

char Buff[BUFFSIZE];
char Names[32][256];

int Depth;

struct SPXMLListEle * cphdXMLlist = NULL;

#define  writeToBufParm(ptr, text, param){                  \
    char *__tmp = (char *)malloc(sizeof(char) * 1024);      \
    size_t __len ;                                          \
    sprintf(__tmp,text,param);                              \
    memcpy(ptr, __tmp, __len = strlen(__tmp));              \
    ptr += __len;                                           \
    free(__tmp);                                            \
}

#define  writeToBufText(ptr, text){                         \
    char *__tmp = (char *)malloc(sizeof(char) * 1024);      \
    size_t __len ;                                          \
    sprintf(__tmp,text);                                    \
    memcpy(ptr, __tmp, __len = strlen(__tmp) ) ;            \
    ptr += __len ;                                          \
    free(__tmp) ;                                           \
}


#endif /* dataio_cphdx_h */
