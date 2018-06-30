/** @file********************************************************************
 *
 *       Module:    dataio_sio.h
 *      Program:    sarclib
 *   Created by:    Matt Nottingham on 17/07/2006.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      Description: Files for reading/writing SIO format files.
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

#ifndef sarclib_DATAIO_SIO_H__
#define sarclib_DATAIO_SIO_H__

#define SIO_MAGIC_BASIC     (0xff017ffe)
#define SIO_MAGIC_EXTENDED  (0xff027ffd)

typedef struct
{
  FILE *  fp;              ///< The file pointer  
  int     magic;           ///< The SIO magic number  
  int     nl;              ///< The number of lines (ny)  
  int     ne;              ///< The number of elements per line (nx)  
  int     d_type;          ///< The data type ID  
  int     d_size;          ///< sizeof(data type)  
  int     do_byte_swap;    ///< Whether we need to swap bytes due to little/big endian issues  
  char ** field_labels;    ///< If an extended SIO file, then an array of strings with the field labels in  
  off_t  data_start;       ///< Position in the file where the actual data starts  
} SIOheader;

SPStatus * im_open_sio(SIOheader * hdr, const char * filename, SPStatus * status);
SPStatus * im_read_header_sio(SIOheader * hdr, SPStatus * status);
SPStatus * im_load_sio_subset(SPImage * data, SIOheader * hdr, int startx, int starty, SPStatus * status);
SPStatus * im_load_sio(SPImage * data, SIOheader * hdr, SPStatus * status);
SPStatus * im_close_sio(SIOheader * hdr, SPStatus * status);
SPStatus * im_destroy_sio(SIOheader * hdr, SPStatus * status);
SPStatus * im_sio_dtype_to_itype(SIOheader * hdr, SPImageType * itype, SPStatus * status);


#endif

