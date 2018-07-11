/** @file********************************************************************
 *
 *       Module:    dataio.h
 *      Program:    sarclib
 *
 *   Created by:    Emma Griffiths on 01/02/2005.
 *                  and Darren Muff 7th Jun 2017
 *                  Copyright (c) 2017 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains functions for reading and writing to and from files to
 *      the variable types used in the time-domain SAR processor
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  7/06/2017
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

#ifndef sarclib_DATAIO_H__
#define sarclib_DATAIO_H__

#include "dataio_au4.h"
#include "dataio_c8.h"
#include "dataio_cphd.h"
#include "dataio_envi.h"
#include "dataio_gdal.h"
#include "dataio_hdf4.h"
#include "dataio_im.h"
#include "dataio_sicd.h"
#include "dataio_sio.h"
#include "dataio_srf.h"

/// Reads nitems of size 'size' from file stream 'stream' and writes them to 'ptr'.
/// If 'do_swap' is true then perform a byte swap on read. Otherwise just read in
///
size_t fread_byte_swap(void * ptr, size_t size, size_t nitem, FILE * stream, int do_swap, SPStatus *status);

/// Writes nitems of size 'size' from 'ptr' to file stream 'stream'.
/// If 'do_swap' is true then perform a byte swap on read. Otherwise just write the data
///
size_t fwrite_byte_swap(void * ptr, size_t size, size_t nitem, FILE * stream, int do_swap, SPStatus *status) ;

/// byte swaps 'nitem' of scalers of size 'size'
///
void im_swap_bytes(void * ptr, size_t size, size_t nitem);

/// allows the user to input an int at a prompt, or select a default
///  \param prompt  : the text to display when asking the user for input (the default value will also be displayed)
///  \param key     : a unique bit of text (such az "AzPixelSize") which the input system will then use to store parameter values
///  \param help    : the text to display to the user when they just enter '?'
///  \param def     : the default, ie the value to take if the user just presses return
///
int input_int (const char *prompt, const char * key, const char * help, int def);

/// allows the user to input an int64 at a prompt, or select a default *
///  \param prompt  : the text to display when asking the user for input (the default value will also be displayed)
///  \param key     : a unique bit of text (such az "AzPixelSize") which the input system will then use to store parameter values
///  \param help    : the text to display to the user when they just enter '?'
///  \param def     : the default, ie the value to take if the user just presses return
///
int64_t input_int64 (const char *prompt, const char * key, const char * help, int64_t def);

/// allows the user to input a yes (=1)/no(=0) at a prompt, or select a default
///  \param prompt  : the text to display when asking the user for input (the default value will also be displayed)
///  \param key     : a unique bit of text (such az "AzPixelSize") which the input system will then use to store parameter values
///  \param help    : the text to display to the user when they just enter '?'
///  \param def     : the default, ie the value to take if the user just presses return
///
int input_yesno (const char *prompt, const char * key, const char * help, int def);

/// allows the user to input an double at a prompt, or select a default
///  \param prompt  : the text to display when asking the user for input (the default value will also be displayed)
///  \param key     : a unique bit of text (such az "AzPixelSize") which the input system will then use to store parameter values
///  \param help    : the text to display to the user when they just enter '?'
///  \param def     : the default, ie the value to take if the user just presses return
///
double input_dbl (const char *prompt, const char * key, const char * help, double def);

/// allows the user to input an SPVector at a prompt, or select a default
///  \param prompt  : the text to display when asking the user for input (the default value will also be displayed)
///  \param key     : a unique bit of text (such az "AzPixelSize") which the input system will then use to store parameter values
///  \param help    : the text to display to the user when they just enter '?'
///  \param def     : the default, ie the value to take if the user just presses return
///
SPVector input_vect (const char *prompt, const char * key, const char * help, SPVector def);

/// allows the user to input a string at a prompt, or select a default
///  \param prompt  : the text to display when asking the user for input (the default value will also be displayed)
///  \param key     : a unique bit of text (such az "AzPixelSize") which the input system will then use to store parameter values
///  \param help    : the text to display to the user when they just enter '?'
///  \param def     : the default, ie the value to take if the user just presses return
///
char * input_string (const char *prompt, const char * key, const char * help, const char *def);      

#endif

