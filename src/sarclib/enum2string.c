/***************************************************************************
 *
 *       Module:    enum2string.h
 *      Program:    sarclib
 *   Created by:    Emma Griffiths on 20/07/2006.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file is here to simplify those times when you need to be able to
 *      access an enum type as a string.
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


#include "enum2string.h"

// converts an enum representing a direction into a string
//
char *dir2string(SPDirection dir)
{
    switch (dir)
    {
        case X_DIR:
            return("x direction");
            break;
        case Y_DIR:
            return("y direction");
            break;
        case Z_DIR:
            return("z direction");
            break;
    }
    return("broken direction in dir2string");
}

// converts a cast mode to a string
//
char *castmode2string(SPCastMode mode)
{
    switch (mode)
    {
        case CAST_MODULUS:
            return("cast modulus");
            break;
        case CAST_PHASE:
            return("cast phase");
            break;
        case CAST_UNIT:
            return("cast unit");
            break;
        case CAST_REAL:
            return("cast real");
            break;
        case CAST_IMAG:
            return("cast imag");
            break;
            
        case SET_REAL:
            return("set real");
            break;
        case SET_IMAG:
            return("set imag");
            break;
            
        case SET_MODULUS:
            return("set modulus");
            break;
        case SET_PHASE:
            return("set phase");
            break;
    }
    return("broken cast mode in castmode2string");
}

// converts an enum representing an image type into a string
//
char *itype2string(SPImageType image_type) {
    switch (image_type)
    {
        case ITYPE_UNKNOWN:
            return("unknown");
            break;
        case ITYPE_POLAR:
            return("polar");
            break;
        case ITYPE_CMPL_FLOAT:
            return("complex float");
            break;
        case ITYPE_FLOAT:
            return("float");
            break;
        case ITYPE_DOUBLE:
            return("double");
            break;
        case ITYPE_INT64:
            return("int64");
            break;
        case ITYPE_INT32:
            return("int32");
            break;
        case ITYPE_INT16:
            return("int16");
            break;
        case ITYPE_INT8:
            return("int8");
            break;
        case ITYPE_UINT64:
            return("uint64");
            break;
        case ITYPE_UINT32:
            return("uint32");
            break;
        case ITYPE_UINT16:
            return("uint16");
            break;
        case ITYPE_UINT8:
            return("uint8");
            break;
        case ITYPE_VECTOR:
            return("vector");
            break;
        case ITYPE_CMPL_INT64:
            return("complex int64");
            break;
        case ITYPE_CMPL_INT32:
            return("complex int32");
            break;
        case ITYPE_CMPL_INT16:
            return("complex int16");
            break;
        case ITYPE_CMPL_INT8:
            return("complex int8");
            break;
        case ITYPE_CMPL_NIBBLE:
            return("complex nibble");
            break;
    }
    return("broken image type in itype2string");
}
