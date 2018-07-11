/***************************************************************************
 *
 *       Module:    banner.c
 *      Program:    tdpocl
 *   Created by:    Darren Muff on 12/03/2013.
 *                  Copyright (c) 2013 [Dstl]. All rights reserved.
 *
 *   Description:
 *   Simple banner for the programme that uses some useful git version output
 *
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  25/04/2014
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

#include <stdio.h>
#include "tdpoclVersion.h"
#include "colourCodes.h"

void banner (){
    printf(" \n");
    printf(DARK GREEN "              Time Domain Processor OpenCL Variant (TDPOCL)\n" NORMAL);
    printf(BLUE "                            Version :"RED" %s "BLUE"(Rev "RED"%s"BLUE")\n", SHORT_VERSION,REVISION);
    printf(BLUE "           Commit: "RED"%s, %s\n",FULL_VERSION, VERSION_DATE);
    printf(BLUE "               Copyright (c) 2014 "WHITE"["BLUE"Dstl"WHITE"]"BLUE". All rights reserved.\n" RESETCOLOR);
    printf(" \n");

    return ;
}