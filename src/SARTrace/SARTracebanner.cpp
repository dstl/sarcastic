/***************************************************************************
 * 
 *           Module :  SARTracebanner.cpp
 *          Program :  SARTrace
 *       Created by :  Darren Muff on 12/03/2013
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *         Prints out useful info about the programme and version
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

#include <stdio.h>
#include "SartraceVersion.h"
#include "colourCodes.h"

void SARTracebanner (){
    printf(" \n");
    printf(DARK GREEN "                   SARTrace - SAR Ray Tracer\n" NORMAL);
    printf(BLUE       "                Version :" RED" %s \n", FULL_VERSION);
    printf(BLUE       "            Revision: " RED"%s, %s \n",REVISION, VERSION_DATE);
    printf(BLUE       "         Copyright (c) 2015 " WHITE"[" BLUE"Dstl" WHITE"]" BLUE". All rights reserved.\n" RESETCOLOR);
    printf(" \n");

    return ;
}
