/***************************************************************************
 *
 *       Module:    colourCodes.h
 *      Program:    tdpocl
 *   Created by:    Darren Muff on 12/03/2013.
 *                  Copyright (c) 2013 [Dstl]. All rights reserved.
 *
 *   Description:
 *   Header file to define terminal colour codes that can be used in fancy
 *   printf output
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

#ifndef GOSS_colourCodes_h
#define GOSS_colourCodes_h

#define RESETCOLOR  "\033[0m"
#define BLACK       "\033[30m"
#define RED         "\033[31m"
#define GREEN       "\033[32m"
#define YELLOW      "\033[33m"
#define BLUE        "\033[34m"
#define MAGENTA     "\033[35m"
#define CYAN        "\033[36m"
#define WHITE       "\033[37m"
#define BLINK       "\033[5m"
#define BLINKOFF    "\033[25m"
#define NORMAL      "\033[22m"
#define BOLD        "\033[1m"
#define DARK        "\033[1m"
#define FEINT       "\033[2m"
#define LIGHT       "\033[2m"
#define ITALIC      "\033[3m"
#define UNDERLINE   "\033[4m"
#define CROSSEDOUT  "\033[9m"


#endif
