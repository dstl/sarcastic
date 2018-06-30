/** @file********************************************************************
 *
 *       Module:    alloc.h
 *      Program:    sarclib
 *   Created by:    Matt Nottingham on 16/06/2013.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *     Simple memory allocation functions
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

#ifndef sarclib_ALLOC_H__
#define sarclib_ALLOC_H__

#include <stdio.h>
#include <stdlib.h>

void * sp_calloc_fn(size_t count, size_t size,  const char * fname, int line);
void * sp_malloc_fn(size_t size,  const char * fname, int line);

#define sp_calloc(n, s) sp_calloc_fn(n, s, __FILE__, __LINE__)
#define sp_malloc(n)    sp_malloc_fn(n, __FILE__, __LINE__)


#endif
