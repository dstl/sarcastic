/***************************************************************************
 * 
 *           Module :  alloc.c
 *          Program :  sarclib
 *       Created by :  Darren Muff on Sat Jun 30 10:34:25 2018 +0100
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *     Simple memory allocation functions
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

#include "alloc.h"

void *
sp_calloc_fn(size_t count, size_t size, const char * fname, int line)
{
    void * ret;
    
    ret = calloc(count,size);
    
    if (! ret) {
      fprintf(stderr, "Failed to alloc %ld bytes of memory at %s:%d\n", count * size, fname, line);
        exit(67);
    }
    
    return ret;
}

void *
sp_malloc_fn(size_t size, const char * fname, int line)
{
    void * ret;
    
    ret = malloc(size);
    
    if (! ret) {
      fprintf(stderr, "Failed to alloc %ld bytes of memory at %s:%d\n", size, fname, line);
        exit(67);
    }
    
    return ret;
}

