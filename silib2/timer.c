/***************************************************************************
 *
 *       Module:    timer.c
 *      Program:    SILib2
 *   Created by:    Darren Muff on 21/07/2013.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      Simple routines for calculating the time of a task
 *
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  14/03/2013
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

#include "Timer.h"
#include <sys/time.h>

void startTimer(Timer *t, SPStatus *status){
    CHECK_STATUS(status) ;

    struct timeval tv;
    gettimeofday(&tv, NULL);
    t->start = (double)tv.tv_sec+((double)tv.tv_usec / 1.0e6);
}

void endTimer(Timer *t, SPStatus * status){
    CHECK_STATUS(status) ;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    t->end = (double)tv.tv_sec+((double)tv.tv_usec / 1.0e6);
}

double timeElapsedInSeconds(Timer *t, SPStatus *status){
    CHECK_STATUS(status) ;
    return (t->end - t->start) ;
}

double timeElapsedInMilliseconds(Timer *t, SPStatus *status){
    CHECK_STATUS(status) ;
    return ((t->end - t->start)*1000.0) ;
}

double timeElapsedInMinutes(Timer *t, SPStatus *status){
    CHECK_STATUS(status) ;
    return ((t->end - t->start)/60) ;
}

double estimatedTimeRemaining(Timer *t, double percentDone, SPStatus *status) {
    CHECK_STATUS(status) ;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double now = (double)tv.tv_sec+((double)tv.tv_usec / 1.0e6) ;
    
    double durationToNow = now - t->start ;
    double totalDuration = durationToNow * 100 / percentDone ;
    double endTime = t->start + totalDuration ;
    double timeRemaining = endTime - now ;
    return timeRemaining ;
}
