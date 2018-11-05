/***************************************************************************
 * 
 *           Module :  Timer.h
 *          Program :  sarclib
 *       Created by :  Darren Muff on 22/07/2013
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      Simple routines for calculating the time of a task
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

#ifndef sarclib_Timer_h
#define sarclib_Timer_h

#include "sarclib.h"

typedef struct Timer {
    double start ;
    double end ;
} Timer ;

/// Start a timer
///
void startTimer(Timer *t, SPStatus *status);
/// End the timer
///
void endTimer(Timer *t, SPStatus * status);
/// return the time is seconds since the timer was started
///
double timeElapsedInSeconds(Timer *t, SPStatus *status);
/// return the time in milliseconds since the timer was started
///
double timeElapsedInMilliseconds(Timer *t, SPStatus *status);
/// return the time in minutes since the timer was started
///
double timeElapsedInMinutes(Timer *t, SPStatus *status);
/// returns the time remaining until the task is completed
/// based upon the percentage of the task that has been completed
///
double estimatedTimeRemaining(Timer *t, double percentDone, SPStatus *status) ;


#endif
