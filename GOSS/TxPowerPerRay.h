/***************************************************************************
 *
 *       Module:    TxPowerPerRay.h
 *      Program:    Sadilac
 *   Created by:    Darren on 16/06/2013.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *     Calculate the power for a given ray when a radar beam is split  into
 *     xRays x yRays number of rays over a solid angle of xBeamUsed x yBeamUsed
 *     specified in radians.
 *
 *     Only part of the antenna beam is used to render the scene. The directivity
 *     and therefore gain however is always defined in terms of a ratio per unit
 *     solid angle.  This defines the power gain in any given direction. The beam
 *     is quantised however into discrete rays. The gain therefore needs to be
 *     quantised into rays. Failure to do this results in higher power values for
 *     more rays which skews RCS values.
 *
 *   CLASSIFICATION        :  <PENDING>
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

#ifndef _GOSS_TxPowerPerRay_h_
#define _GOSS_TxPowerPerRay_h_

const double AntLen             = 11.0;         // Antenna Length in metres
const double AntHei             = 11.0;         // Antenna height in metres
const double Pt                 = 150;          // Peak power in Watts
const double ApEff              = 70.0;         // Aperture efficiency in %
const double lambda             = 0.04;         // wavelength for gain purposes

double TxPowerPerRay(int xRays, int yRays, double xBeamUsed, double yBeamUsed, double * raySolidAngle, double *effectiveArea);

#endif
