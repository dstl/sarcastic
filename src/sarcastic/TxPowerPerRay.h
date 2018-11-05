/***************************************************************************
 * 
 *           Module :  TxPowerPerRay.h
 *          Program :  sarcastic
 *       Created by :  Darren Muff on 16/06/2013
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
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

#ifndef _SARC_TxPowerPerRay_h_
#define _SARC_TxPowerPerRay_h_

// const int monoStatic            = 1;            // 1 if monostatic, 0 if bistatic
// const int TxDishAntenna         = 1;            // Is Tx a dish ?
// const int RxDishAntenna         = 1;            // Is Rx a dish ?
// const double TxAntLen           = 11.0;         // Tx Antenna Length in metres. If Tx is dish just uses length
// const double TxAntHei           = 11.0;         // Tx Antenna height in metres. Ignored if Tx is dish
// const double RxAntLen           = 11.0;         // Antenna Length in metres. If Rx is dish just uses length
// const double RxAntHei           = 11.0;         // Antenna height in metres. Ignored if Rx is dish
// const double Pt                 = 300;          // Peak power in Watts
// const double ApEffRx            = 70.0;         // Aperture efficiency in %
// const double ApEffTx            = 70.0;         // Aperture efficiency in %
// const double lambda             = 0.04;         // wavelength for gain purposes

// double TxPowerPerRay(double rayWidthRadians, double rayHeightRadians, double *receiverGain);
double TxPowerPerRay(CPHDHeader *hdr, double *receiverGain);

#endif
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
