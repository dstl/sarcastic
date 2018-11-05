/***************************************************************************
 * 
 *           Module :  GeoConst.h
 *          Program :  sarclib
 *       Created by :  Darren Muff on 29/01/2012
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      Header file for Earth related constants
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

#ifndef sarclib_DEMdump_GeoConsts_h
#define sarclib_DEMdump_GeoConsts_h
#include <math.h>

#define GeoConsts_NODATAVALUE   (-32767.0)
#define GeoConsts_Re            (6378137.0)             ///< WGS84(G873) Specification
#define GeoConsts_Pff           (298.257223563)			///< WGS84(G873) Specification
#define GeoConsts_GM            (0.3986004418e15)		///< units are m^3 s^-2.
#define GeoConsts_zetaZ         (-0.53)					///< 'Zero degree term'
#define GeoConsts_Geqt          (9.7803253359)			///< Gravitation accel at equator in m/sec^2
#define GeoConsts_RADTODEG      (57.295779513082325)    ///< Convert radians to degrees
#define GeoConsts_DEGTORAD      (0.017453292519943)     ///< Convert degrees to radians
#define GeoConsts_a             (GeoConsts_Re)
#define GeoConsts_f             (1./GeoConsts_Pff)
#define GeoConsts_e             (sqrt(GeoConsts_f * (2.0-GeoConsts_f)))
#define GeoConsts_e2            (GeoConsts_e*GeoConsts_e)
#define GeoConsts_ep            (sqrt(GeoConsts_e2 / (1.0-GeoConsts_e2) ))
#define GeoConsts_ep2           (GeoConsts_ep*GeoConsts_ep)
#define GeoConsts_n             (GeoConsts_f / (2.0-GeoConsts_f))
#define GeoConsts_n2            (GeoConsts_n * GeoConsts_n)
#define GeoConsts_n3            (GeoConsts_n2 * GeoConsts_n)
#define GeoConsts_n4            (GeoConsts_n3 * GeoConsts_n)
#define GeoConsts_n5            (GeoConsts_n4 * GeoConsts_n)
#define GeoConsts_Ap            (GeoConsts_a*(1 - GeoConsts_n + (1.25*(GeoConsts_n2-GeoConsts_n3)) + (81./64.*(GeoConsts_n4-GeoConsts_n5))))
#define GeoConsts_Bp            (1.5*GeoConsts_a*(GeoConsts_n - GeoConsts_n2 + (0.875*(GeoConsts_n3-GeoConsts_n4)) + (55./64.*GeoConsts_n5)))
#define GeoConsts_Cp            (0.9375*GeoConsts_a*(GeoConsts_n2 - GeoConsts_n3 + 0.75*(GeoConsts_n4-GeoConsts_n5)))
#define GeoConsts_Dp            (35./48.*GeoConsts_a*(GeoConsts_n3 - GeoConsts_n4 + 11./16.*GeoConsts_n5))
#define GeoConsts_Ep            (315./512.*GeoConsts_a*(GeoConsts_n4-GeoConsts_n5))
#define GeoConsts_lat0          (0.0)
#define GeoConsts_k0            (0.9996)
#define GeoConsts_k02           (GeoConsts_k0  * GeoConsts_k0)
#define GeoConsts_k03           (GeoConsts_k02 * GeoConsts_k0)
#define GeoConsts_k04           (GeoConsts_k03 * GeoConsts_k0)
#define GeoConsts_k05           (GeoConsts_k04 * GeoConsts_k0)
#define GeoConsts_k06           (GeoConsts_k05 * GeoConsts_k0)
#define GeoConsts_k07           (GeoConsts_k06 * GeoConsts_k0)
#define GeoConsts_k08           (GeoConsts_k07 * GeoConsts_k0)
#define GeoConsts_S             (GeoConsts_Ap*GeoConsts_lat0 - GeoConsts_Bp*sin(2.0*GeoConsts_lat0)                     \
                                    + GeoConsts_Cp*sin(4.0*GeoConsts_lat0) - GeoConsts_Dp*sin(6.0*GeoConsts_lat0)       \
                                    + GeoConsts_Ep*sin(8.0*GeoConsts_lat0))
#define GeoConsts_p             (GeoConsts_a*(1-GeoConsts_e2) / pow( (1 - GeoConsts_e2 * pow(sin(GeoConsts_lat0), 2)), 1.5))
#define GeoConsts_ep4           (GeoConsts_ep2 * GeoConsts_ep2)
#define GeoConsts_ep6           (GeoConsts_ep4 * GeoConsts_ep2)
#define GeoConsts_ep8           (GeoConsts_ep4 * GeoConsts_ep4)
#define GeoConsts_k02_2         (GeoConsts_k02 * 2.0)
#define GeoConsts_ep4_4         (GeoConsts_ep4 * 4.0)
#define GeoConsts_ep2_9         (GeoConsts_ep2 * 9.0)
#define GeoConsts_k04_24        (GeoConsts_k04 * 24.0)
#define GeoConsts_ep2_46        (GeoConsts_ep2 * 46.0)
#define GeoConsts_ep2_252       (GeoConsts_ep2 * 252.0)
#define GeoConsts_ep4_n3        (GeoConsts_ep4 * (-3.0))
#define GeoConsts_ep6_100       (GeoConsts_ep6 * 100.0)
#define GeoConsts_ep4_66        (GeoConsts_ep4 * 66.0)
#define GeoConsts_ep2_n90       (GeoConsts_ep2 * (-90.0))
#define GeoConsts_ep8_88        (GeoConsts_ep8 * 88.0)
#define GeoConsts_k06_720       (GeoConsts_k06 * 720.0)
#define GeoConsts_k08_40320     (GeoConsts_k08 * 40320.0)
#define GeoConsts_k03_6         (GeoConsts_k03 * 6.0)
#define GeoConsts_ep2_6         (GeoConsts_ep2 * 6.0)
#define GeoConsts_ep4_3         (GeoConsts_ep4 * 3.0)
#define GeoConsts_ep2_8         (GeoConsts_ep2 * 8.0)
#define GeoConsts_ep6_4         (GeoConsts_ep6 * 4.0)
#define GeoConsts_ep6_24        (GeoConsts_ep6 * 24.0)
#define GeoConsts_k05_120       (GeoConsts_k05 * 120.0)
#define GeoConsts_k07_5040      (GeoConsts_k07 * 5040.0)
#define GeoConsts_FE            (500000)

#endif
