/***************************************************************************
 *  materialProperties.h
 *  triDecomp
 *
 *  Created by Muff Darren on 25/09/2014.
 *  Copyright (c) 2014 [dstl]. All rights reserved.
 *
 *   Description:
 *      Header file defining different radar material properties
 *
 *  CLASSIFICATION       :   UNCLASSIFIED
 *  Date of CLASSN       :   25/09/2014
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
 ***************************************************************************/

#ifndef sarcastic_materialProperties_h
#define sarcastic_materialProperties_h
#define MATBYTES 128

typedef struct scatProps {
    char   matname[MATBYTES] ;
    float  corlen      ;
    float  roughness   ;
    float  resistivity ;
    float  specular    ;
    float  diffuse     ;
    float  shinyness   ;
} scatProps ;

#define NMATERIALS 9
scatProps materialProperties[NMATERIALS] = {
//   Name          corrLen      Roughness   Resistivity Specular    Diffuse     Shinyness
    {"ASPHALT",     0.5,        0.005,      1.0e18,     0.8,        0.2,        30.0        },
    {"BRICK",       0.1,        0.001,      1.0e18,     0.7,        0.3,        20.0        },
    {"CONCRETE",    0.2,        0.01,       120.0,      0.3,        0.7,        10.0        },
    {"METAL",       100.0,      0.0,        1.0e-8,     1.0,        0.0,        50.0        },
    {"MATERIAL",    100.0,      0.0,        0.0,        1.0,        0.0,        50.0        },
    {"ROOFING",     0.1,        0.1,        1.0e18,     0.6,        0.4,        40.0        },
    {"VEGETATION",  0.01,       0.1,        2000.0,     0.2,        0.8,        5.0         },
    {"WATER",       0.01,       0.1,        2.0e1,      1.0,        0.0,        50.0        },
    {"WOOD",        0.1,        0.001,      1.0e14,     0.6,        0.4,        10.0        }
} ;


#endif
