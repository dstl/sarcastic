/***************************************************************************
 * 
 *           Module :  materialProperties.h
 *          Program :  Sarcastic
 *       Created by :  Darren Muff on 25/09/2014
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *   Description:
 *      Header file defining different radar material properties
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

#ifndef sarcastic_materialProperties_h
#define sarcastic_materialProperties_h
#define MATBYTES 128

typedef struct scatProps {
    char   matname[MATBYTES] ;
    float  corlen      ;
    float  roughness   ;
    float  Rs          ;
    float  Rm          ;
    float  specular    ;
    float  diffuse     ;
    float  shinyness   ;
} scatProps ;
#define NMATERIALS 10
extern int         globalNumMats ;
extern scatProps   *globalMatProps;
extern int         *globalMatColours;


// Resistivities taken from http://chemistry.about.com/od/moleculescompounds/a/Table-Of-Electrical-Resistivity-And-Conductivity.htm
//
static scatProps materialProperties[NMATERIALS] = {
//   Name           corrLen     Roughness   Rs          Rm     Specular    Diffuse     Shinyness
    {"MATERIAL",    0.05,        0.0,        0.0,        9e9,  0.7,        0.0,        50.0        },
    {"ASPHALT",     0.5,        0.05,       1.0e18,     9e9,   0.8,        0.2,        30.0        },
    {"BRICK",       0.1,        0.001,      1.0e18,     9e9,   0.7,        0.3,        20.0        },
    {"CONCRETE",    0.2,        0.01,       120.0,      9e9,   0.3,        0.7,        10.0        },
    {"METAL",       0.6,        0.0,        1.0e-8,     9e9,   1.0,        0.0,        50.0        },
    {"ROOFING",     0.1,        0.1,        1.0e18,     9e9,   0.6,        0.4,        40.0        },
    {"VEGETATION",  0.1,        0.1,        2000.0,     9e9,   0.2,        0.8,        5.0         },
    {"WATER",       0.01,       0.1,        2.0e1,      9e9,   1.0,        0.0,        50.0        },
    {"WOOD",        0.1,        0.001,      1.0e14,     9e9,   0.6,        0.4,        10.0        },
    {"GRASS",       0.05,       0.1,        2000.0,     9e9,   0.2,        0.8,        5.0         }
} ;

static int materialColours[NMATERIALS][3] = {
    {255, 255, 255},        // 0 = Material
    {128, 128, 128},        // 1 = Asphalt
    {224, 224, 224},        // 2 = Brick
    {176,  19,  35},        // 3 = Concrete
    {176,  94,  41},        // 4 = Metal
    {214, 186, 062},        // 5 = Roofing
    {166, 214, 054},        // 6 = Vegetation
    {056, 125, 214},        // 7 = Water
    {137,  46, 014},        // 8 = Wood
    {100, 214, 050},        // 9 = Grass
};


#endif
