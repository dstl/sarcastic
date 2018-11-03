/***************************************************************************
 * 
 *           Module :  colladainterface.h
 *          Program :  colladaToPlyFile
 *       Created by :  Darren Muff on Wed Jul 11 08:12:36 2018 -0400
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *      Description :  A basic C++ class for reading collada (.dae) files
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
#ifndef COLLADAINTERFACE_H
#define COLLADAINTERFACE_H

#include <iostream>
#include <vector>
#include <map>
#include <sstream>
#include <iterator>

#include "tinyxml/tinyxml.h"

// The following defines are taken from <OpenGL/gl3.h> and
// are used here is collada is needed without using OpenGL
//
#ifndef __gl3_h_
typedef uint32_t GLenum;
/* BeginMode */
#define GL_POINTS                         0x0000
#define GL_LINES                          0x0001
#define GL_LINE_LOOP                      0x0002
#define GL_LINE_STRIP                     0x0003
#define GL_TRIANGLES                      0x0004
#define GL_TRIANGLE_STRIP                 0x0005
#define GL_TRIANGLE_FAN                   0x0006
/* DataType */
#define GL_BYTE                           0x1400
#define GL_UNSIGNED_BYTE                  0x1401
#define GL_SHORT                          0x1402
#define GL_UNSIGNED_SHORT                 0x1403
#define GL_INT                            0x1404
#define GL_UNSIGNED_INT                   0x1405
#define GL_FLOAT                          0x1406
#define GL_DOUBLE                         0x140A
#endif

struct SourceData {
    GLenum type;
    unsigned int size;
    unsigned int stride;
    void* data;
};

typedef std::map<std::string, SourceData> SourceMap;

struct ColGeom {
    std::string name ;
    SourceMap map ;
    GLenum primitive ;
    int index_count ;
    unsigned short* indices ;
    float scaling ;
    double transform[16];
    std::string materialSide1 = std::string("") ;
    std::string materialSide2 = std::string("") ;
};

SourceData readSource(TiXmlElement*);

class ColladaInterface {
    
public:
    ColladaInterface() {};
    static void readGeometries(std::vector<ColGeom>*, const char*);
    static void freeGeometries(std::vector<ColGeom>*);
};

#endif

