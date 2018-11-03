/***************************************************************************
 * 
 *           Module :  TriangleFile.h
 *          Program :  Sarcastic
 *       Created by :  Darren Muff on 21/07/2015
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *   Description:
 *      Class header for TriangleFile class
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


#ifndef __sarcastic__TriangleFile__
#define __sarcastic__TriangleFile__

#include <stdio.h>
#include <iostream>
#include <vector>
#include <sarclib/sarclib.h>
#include "trianglecpp.h"

class TriangleFile {
    static std::string version ;

public:
    std::vector<Triangle> triangles ;
    std::vector<triangleReference> triangleReferences ;
    std::vector<Point> triangleVertices ;

    TriangleFile( std::vector<Triangle> tris ) ;
    TriangleFile( std::string filename) ;
    void WriteFile( std::string filename ) ;
    void WritePLYFile( std::string filename, bool binary) ;
    void NewWritePLYFile( std::string filename, bool binary) ;
    void sortTrianglesAndPoints() ;

};

#endif /* defined(__sarcastic__TriangleFile__) */
