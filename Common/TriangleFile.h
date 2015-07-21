//
//  TriangleFile.h
//  sarcastic
//
//  Created by Darren Muff on 21/07/2015.
//  Copyright (c) 2015 Dstl. All rights reserved.
//

// Usage : ostream << TriangleFile

#ifndef __sarcastic__TriangleFile__
#define __sarcastic__TriangleFile__

#include <stdio.h>
#include <iostream>
#include <vector>
#include <SIlib/SIlib.h>
#include "trianglecpp.h"


class TriangleFile {
    static std::string version ;
    std::vector<Triangle> triangles ;

public:
    TriangleFile( std::vector<Triangle> tris ) ;
    void WriteFile( std::string filename ) ;
    void WritePLYFile( std::string filename, bool binary) ;
//    friend std::ostream & operator<<(std::ostream & os, const TriangleFile & triFile);
    


};

#endif /* defined(__sarcastic__TriangleFile__) */
