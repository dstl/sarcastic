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
    float transform[16];
    std::string material ;
};

SourceData readSource(TiXmlElement*);

class ColladaInterface {
    
public:
    ColladaInterface() {};
    static void readGeometries(std::vector<ColGeom>*, const char*);
    static void freeGeometries(std::vector<ColGeom>*);
};

#endif

