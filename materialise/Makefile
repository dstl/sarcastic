# /***************************************************************************
# *
# *       Module:    Makefile
# *      Program:    materialise
# *   Created by:    Darren Muff on 11/06/2017.
# *                  Copyright (c) 2017 Dstl. All rights reserved.
# *
# *   Description:
# *
# *      This programme reads in a .ply file and breaks it down into Delaunay
# *      triangles. The input .ply file can be a mesh or a triangle soup. In the
# *      latter case (where vertex points can be unique to a triangle) the
# *      programme searches for them and splits the soup up so that each triangle
# *      has only onle pertner for each of its edges. The final trianglation has
# *      a roughness applied that is determined by the input 'material'
# *      characteristics and the surface roughness parameters specified in the
# *      materialProperties.h header file. if any of the surfaces of the model
# *      has a surface roughness that is not zero then the triangulation will
# *      perform a post process pass to randomly 'roughen' the surface such that
# *      the standard deviation is equivalent to the specified surface roughness.
# *
# *      The progremme makes extensive use of teh CGAL library [1] and in
# *      particular the Delaunay triangulation routines [2]
# *
# *          [1] CGAL, Computational Geometry Algorithms Library, http://www.cgal.org
# *          [2] Laurent Rineau. 2D Conforming Triangulations and Meshes.
# *              In CGAL User and Reference Manual. CGAL Editorial Board, 4.10 edition, 2017
# *
# *
# *   CLASSIFICATION        :  OFFICIAL
# *   Date of CLASSN        :  11/06/2017
# *
# * Permission is hereby granted, free of charge, to any person obtaining a
# * copy of this software and associated documentation files (the "Software"),
# * to deal in the Software without restriction, including without limitation
# * the rights to use, copy, modify, merge, publish, distribute, sublicense,
# * and/or sell copies of the Software, and to permit persons to whom the
# * Software is furnished to do so, subject to the following conditions:
# *
# * The above copyright notice and this permission notice shall be included
# * in all copies or substantial portions of the Software.
# *
# * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# * DEALINGS IN THE SOFTWARE.
# *
# * THE SOFTWARE IN ITS ENTIRETY OR ANY SUBSTANTIAL PORTION SHOULD NOT BE
# * USED AS A WHOLE OR COMPONENT PART OF DERIVATIVE SOFTWARE THAT IS BEING
# * SOLD TO THE GOVERNMENT OF THE UNITED KINGDOM OF GREAT BRITAIN AND NORTHERN
# * IRELAND.
# *
# ***************************************************************************/

EXE_NAME = materialise
INSTALL_DIR = /usr/local/dstl/bin
VERSION    := $(shell `pwd`/GenGitVersion.sh)

CPPSRC = \
../colladaToPlyFile/colladainterface.cpp \
main.cpp \
../Common/tinyxml/tinystr.cpp \
../Common/tinyxml/tinyxml.cpp \
../Common/tinyxml/tinyxmlerror.cpp \
../Common/tinyxml/tinyxmlparser.cpp \
../Common/TriangleMesh.cpp \
../Common/readMaterialFile.cpp

SRC = \
../Common/boxMullerRandom.c \
../Common/matrixMultiplication.c \
../Common/pmath_rng.c \
../Common/printProgress.c \
../Common/ranf.c


LIB = -lreadline -lSIlib2 -lboost_filesystem -lboost_system -lboost_thread-mt -lgmp.10 -lmpfr.4 -lCGAL_Core -lCGAL
CPP  = clang -O3 -x c++ -std=gnu++11 -stdlib=libc++
CC   = clang

OBJ    = $(patsubst %.c, %.o, $(SRC))
CPPOBJ = $(patsubst %.cpp, %.o, $(CPPSRC))

CFLAGS = -m64 -c `SIlib2-config --cflags` -I../sarcastic -I../Common -I../Common/tinyxml -I../fastKdTree -I../sarcastic2 -I/usr/local/include -I../colladaToPlyFile
LFLAGS = -arch x86_64 `SIlib2-config --libs` -L/usr/local/dstl/lib -L/usr/local/lib 

%.o:	%.cpp
		$(CPP) $< $(CFLAGS) -o $@

%.o:	%.c
		$(CC) $< $(CFLAGS) -o $@

$(EXE_NAME): $(OBJ) $(CPPOBJ)
	mkdir -p ../bin
	clang++ -stdlib=libc++ $(OBJ) $(CPPOBJ) $(LFLAGS) $(LIB) -o ../bin/$(EXE_NAME)

clean:
	rm -fr $(OBJ) $(CPPOBJ)
	rm -fr ../bin/$(EXE_NAME)

install:
	cp -p ../bin/$(EXE_NAME) $(INSTALL_DIR)

	


