# /***************************************************************************
# *
# *       Module:    Makefile
# *      Program:    sarcastic2
# *   Created by:    Darren Muff on 11/06/2017.
# *                  Copyright (c) 2017 Dstl. All rights reserved.
# *
# *   Description:
# *     Makefile for sarcastic
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
EXE_NAME = sarcastic
INSTALL_DIR = /usr/local/dstl/bin
VERSION    := $(shell `pwd`/GenGitVersion.sh)

CPPSRC = \
../fastKdTree/AABB.cpp \
banner.cpp \
buildRays.cpp \
../fastKdTree/buildTree.cpp \
../fastKdTree/clipToAABB.cpp \
../colladaToPlyFile/colladainterface.cpp \
../fastKdTree/kdTreeNode.cpp \
../fastKdTree/rayTrace.cpp \
reflect.cpp \
../fastKdTree/splitCandidate.cpp \
../Common/tinyxml/tinystr.cpp \
../Common/tinyxml/tinyxml.cpp \
../Common/tinyxml/tinyxmlerror.cpp \
../Common/tinyxml/tinyxmlparser.cpp \
../Common/TriangleMesh.cpp \
getUserInput.cpp \
shadowRays.cpp \
accelerateTriangles.cpp \
cpuPOField.cpp \
buildRopesAndBoxes.cpp \
sarcastic2.cpp \
tryReadFile.cpp \
threadCore.cpp \
../Common/readMaterialFile.cpp

SRC = \
../Common/boxMullerRandom.c \
../Common/printProgress.c \
../sarcastic/ecef2SceneCoords.c \
../Common/matrixMultiplication.c \
../Common/ranf.c \
../sarcastic/RCS.c \
../sarcastic/TxPowerPerRay.c \
../Common/pmath_rng.c


LIB = -lreadline -lSIlib2


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

	


