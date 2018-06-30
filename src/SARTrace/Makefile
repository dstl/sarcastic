# ***************************************************************************
# *
# *       Module:    Makefile
# *      Program:    SARcast
# *   Created by:    Darren Muff 15/04/2015
# *                  Copyright (c) 2015 [dstl]. All rights reserved.
# *
# *   Description:
# *      Makefile
# *
# *   CLASSIFICATION        :  UNCLASSIFIED
# *   Date of CLASSN        :  15/04/2015
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
EXE_NAME        = SARTrace
INSTALL_DIR     = /Users/darren
PROJ_DIR        = /Users/darren/Development/sarcastic
VERSION        := $(shell `pwd`/GenGitVersion.sh)
OCLKERNELSPATH := $(PROJ_DIR)/src/kernels

CPPSRC = \
../fastKdTree/AABB.cpp \
../fastKdTree/clipToAABB.cpp \
../fastKdTree/rayTrace.cpp \
../fastKdTree/kdTreeNode.cpp \
../fastKdTree/splitCandidate.cpp \
../Common/tinyxml/tinystr.cpp \
../Common/tinyxml/tinyxml.cpp \
../Common/tinyxml/tinyxmlerror.cpp \
../Common/tinyxml/tinyxmlparser.cpp \
../fastKdTree/buildTree.cpp \
../colladaToPlyFile/colladainterface.cpp \
../Common/TriangleMesh.cpp \
../sarcastic/buildRays.cpp \
../sarcastic/accelerateTriangles.cpp \
../sarcastic/buildRopesAndBoxes.cpp \
../bircs/buildKernel.cpp \
../sarcastic/reflect.cpp \
getSARTraceUserInput.cpp \
SARTrace.cpp \
SARTracebanner.cpp \
traceThreadCore.cpp \
../Common/readMaterialFile.cpp


SRC = \
../Common/matrixMultiplication.c \
../Common/OpenCLUtils.c \
../sarcastic/ecef2SceneCoords.c \
../Common/pmath_rng.c \
../Common/boxMullerRandom.c \
../Common/ranf.c


LIB = -lreadline -lsarclib

CPP  = clang -O3 -x c++ -std=gnu++11 -stdlib=libc++ -DOCLKERNELSPATH=$(OCLKERNELSPATH)
CC   = clang -O3 -DOCLKERNELSPATH=$(OCLKERNELSPATH)

OBJ    = $(patsubst %.c, %.o, $(SRC))
CPPOBJ = $(patsubst %.cpp, %.o, $(CPPSRC))

CFLAGS = -m64 -c `../scripts/sarclib-config --cflags` -I../sarcastic -I../Common -I/usr/local/include -I../fastKdTree -I../bircs -I../colladaToPlyFile
LFLAGS = -framework OpenCL -arch x86_64 `../scripts/sarclib-config --libs` -L/usr/local/dstl/lib -L/usr/local/lib

%.o:	%.cpp
		$(CPP) $< $(CFLAGS) -o $@

%.o:	%.c
		$(CC) $< $(CFLAGS) -o $@

$(EXE_NAME): $(OBJ) $(CPPOBJ)
		clang++ -stdlib=libc++ $(OBJ) $(CPPOBJ) $(LFLAGS) $(LIB) -o $(EXE_NAME)

clean:
		rm -fr $(OBJ) $(CPPOBJ)
		rm -fr $(EXE_NAME)

install:
		cp -p $(EXE_NAME) $(PROJ_DIR)/bin



