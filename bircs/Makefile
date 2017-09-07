# /***************************************************************************
# *
# *       Module:    Makefile
# *      Program:    bircs
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
EXE_NAME = bircs
INSTALL_DIR = /usr/local/dstl/bin
VERSION    := $(shell `pwd`/GenGitVersion.sh)
OCLKERNELSPATH = \"/users/darren/Dropbox/Development/sarcastic/bircs/kernels\"

CPPSRC = \
../sarcastic2/buildRays.cpp \
../fastKdTree/rayTrace.cpp \
../Common/tinyxml/tinyxmlparser.cpp \
../Common/tinyxml/tinystr.cpp \
../Common/tinyxml/tinyxml.cpp \
../Common/tinyxml/tinyxmlerror.cpp \
../colladaToTriFile/colladainterface.cpp \
../fastKdTree/splitCandidate.cpp \
../fastKdTree/AABB.cpp \
../fastKdTree/clipToAABB.cpp \
../fastKdTree/kdTreeNode.cpp \
../fastKdTree/buildTree.cpp \
../Common/TriangleMesh.cpp \
../sarcastic2/tryReadFile.cpp \
../sarcastic2/accelerateTriangles.cpp \
buildKernel.cpp \
../sarcastic2/buildRopesAndBoxes.cpp \
../sarcastic2/cpuPOField.cpp \
main.cpp \
threadCore.cpp \
../sarcastic2/shadowRays.cpp \
../sarcastic2/reflect.cpp

SRC = \
bircsBanner.c \
../sarcastic/TxPowerPerRay.c \
../Common/matrixMultiplication.c \
../Common/ranf.c \
../Common/boxMullerRandom.c \
../Common/pmath_rng.c \
../sarcastic/RCS.c \
../sarcastic/OpenCLUtils.c


LIB = -lreadline -lSIlib2


CPP  = clang -O3 -x c++ -std=gnu++11 -stdlib=libc++ -DOCLKERNELSPATH=$(OCLKERNELSPATH)
CC   = clang -O3 -DOCLKERNELSPATH=$(OCLKERNELSPATH)

OBJ    = $(patsubst %.c, %.o, $(SRC))
CPPOBJ = $(patsubst %.cpp, %.o, $(CPPSRC))

CFLAGS = -m64 -c `SIlib2-config --cflags` -I../sarcastic -I../Common -I../Common/tinyxml -I../fastKdTree -I../sarcastic2 -I/usr/local/include -I../colladaToTriFile
LFLAGS = -framework OpenCL -arch x86_64 `SIlib2-config --libs` -L/usr/local/dstl/lib -L/usr/local/lib

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

	


