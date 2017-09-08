# ***************************************************************************
# *
# *       Module:    Makefile
# *      Program:    colladaToTriFile
# *   Created by:    Darren Muff 15/04/2015
# *                  Copyright (c) 2015 [dstl]. All rights reserved.
# *
# *   Description:
# *      Makefile for colladaToTriFile
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

EXE         = colladaToPlyFile
INSTALL_DIR = /usr/local/dstl/bin
VERSION    := $(shell `pwd`/GenGitVersion.sh)

CPPSRC = \
colladainterface.cpp \
main.cpp \
banner.cpp \
../Common/tinyxml/tinystr.cpp \
../Common/tinyxml/tinyxml.cpp \
../Common/tinyxml/tinyxmlerror.cpp \
../Common/tinyxml/tinyxmlparser.cpp \
../Common/TriangleMesh.cpp \
../Common/readMaterialFile.cpp

SRC = ../Common/matrixMultiplication.c

LIB = -lreadline -lSIlib2
CPP  = clang -O3 -x c++ -std=gnu++11 -stdlib=libc++
CC   = clang -O3

OBJ    = $(patsubst %.c, %.o, $(SRC))
CPPOBJ = $(patsubst %.cpp, %.o, $(CPPSRC))

CFLAGS = -m64 -c `SIlib2-config --cflags` -I. -I../sarcastic -I../sarcastic2 -I../Common -I../Common/tinyxml -I../fastKdTree -I/usr/local/include
LFLAGS = -arch x86_64 -L/usr/local/dstl/lib `SIlib2-config --libs` -L/usr/local/lib 

%.o:	%.cpp
		$(CPP) $< $(CFLAGS) -o $@

%.o:	%.c
		$(CC) $< $(CFLAGS) -o $@

$(EXE): $(OBJ) $(CPPOBJ)
	clang++ -stdlib=libc++ $(OBJ) $(CPPOBJ) $(LFLAGS) $(LIB) -o ../bin/$(EXE)

clean:
	rm -fr $(OBJ) $(CPPOBJ) ../bin/$(EXE)

install:
	cp -p ../bin/$(EXE) $(INSTALL_DIR)

	


