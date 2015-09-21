# ***************************************************************************
# *
# *       Module:    Makefile
# *      Program:    materialise
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

PROJECT_DIR = /Users/darren/Development/sarcastic
INSTALL_DIR = /usr/local/dstl/bin

CPPSRC = \
main.cpp \
TriangleFile.cpp \
trianglecpp.cpp

SRC = \
boxMullerRandom.c \
JRStriangle.c \
matrixMultiplication.c \
pmath_rng.c \
ranf.c

LIB = -lreadline -lSIlib
CPP  = clang -x c++ -std=gnu++11 -stdlib=libc++
CC   = clang

OBJ    = $(patsubst %.c, objs/%.o, $(SRC))
CPPOBJ = $(patsubst %.cpp, objs/%.o, $(CPPSRC))

CFLAGS = -m64 -c `SIlib-config --cflags` -I../sarcastic -I../Common -I../Common/tinyxml
LFLAGS = -arch x86_64 -L/usr/local/dstl/lib

objs/%.o:	%.cpp
		$(CPP) $< $(CFLAGS) -o $@

objs/%.o:	%.c
		$(CC) $< $(CFLAGS) -o $@

materialise: $(OBJ) $(CPPOBJ)
	clang++ -stdlib=libc++ $(OBJ) $(CPPOBJ) $(LFLAGS) $(LIB) -o ../bin/materialise

clean:
	rm -fr objs/*.o

install:
	cp -p ../bin/materialise $(INSTALL_DIR)

	


