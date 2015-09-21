# ***************************************************************************
# *
# *       Module:    Makefile
# *      Program:    sarcastic
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
UNAME      := $(shell uname)
VERSION    := $(shell `pwd`/GenGitVersion.sh)
INSTALLDIR := /usr/local/dstl
CURRDIR    := $(shell pwd)
OCLKERNELSPATH = $(PROJECT_DIR)/sarcastic/kernels

CPPSRC =


SRC = \
banner.c \
buildKernel.c \
BuildRopesAndBoxes.c \
ecef2SceneCoords.c \
generate_gaussian_ray_distribution.c \
getUserInput.c \
oclBuildShadowRays.c \
oclKdTreeHits.c \
oclPOField.c \
oclRandomRays.c \
oclReflect.c \
oclReflectionPower.c \
matrixMultiplication.c \
OpenCLUtils.c \
RCS.c \
readKdTree.c \
sarcastic.c \
TxPowerPerRay.c \
threadCore.c


LIB = -lreadline -lSIlib -lfftw3f_threads -lfftw3f -lfftw3_threads -lfftw3 -framework OpenCL -lexpat
CPP  = clang -x objective-c -fobjc-arc
CC   = clang

OBJ    = $(patsubst %.c, objs/%.o, $(SRC))
CPPOBJ = $(patsubst %.m, objs/%.o, $(CPPSRC))

CFLAGS = -m64 -c `SIlib-config --cflags` -I../sarcastic -I../Common -DOCLKERNELSPATH=\"$(OCLKERNELSPATH)\"
LFLAGS = -arch x86_64 -L/usr/local/dstl/lib

objs/%.o:	%.m
		$(CPP) $< $(CFLAGS) -o $@

objs/%.o:	%.c
		$(CC) $< $(CFLAGS) -o $@

sarcastic: $(OBJ) $(CPPOBJ)
	clang  $(OBJ) $(CPPOBJ) $(LFLAGS) $(LIB) -o ../bin/sarcastic

clean:
	rm -fr objs/*.o

install:
	cp -p ../bin/sarcastic $(INSTALL_DIR)

	


