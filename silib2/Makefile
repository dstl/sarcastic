#/***************************************************************************
# *
# *       Module:    Makefile.in
# *      Program:    SIlib2
# *   Created by:    Darren Muff 21/07/2013
# *                  Copyright (c) 2016 [dstl]. All rights reserved.
# *
# *   Description:
# *      Makefile template
# *
# *   CLASSIFICATION        :  UNCLASSIFIED
# *   Date of CLASSN        :  19/07/2013
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

SRC = alloc.c dataio.c dataio_au4.c dataio_c8.c dataio_ce.c dataio_cphd.c dataio_cphd3.c dataio_cphdx.c dataio_envi.c dataio_gdal.c dataio_gdal_stubs.c dataio_hdf4.c dataio_im.c dataio_input.c dataio_sicd.c dataio_sio.c dataio_srf.c enum2string.c fft_c8.c GeoVector.c image.c image_average.c image_conversion.c image_extract.c image_fft.c image_fftw.c image_interp.c image_math.c image_matrix.c image_pad.c image_random.c image_recip.c image_stats.c image_transpose.c image_unwrap.c interpolate.c kaiser.c latlon.c timer.c vector.c collectionGeometry.c printCPHDCollectionInfo.c


INCLUDE = SIlib2.h alloc.h cmplx.h dataio.h dataio_au4.h dataio_c8.h dataio_ce.h dataio_cphd.h dataio_envi.h dataio_gdal.h dataio_hdf4.h dataio_im.h dataio_sicd.h dataio_sio.h dataio_srf.h enum2string.h error_defn.h fft.h GeoConst.h GeoVector.h image.h image_interp.h image_stats.h interpolate.h latlon.h math_macro.h matrix.h SIlib2_xml.h mtmaths.h  physical.h Timer.h vector.h window.h collectionGeometry.h Version.h


CC	:= /usr/bin/clang
CFLAGS = -O2 -Wall -c -g -D_REENTRANT -Wmissing-prototypes 

# The following line is expanded by the configure script.
have_gdal=Yes

ifeq ($(have_gdal), No)
	SRC += dataio_gdal_stubs.c
	HAVE_GDAL=NO
	GDAL_LIBS=
	GDAL_CFLAGS=
else
	SRC += dataio_gdal.c
	CFLAGS += -DSILIB2_HAVE_GDAL `gdal-config --cflags`
	HAVE_GDAL=YES
	GDAL_CFLAGS=-DSILIB2_HAVE_GDAL
	GDAL_LIBS := $(shell gdal-config --libs)
endif

OBJ = $(patsubst %.c, %.o, $(SRC))

VERSION := $(shell `pwd`/GenGitVersion.sh)
CFLAGS+=-DVERSION="\"$(VERSION)\""

uname :=  $(shell uname)

ifeq ($(uname),SunOS)
	CFLAGS += -m64 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DHAVE_READLINE=1 -I/usr/local/dstl/include
endif

ifeq ($(uname),Linux)
	CFLAGS += -m64 -std=gnu99 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DHAVE_READLINE=1 -I/usr/local/include
# Not need anymore?  -I/usr/local/cuda/include -I/opt/AMDAPP/include
endif

ifeq ($(uname),Darwin)
	CFLAGS += -arch x86_64 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DHAVE_READLINE=1 -I/usr/local/include
endif

lib:	libSIlib2.a

libSIlib2.a:	$(OBJ) $(INCLUDE)
		@echo Making library libSIlib2.a
		ar -r libSIlib2.a $(OBJ)

install:	lib
		mkdir -p /usr/local/dstl/lib
		cp libSIlib2.a /usr/local/dstl/lib/libSIlib2-$(VERSION).a
		rm -f /usr/local/dstl/lib/libSIlib2.a
		ln -s /usr/local/dstl/lib/libSIlib2-$(VERSION).a /usr/local/dstl/lib/libSIlib2.a

		mkdir -p /usr/local/dstl/include/SIlib2-$(VERSION)
		cp $(INCLUDE) /usr/local/dstl/include/SIlib2-$(VERSION)
		rm -f /usr/local/dstl/include/SIlib2
		ln -s /usr/local/dstl/include/SIlib2-$(VERSION) /usr/local/dstl/include/SIlib2
		rm -f SIlib2-config
		@echo >SIlib2-config "#!/bin/sh"
		@echo >>SIlib2-config "CONFIG_LIBS=\"-L/usr/local/dstl/lib -lSIlib2 $(GDAL_LIBS) -lreadline -lexpat -lfftw3f -lfftw3f_threads\""
		@echo >>SIlib2-config "CONFIG_CFLAGS=\"-I/usr/local/dstl/include $(GDAL_CFLAGS)\""
		@echo >>SIlib2-config "CONFIG_VERSION=\"$(VERSION)\""
		@echo >>SIlib2-config "CONFIG_PREFIX=\"/usr/local/dstl\""
		@echo >>SIlib2-config "CONFIG_GDAL=\"$(HAVE_GDAL)\""
		cat >>SIlib2-config SIlib2-config.in
		cp  SIlib2-config /usr/local/dstl/bin
		chmod +x /usr/local/dstl/bin/SIlib2-config

%.o:	%.c $(INCLUDE)
	$(CC) $< $(CFLAGS) -o $@ 

fft_c8.o:	fft_c8.c
	$(CC) $< $(CFLAGS) -DTRUNCATE -o $@

clean:	
	rm -f *.o

realclean:
	rm -f *.o libSIlib2.a

docs:
	doxygen doc/Doxyfile

