#!/bin/bash
#***************************************************************************
#
#       Module:    GenGitVersion
#      Program:    SIlib2
#   Created by:    Darren Muff on 19/07/2013.
#                  Copyright (c) 2013 [dstl]. All rights reserved.
#
#   Description:
#       This script automatically sets the revision, full and short version
#       string information into the header file Version.h.
#       It should be called by automatically before C file compilation
#       (or Version.h will not exist)
#       It uses git tags which must be created with the -a option (git tag -a <tagname>)
#       and be in the format v(or V)x.x.x 
#
#
#   CLASSIFICATION        :  UNCLASSIFIED
#   Date of CLASSN        :  19/07/2013
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
#
# THE SOFTWARE IN ITS ENTIRETY OR ANY SUBSTANTIAL PORTION SHOULD NOT BE
# USED AS A WHOLE OR COMPONENT PART OF DERIVATIVE SOFTWARE THAT IS BEING
# SOLD TO THE GOVERNMENT OF THE UNITED KINGDOM OF GREAT BRITAIN AND NORTHERN
# IRELAND.
#
#**************************************************************************

# Get git tag and hash in the FulLVersion
FULL_VERSION=$(git describe --dirty | sed -e 's/^v//' -e 's/^V//' -e 's/g//')

# Use the latest tag for short version
# Tags should be created as annotaed tags rather than lightweight tags (ie git tag -a <tagname>)
SHORT_VERSION=$(git describe --abbrev=0 | sed -e 's/^v//' -e 's/^V//')

# In order to create a monotonically increasing number for versions rather
# than a hash code I've included REVISION as well.
REVISION=$(git rev-list master | wc -l | awk '{print $1}')

echo "#define REVISION \"${REVISION}\""            > Version.h.tmp
echo "#define SHORT_VERSION \"${SHORT_VERSION}\""  >>Version.h.tmp
echo "#define FULL_VERSION \"${FULL_VERSION}\""    >>Version.h.tmp
echo "${FULL_VERSION}"

# Only copy Version.h.tmp on top of Version.h if they are not the same - this is
# to ensure that the timestamp on Version.h only gets modified when it really 
# needs to be, to ensure Make works correctly.
cmp -s Version.h Version.h.tmp || cp Version.h.tmp Version.h

# Remove our temp file
rm Version.h.tmp

