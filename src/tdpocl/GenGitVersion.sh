#!/bin/bash
#***************************************************************************
#
#       Module:    GenGitVersion
#      Program:    TDPOCL
#   Created by:    Darren Muff on 19/07/2013.
#
#    (c) Crown Copyright 2018 Defence Science and Technology Laboratory
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
#   CLASSIFICATION        :  OFFICIAL
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
#**************************************************************************

# Get git tag and hash in the FullVersion
FULL_VERSION=$(git describe --dirty | sed -e 's/^v//' -e 's/^V//' -e 's/g//')

# Use the latest tag for short version
# Tags should be created as annotaed tags rather than lightweight tags (ie git tag -a <tagname>)
SHORT_VERSION=$(git describe --abbrev=0 | sed -e 's/^v//' -e 's/^V//')

# In order to create a monotonically increasing number for versions rather
# than a hash code I've included REVISION as well.
REVISION=$(git rev-list master | wc -l | awk '{print $1}')

# Also print the version date
VERSION_DATE=$(git show --format="%ci" | head -1)

echo "#ifdef REVISION"                             > tdpoclVersion.h
echo "#undef REVISION"                             >>tdpoclVersion.h
echo "#endif"                                      >>tdpoclVersion.h
echo "#ifdef SHORT_VERSION"                        >>tdpoclVersion.h
echo "#undef SHORT_VERSION"                        >>tdpoclVersion.h
echo "#endif"                                      >>tdpoclVersion.h
echo "#ifdef FULL_VERSION"                         >>tdpoclVersion.h
echo "#undef FULL_VERSION"                         >>tdpoclVersion.h
echo "#endif"                                      >>tdpoclVersion.h
echo "#ifdef VERSION_DATE"                         >>tdpoclVersion.h
echo "#undef VERSION_DATE"                         >>tdpoclVersion.h
echo "#endif"                                      >>tdpoclVersion.h
echo "#define REVISION \"${REVISION}\""            >>tdpoclVersion.h
echo "#define SHORT_VERSION \"${SHORT_VERSION}\""  >>tdpoclVersion.h
echo "#define FULL_VERSION \"${FULL_VERSION}\""    >>tdpoclVersion.h
echo "#define VERSION_DATE \"${VERSION_DATE}\""    >>tdpoclVersion.h
echo "${FULL_VERSION}"
