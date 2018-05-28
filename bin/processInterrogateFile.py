#!/usr/bin/env python
# /***************************************************************************
# *
# *       Module:    ProcessInterrogateFile.py
# *      Program:    ProcessInterrogateFile.py
# *   Created by:    Darren Muff on 12/07/2017.
# *                  Copyright (c) 2017 Dstl. All rights reserved.
# *
# *   Description:
# *     This script is used to extract information about individual rays from 
# *     an interrogate.txt file produced by SARCASTIC.
# *     The workflow goes like this:
# *      1) Generate a simulated cphd file using SARCASTIC
# *      2) Process the file into an image using a SAR image formation processor (eg tdpocl)
# *      3) Look at the image (eg using the dstl .dat file viewer)
# *      4) wonder what a blob is in the image
# *      5) find the location of the blob in x and y in the image
# *      6) run SARCASTIC again but for only one pulse. This time output 
# *         bounce 1 and 'Y' for interrogate
# *      7) enter your x & y coords of the weird artefact and the name/place to store 
# *         your interrogate file
# *      8) Load up Grapher on OSX (or some similar package that can input x,y,z coords)
# *      9) Copy and paste the sarcastic bounce info that came up in the terminal. You 
# *         should see a 3D model of the scene being illuminated
# *     10) run this script on the interrogation file produced by SARCASTIC. enter the 
# *         bounce info you are interested in seeing for the point you are investigating
# *     11) Plot the output on Grapher. This should be in teh form of rays
# *
# *   CLASSIFICATION        :  OFFICIAL
# *   Date of CLASSN        :  12/07/2017
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

import os
import os.path as osp
import argparse

class Point( object ):
    def __init__( self, x, y, z, data ):
        self.x, self.y, self.z = x, y, z
        self.data = data


parser = argparse.ArgumentParser(description="Program to parse an interrogate.txt file produced by SARCASTIC so that rays that cause a feature in an image can be individually identified using a plotting tool (like Grapher on OSX)")
parser.add_argument("filename", help="The name of an \'interrogate.txt\' file generated by SARCASTIC")
parser.add_argument("-b", "--bounce", help="The number of the bounce to process", type=int)
args = parser.parse_args()

question = -1
if args.bounce :
	question = args.bounce

print 'input file is ', args.filename

file = open(args.filename,'r')
for i in range(6):
	line = file.readline()
	print line[:-1]

part = line.split(":")[1]
TxPos = Point(float(part.split(',')[0]),float(part.split(',')[1]),float(part.split(',')[2]),"TxPos")

# Search through the file to find the approximate scattered energy back to the receiver for
# each bounce. (Each bounce in the interrogate file is for a point in the image)
#
EsDict = dict()
while True:
	line = file.readline()
	if not line: break
	if line.split()[0] == 'Bounce':
		ibounce = int(line.split(':')[1])
		line = file.readline()
		E = float(line.split(':')[1])
		EsDict[ibounce] = E


Es = sorted(EsDict, key=EsDict.__getitem__, reverse=True)
maxBounce = max(k for k,v in EsDict.iteritems())

print 'Summary of energy returned for each bounce for the interrogated image point (v/m^2)'
for i in Es:
	print '\tBounce ',i,' : ',EsDict[i]

if question < 0 or question > maxBounce:
	while True:
	    try:
	    	question = int(raw_input('Which bounce would you like to see ?'))
	    	if question > -1 and question <= maxBounce:
		        break
	    except:
	        print("That's not a valid option!")
	
print 'Output for bounce ',question
file.seek(0)
ind = 0
while True:
	line = file.readline()
	if not line: break
	if line.split()[0] == 'Bounce':
		ibounce = int(line.split(':')[1])
		# bounce = Es[ind]
		bounce = question
		if ibounce == bounce :
			line = file.readline()
			line = file.readline()
			ninters = int(line.split(':')[1])
			line = file.readline()
			for i in range(ninters) :
				line = file.readline()
				print TxPos.x,',',TxPos.y,',',TxPos.z
				for j in range(bounce+1) :
					coord = line.split()[3+j]
					pt = Point(float(coord.split(',')[0]),float(coord.split(',')[1]),float(coord.split(',')[2]),"point")
					print pt.x,',',pt.y,',',pt.z
				for j in range(bounce,-1,-1):
					coord = line.split()[3+j]
					pt = Point(float(coord.split(',')[0]),float(coord.split(',')[1]),float(coord.split(',')[2]),"point")
					print pt.x,',',pt.y,',',pt.z
				print TxPos.x,',',TxPos.y,',',TxPos.z
			ind = ind + 1

