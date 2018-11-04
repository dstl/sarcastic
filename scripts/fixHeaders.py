#!/usr/bin/env python3

import argparse
import re
import subprocess
import datetime
import os

def processfile(ip, prog, classn):
	regexpdesc = re.compile(r'description\s*:[a-z\s]*',re.IGNORECASE)
	regexpblank = re.compile(r'^\s*\*\s*$',re.IGNORECASE)
	regexpcreated = re.compile(r'created\s*by\s*.*\s([0-9]*/.*/[0-9]*).*',re.IGNORECASE)
	desc = ""
	date = datetime.datetime.today().strftime('%d-%b-%Y')
	ps = subprocess.Popen(('git', 'log', '--format=%ad',ip), stdout=subprocess.PIPE, cwd=currdir)
	createdate = subprocess.check_output(('tail', '-1'), stdin=ps.stdout).decode().rstrip()
	ps.wait()
	headerfound = False
	createdLine = ""
	now = datetime.datetime.now()
	f = open(ip)
	line = f.readline()
	while line:
		if line.startswith('/**') :
			headerfound = True
			while line and not line.endswith('*/\n'):
				line = f.readline()
				if regexpdesc.search(line) :
					# Found the description now keep reading lines until two blank ones are found
					blankcnt=0
					while line and blankcnt != 2 :
						desc = desc + line
						line = f.readline()
						if regexpblank.search(line) :
							blankcnt = blankcnt + 1
						else:
							blankcnt = 0

				if regexpcreated.search(line) :
					createdate = regexpcreated.search(line).group(1)

			print('/***************************************************************************')
			print(' * ')
			print(' *           Module : ',ip)
			print(' *          Program : ',prog)
			print(' *       Created by :  Darren Muff on',createdate)
			print(' *   CLASSIFICATION : ',classn)
			print(' *   Date of CLASSN : ',now.strftime("%d-%b-%Y"))
			print(desc.rstrip())
			print(' * ')
			print(' *   (c) Crown Copyright 2018 Defence Science and Technology Laboratory')
			print(' * ')
			print(' * Permission is hereby granted, free of charge, to any person obtaining a')
			print(' * copy of this software and associated documentation files (the "Software")')
			print(' * to deal in the Software without restriction, including without limitation')
			print(' * the rights to use, copy, modify, merge, publish, distribute, sublicense,')
			print(' * and/or sell copies of the Software, and to permit persons to whom the')
			print(' * Software is furnished to do so, subject to the following conditions:')
			print(' *')
			print(' * The above copyright notice and this permission notice shall be included')
			print(' * in all copies or substantial portions of the Software.')
			print(' * ')
			print(' * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS')
			print(' * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,')
			print(' * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL')
			print(' * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER')
			print(' * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING')
			print(' * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER')
			print(' * DEALINGS IN THE SOFTWARE.')
			print(' * ')

	
		if not headerfound :
			print('ERROR: file',ip,'does not have a header')
			exit(1)
		print(line, end='')
		line = f.readline()

	f.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Fix up the header info of a source file")
	parser.add_argument("-c", "--classification", help="Classification of module. [Official]")
	parser.add_argument("-p", "--program", help="Program name for module [current dirname]")
	parser.add_argument("filename", help="The name of the input source file")

	args = parser.parse_args()
	currdir = os.getcwd()

	ip   = args.filename
	if args.program :
		prog = args.program
	else:
		prog = os.path.basename(currdir)

	if args.classification :
		classn = args.classification
	else :
		classn ="Official"

	processfile(ip,prog,classn)

