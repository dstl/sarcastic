#!/usr/bin/env python

# Use GDAL to read the NITF data
#
from osgeo import gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
from DstlDat import DstlDat

# and of course numpy and opencv
#
import numpy as np
import cv2

# import Argument parser
#
import argparse

# ALlow deepcopies for 2d arrays
#
from copy import copy, deepcopy

def loadDatFile(fname):
        f = DstlDat(fname)
        dat = f.load((0,f.ny,0,f.nx))
        return dat

def xoversample(amparr, currentOversample, xover):
	nx = amparr.shape[1]
	cutoff = int(floor(nx / currentOversample / 2))
	IM = np.fft.fft(amparr,axis=1)
	A = IM[:,range(nx-cutoff,nx)]
	B = IM[:,range(cutoff)]
	Z = zeros((amparr.shape[0],int(cutoff*2*xover - (cutoff*2))), dtype=complex )
	IM = np.concatenate([B,Z],axis=1)
	C = np.concatenate([IM,A], axis=1)
	D = np.fft.fftshift(C, axes=(1,))
	im = np.fft.ifft(D, axis=1)
	amparr = np.abs(im)
	return amparr

def yoversample(amparr, currentOversample, yover):
	amparr = np.transpose(amparr)
	amparr = xoversample(amparr, currentOversample, yover)
	amparr = np.transpose(amparr)
	return amparr

# define a simple function that averages down an image based on an input shape
# (we do this because we want to display the nitf image and the datasets are quite large)
#
def aveShape(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

# 
mouseX,mouseY = -1,-1 
def draw_circle(event,x,y,flags,param):
	global mouseX,mouseY
	if event == cv2.EVENT_LBUTTONUP:
		mouseX,mouseY = x,y

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Program to view complex NITF files")
	parser.add_argument("filename", help="The name of a complex NITF file display")
	parser.add_argument("-x", "--xcentre", help="The x coord of the centre of the display region", type=int)
	parser.add_argument("-y", "--ycentre", help="The y coord of the centre of the display region", type=int)
	parser.add_argument("-nx", "--xsize", help="The x dimension of the display region", type=int)
	parser.add_argument("-ny", "--ysize", help="The y dimension of the display region", type=int)
	parser.add_argument("-xs", "--xshrink", help="The x shrink factor", type=int)
	parser.add_argument("-ys", "--yshrink", help="The y shrink factor", type=int)
	parser.add_argument("-a", "--alpha", help="Contrast scaling factor out=(alpha*in)+beta. Default 1", type=float)
	parser.add_argument("-b", "--beta", help="Contrast scaling factor out=(alpha*in)+beta. Default 0", type=float)
	parser.add_argument("-p", "--phase", help="Display phase of image", action="store_true")
	parser.add_argument("-xo", "--xoversample", help="Amount of oversampling in data in X", type=float)
	parser.add_argument("-yo", "--yoversample", help="Amount of oversampling in data in Y", type=float)

	args = parser.parse_args() 

	if args.phase :
		dispPhase = True
	else:
		dispPhase = False

	fileName = args.filename
	if fileName.endswith('.dat'):
		cImg = loadDatFile(fileName)
		amplitudeArr = np.abs(cImg)
	elif fileName.endswith('.ntf'):
		dataset  = gdal.Open(fileName, GA_ReadOnly)
		amplitudeBand = dataset.GetRasterBand(1)
		phaseBand     = dataset.GetRasterBand(2)
		if dispPhase:
			amplitudeArr = BandReadAsArray(phaseBand)
		else:
			amplitudeArr = BandReadAsArray(amplitudeBand)
	else:
		print 'Unknown File Type'
		exit(1)
	
	if args.xsize :
		xsize = args.xsize
	elif amplitudeArr.shape[1] > 1024 :
		xsize = 1024
	else :
		xsize = amplitudeArr.shape[1]

	if args.ysize :
		ysize = args.ysize
	elif amplitudeArr.shape[0] > 1024 :
		ysize = 1024
	else :
		ysize = amplitudeArr.shape[0]

	if args.xcentre :
		xc = args.xcentre
		if xc < xsize / 2:
			xc = xsize / 2
	else :
		xc = amplitudeArr.shape[1] / 2

	if args.ycentre :
		yc = args.ycentre
		if yc < ysize / 2:
			yc = ysize / 2
	else :
		yc = amplitudeArr.shape[0] / 2

	if args.xshrink :
		xshrink = args.xshrink
	else :
		xshrink = 1

	if args.yshrink :
		yshrink = args.yshrink
	else :
		yshrink = 1

	if args.alpha :
		alp = args.alpha
	else:
		alp = 1

	if args.beta :
		bet = args.beta
	else :
		bet = 0

	xstart = xc - (xsize/2)
	if xstart < 0:
		xstart = 0
	xend = xc + (xsize/2)
	if xend > amplitudeArr.shape[1] :
		xend = amplitudeArr.shape[1]
	ystart = yc - (ysize/2)
	if ystart < 0 :
		ystart = 0
	yend = yc + (ysize/2)
	if yend > amplitudeArr.shape[0] :
		yend = amplitudeArr.shape[0]


	amparr = amplitudeArr[ystart:yend, xstart:xend]
	
	if dispPhase:
		amparr = amparr * 2 * np.pi / 65536

	if args.xoversample :
		amparr = xoversample(amparr, 1, args.xoversample)

	if args.yoversample :
		amparr = yoversample(amparr, 1, args.yoversample)

	# OSX = 1.0
	# nx = amparr.shape[1]
	# oversamp = 1.33
	# cutoff = int(floor(nx / oversamp / 2))

	# IM = np.fft.fft(amparr,axis=1)
	# A = IM[:,range(nx-cutoff,nx)]
	# B = IM[:,range(cutoff)]
	# Z = zeros((amparr.shape[0],int(cutoff*2*OSX - (cutoff*2))), dtype=complex )
	# IM = np.concatenate([B,Z],axis=1)
	# C = np.concatenate([IM,A], axis=1)
	# D = np.fft.fftshift(C, axes=(1,))
	# im = np.fft.ifft(D, axis=1)
	# amparr = np.abs(im)

	sh     = amparr.shape
	imgsm  = aveShape(amparr,(sh[0]//yshrink,sh[1]//xshrink)) 
	imguint8 = cv2.convertScaleAbs(imgsm, alpha=alp, beta=bet)

	print "Original file dimensions: ",amplitudeArr.shape[1],",",amplitudeArr.shape[0]
	print "New shape is :",imguint8.shape[1],",",imguint8.shape[0]
	print "X range is ",xstart,"--",xend
	print "Y range is ",ystart,"--",yend
	print "Original file statistics (ROI)"
	print "min   :",np.min(amparr)
	print "max   :",np.max(amparr)
	print "mean  :",np.mean(amparr)
	print "stdev :",np.std(amparr)
	print "Display statistics "
	print "min   :",np.min(imguint8)
	print "max   :",np.max(imguint8)
	print "mean  :",np.mean(imguint8)
	print "stdev :",np.std(imguint8)

	print "hit \'q\' key in the image window to close it..."

	# Use openCV to display the image
	#
	cv2.namedWindow(fileName)
	cv2.setMouseCallback(fileName,draw_circle)
	while(True):
		cv2.imshow(fileName,imguint8)
		k = cv2.waitKey(20)
		if k == ord('q'):
			break
		elif k == ord('a'):
			print "mouse coords are :",(mouseX*xshrink)+xstart,(mouseY*yshrink)+ystart
			print amparr[(mouseY*yshrink)+ystart][(mouseX*xshrink)+xstart]
		# cv2.waitKey(0)

	cv2.destroyAllWindows()

	# # extract region from within image
	# #

	# # convert the data into an array of numbers that we can do something with
	# #
	# amplitudeArr_full  = BandReadAsArray(amplitudeBand)
	# phaseArr      = BandReadAsArray(phaseBand)

	# amplitudeArr = amplitudeArr_full[3250:4250, 3000:6000]
	# # amplitudeArr = amplitudeArr_full

	# # Average down the data so that it will fit on the screen
	# #
	# sh     = amplitudeArr.shape
	# factor = 1
	# imgsm  = aveShape(amplitudeArr,(sh[0]//factor,sh[1]//factor))

	# # imgsm is now of type float64. We need to convert it to 8bit.
	# #
	# # imguint8 = cv2.convertScaleAbs(imgsm)
	# imguint8 = cv2.convertScaleAbs(imgsm, alpha=20, beta=0)

	# print "Original file dimensions: ",sh[0],",",sh[1]
	# print "New shape is :",imguint8.shape[0],",",imguint8.shape[1]
	# print "hit a key in the image window to close it..."

	# # Use openCV to display the image
	# #
	# cv2.imshow(fileName,imguint8)
	# cv2.waitKey(0)
	# cv2.destroyAllWindows()