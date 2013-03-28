#!/usr/bin/python
# batchcorrNP
# Reproduces the functionality of batchcorr.exe for dark-correcting files produced
# using the a-Si detector at 1-ID at APS, using Python.
# Current version will run the correction on ALL *.ge2 files in the current working directory.
# Written by Chris Cochrane (cochranec@gmail.com), Dec. 2012

# Functionality presently requires the NumPy library - doing the file loading and 
# corrections would take up to 100 times longer using the native Python functions.
# NumPy can be downloaded from: http://www.numpy.org/
# Written using NumPy version 1.5.1, on Python 2.7.1+
# Not tested on Python 3.x; send any feedback to cochranec@gmail.com

import numpy
import os
import re
import glob
import sys
import argparse

# Command-line parser arguments - make everything more user friendly
parser = argparse.ArgumentParser(
  description='Dark correction and summing of GE2 files.',
  epilog='Written by Chris Cochrane, Dec. 2012. E-mail: cochranec@gmail.com')
parser.add_argument('--lo', type=int, nargs=1, default=[-1], help='Lower bound of run numbers.')
parser.add_argument('--hi', type=int, nargs=1, default=[-1], help='Upper bound of run numbers.')
parser.add_argument('--all','-a', action='store_true', default=True, help='Flag to perform dark correction on all GE2 files in present directory.')
parser.add_argument('--ndel', action='store_true', default=False, help='Save dark corrected output file for each frame in each GE2 file, as well as the sum files.  (WARNING: May use a LOT of space.)')
parser.add_argument('--drk', type=str, nargs=1, default=['dark'], help='Dark stub.  Some string that is unique to dark files.  Need not be the ENTIRE stub.  Default = "dark"')
parser.add_argument('--inpath', type=str, nargs=1, default=['.'], help='Path where the image files live. Default = ".".')
parser.add_argument('--outpath', type=str, nargs=1, default=['.'], help='Path where the corrected image files are saved. Default = ".".')
parser.add_argument('--genum', type=int, nargs=1, default=[2], help='GE number. Default = 2. NOT YET IMPLEMENTED')
clargs = parser.parse_args()

num_X = 2048
num_Y = 2048

lo = int(clargs.lo[0])
hi = int(clargs.hi[0])

genum = int(clargs.genum[0])

drk = clargs.drk[0]
inpath = clargs.inpath[0]
outpath = clargs.outpath[0]

if (lo != -1) or (hi != -1):
  clargs.all = False 

allfiles = glob.glob(inpath + os.sep + '*[0-9].ge' + str(genum))

#Recast the file reading function, which reduces runtime
fread = numpy.fromfile

#Pre-allocate memory space for arrays
sumvalues = numpy.zeros(num_X*num_Y,numpy.float32)
binvalues = numpy.zeros(num_X*num_Y,numpy.float32)
corrected = numpy.array(num_X*num_Y,numpy.float32)
darkvalues= numpy.array(num_X*num_Y,numpy.float32)
badPixels = numpy.array(num_X*num_Y,numpy.float32)

#Read in bad pixel data
# USER MUST DEFINE PATH TO BAD PIXEL INFORMATION
# Typically 'C:\DetectorData\1339.6\Full\1339.6Full_BadPixel_d.txt.img'
# Hard-coding this location in is inadvisable, as it makes the code less portable.
# Either a relative or an absolute path can be used.  It's probably safer to use an absolute path.
# (NB: the r prior to the string indicates a raw string, and must be included)
badPixFile = r'C:\DetectorData\1339.6\Full\1339.6Full_BadPixel_d.txt.img'
try:
  with open(badPixFile, mode='rb') as badPxobj:
    badPxobj.seek(8192)
    badPixels = fread(badPxobj, numpy.uint16, num_X * num_Y)
except IOError as e:
   print '\nUnable to access bad pixel information at ' + badPixFile
   print 'Ensure that the file exists, or change the "badPixFile" variable on line 26 to direct to the file location.\n'
   sys.exit()

# Find dark files.
# The 'dark' stub can be anywhere in the filename (not necessarily at the start)
darks = [x for x in allfiles if drk.lower() in x.lower()]

print len(darks), "candidate(s) for dark file found."
if len(darks) == 0:
  print "Double-check that dark file is properly located."
  sys.exit()
for i in range(len(darks)):
  print "   (" + str(i+1) + ") " + darks[i]
print "Which dark file should we use?"
dkI = 0
while dkI < 1 or dkI > len(darks):
  try:
    dkI = int(raw_input().strip())
  except:
    dkI = 0
    print('Choose one of the available options. [1-' + str(len(darks)) + ']')
darkfile = darks[dkI-1]
print "Using", darkfile

# Read in dark file
# Average over all the exposures in the file 
# This reduces the number of 'over reduced' pixels.
with open(darkfile, mode='rb') as darkobj:
  statinfo = os.stat(darkfile)
  nFrames = (statinfo.st_size - 8192) / (2 * num_X * num_Y)
  darkobj.seek(8192)
  for i in range(nFrames):
    binvalues = fread(darkobj, numpy.uint16, num_X * num_Y).astype('float32')
    sumvalues = sumvalues + binvalues

darkvalues = sumvalues / nFrames
sumvalues[:] = 0

# Pixel data is stored as 0, 1, 2, 3
# Any pixel with a non-zero value is deemed 'bad'
badInd = numpy.array(numpy.where(badPixels == 2))
badInd1= numpy.array(numpy.where(badPixels%2== 1))

print "Dark file and bad pixel data read successfully."

# Produce list of files to be dark corrected
if not clargs.all:
  print "High val is",hi,"; low val is",lo
  files = [x for x in allfiles if lo <= int(re.findall('_([0-9]*).ge2',x)[0]) <= hi and drk.lower() not in x.lower()]
else:
  print "Correcting all files"
  files = [x for x in allfiles if drk.lower() not in x.lower()]

print len(files), "of", len(allfiles), "GE files in directory are in range. ", len([x for x in allfiles if drk.lower() in x.lower()]), "dark files ignored."

# Maybe use a try: construct here
c = raw_input('Perform dark correction on all available files? ([y]/n)').strip()
if c.lower() == 'n':
  print "No dark correction will be performed.  Terminating script."
  sys.exit()
else:
  print "Proceeding with dark correction."

#Perform a loop over all files
for f in files:
  statinfo = os.stat(f)
  nFrames = (statinfo.st_size - 8192) / (2 * num_X * num_Y)
  print "\nReading:",f, "\nFile contains", nFrames," frames.  Summing and dark correcting."

  # Sum all values in this file
  (d, fout) = os.path.split(f)
  with open(f, mode='rb') as fileobj:
    fileobj.seek(8192)
    for i in range(nFrames):
      binvalues = fread(fileobj, numpy.uint16,num_X * num_Y).astype('float32')

      ### SUBTRACT BACKGROUND
      binvalues = binvalues - darkvalues

      ### OLD CORRECTION METHOD
      ### MIGHT WANT TO DO DIAGONAL NEIGHBORS
      binvalues[badInd] = (binvalues[badInd + 1] + binvalues[badInd - 1] + binvalues[badInd + num_X] + binvalues[badInd - num_X]) / 4

      # Set border region and negative pixels to 0
      binvalues[badInd1] = 0
      binvalues[numpy.where(binvalues<0)] = 0

      ### SUM OF ALL FRAMES
      sumvalues = sumvalues + binvalues

      ### DUMP IF WANT INDIVIDUAL FRAMES
      if clargs.ndel:
        frameName = outpath + os.sep + fout + '.frame.' + str(i+1) + '.cor'
        print 'saving frame number ' + str(i) + ' to ' + frameName
        with open(frameName, mode='wb') as outFile:
          binvalues.tofile(outFile)
          
  ### DUMP SUM OF FRAMES
  sumName = outpath + os.sep + fout + '.sum'
  print "saving sum to " + sumName
  with open(sumName, mode='wb') as outFile:
    sumvalues.tofile(outFile)
    
  ### DUMP IF WE WANT AVE OF FRAMES
  sumvalues = sumvalues/nFrames
  aveName = outpath + os.sep + fout + '.ave'
  print "saving ave to " + aveName
  with open(aveName, mode='wb') as outFile:
    sumvalues.tofile(outFile)
  
  sumvalues[:] = 0

print("Done")
