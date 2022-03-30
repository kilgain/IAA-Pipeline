import vaex
import pandas as pd
import sys
import os
import numpy
#import pathlib
#import functools
#remove yaml, getopt, functools, pathlib, numpy



fileRead = sys.argv[1]
fileRead = fileRead.split(',')


testMasses = pd.read_csv(sys.argv[2])

polarity = sys.argv[4]




testMasses = testMasses['x'].values


#fileRead = str(fileRead)
dv = vaex.open_many(fileRead)

#loop through the array of scans
#dv = dv[dv.polarity == str(polarity)]
for i in testMasses:
	#print(i)
	i = float(i)
	i = round(i, ndigits = 5)
	mass = i
	ppm = 10
	myMassRange = .000001 * ppm * mass
	subsetChromatogramI = dv[dv.newMass > (mass - myMassRange)]
	subsetChromatogramII = subsetChromatogramI[subsetChromatogramI.newMass < (mass + myMassRange)]

	tempPath = sys.argv[3]

	if not os.path.exists(tempPath+'/'+str(i)):
		os.makedirs(tempPath+'/'+str(i), exist_ok=True)

  
	base=os.path.basename(fileRead[0])
	fullName = tempPath + "/" + str(i) + '/' + str(i) + base+ ".txt"
	#print(fullName)
	#print("printing out now")
	#print(subsetChromatogramII)
	#print("is what we are printing")
	#print(len(subsetChromatogramII))
	if len(subsetChromatogramII) > 0:
		subsetChromatogramII.export_csv(path=fullName, parallel=False)
