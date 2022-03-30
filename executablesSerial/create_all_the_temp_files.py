import vaex
import pandas as pd
import sys
import os
import numpy
import pathlib
import functools



fileRead = sys.argv[1]

testMasses = pd.read_csv(sys.argv[2])

#polarity = sys.argv[4]




testMasses = testMasses['x'].values


fileRead = str(fileRead)
dv = vaex.open(fileRead)

#loop through the array of masses and subset out all data within 20 ppm. Export the file.
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

  
	base=os.path.basename(fileRead)
	fullName = tempPath + "/" + str(i) + '/' + str(i) + base+ ".txt"
	#print(fullName)
	#print("printing out now")
	#print(subsetChromatogramII)
	#print("is what we are printing")
	#print(len(subsetChromatogramII))
	if len(subsetChromatogramII) > 0:
		subsetChromatogramII.export_csv(path=fullName, parallel=False)
