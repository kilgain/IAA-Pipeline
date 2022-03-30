import vaex
import pandas as pd
import sys
import os
import numpy
import pathlib
import functools

print('hello')

fPath = sys.argv[1]
print(fPath)
filepath = open(fPath, 'r')
fileRead = filepath.readline()
filepath.close()
fileRead = fileRead.split(',')



mPath = sys.argv[2]
print(mPath)
massesPath = open(mPath, 'r')
testMasses = massesPath.readline()
massesPath.close()
testMasses = testMasses.split(',')

#os.remove(mPath)

#polarity = sys.argv[4]




#testMasses = testMasses['x'].values



print(fileRead[0:3])
print(testMasses[0:3])

dv = vaex.open_many(fileRead)

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

	#if not os.path.exists(tempPath+'/'+str(i)):
		#os.makedirs(tempPath+'/'+str(i), exist_ok=True)

  
	#base=os.path.basename(fileRead)
	fullName = tempPath + "/" + str(i) + "merged.txt"
	#print(fullName)
	#print("printing out now")
	#print(subsetChromatogramII)
	#print("is what we are printing")
	#print(len(subsetChromatogramII))
	if len(subsetChromatogramII) > 0:
		subsetChromatogramII.export_csv(path=fullName, parallel=False, chunk_size=12500)
