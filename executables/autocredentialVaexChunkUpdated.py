import vaex
import numpy as np
import sys
from functools import reduce
import pandas as pd
import os

#Data import
path = sys.argv[2]

inputFiles = os.listdir(path)
inputFiles = [path + '/' + f for f in inputFiles]


polarity = sys.argv[1]


iterator = int(sys.argv[5])

start = iterator - 101

stop = iterator + 1


search = 'chunk' + str(iterator) + 'iso'

filter = [search in f for f in inputFiles]

inputFiles = [i for (i,v) in zip(inputFiles, filter) if v]

path = inputFiles[0]

inputFiles = os.listdir(inputFiles[0])

inputFiles = [path + '/' + f for f in inputFiles]

if len(inputFiles) == 0:
	print('no files detected for chunk', str(iterator))
	exit()

dvI = vaex.open_many(inputFiles)


#removes the polarity that is not of interest
if str(polarity) == "positive":
	dvI = dvI[dvI.polarity == "Positive"]
if str(polarity) == "negative":
	dvI = dvI[dvI.polarity == "Negative"]

dvI = dvI[dvI.newMass >= start]
dvI = dvI[dvI.newMass < stop]


#need to round masses to fifth decimal places to find isopairs
toRound = dvI['newMass']
toRound = np.round(toRound, decimals = 5)
dvI['Rounded'] = toRound


#Need to pull rounded masses into memory for set intersection
massesToTest = dvI.to_pandas_df(['Rounded'])
#need to add duplicate of column in order for aggregate to work properly
massesToTest['Rounded2'] = massesToTest['Rounded']

testSum = massesToTest.groupby(by=["Rounded2"]).count()

cutoff  = 0

outputPathPositive = sys.argv[3]
outputPathNegative = sys.argv[4]

#print(outputPathPositive)

if str(polarity) != "negative":
	for x in range(20):
		sd =  1
		cutoff = cutoff + 10 * sd
		#print(cutoff)
		#print("is cutoff")
		allSumsPassThreshold = testSum[testSum['Rounded'] > (cutoff) ]
		arrayOfMassesI = allSumsPassThreshold.index.values
		arrayOfMassesII = arrayOfMassesI - 1.00336
		isoIntersect = reduce(np.intersect1d, (arrayOfMassesI,arrayOfMassesII))
		Outputname = 'DOE_autocredential_isopairs_polarity_positive_' + str(cutoff) + '_chunk' + str(iterator) + '.txt' 
		with open(outputPathPositive+ '/' + Outputname, 'w') as f:
			for item in isoIntersect:
				f.write("%s\n" % item)



if str(polarity) != "positive":
	for x in range(20):
		sd =  1
		cutoff = cutoff + 10 * sd
		#print(cutoff)
		#print("is cutoff")
		allSumsPassThreshold = testSum[testSum['Rounded'] > (cutoff) ]
		arrayOfMassesI = allSumsPassThreshold.index.values
		arrayOfMassesII = arrayOfMassesI - 1.00336
		isoIntersect = reduce(np.intersect1d, (arrayOfMassesI,arrayOfMassesII))
		Outputname = 'DOE_autocredential_isopairs_polarity_negative_' + str(cutoff) +  '_chunk' + str(iterator) + '.txt'
		with open(outputPathNegative+ '/' + Outputname, 'w') as f:
			for item in isoIntersect:
				f.write("%s\n" % item)

