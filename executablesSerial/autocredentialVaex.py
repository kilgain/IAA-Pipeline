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

#print(inputFiles)


polarity = sys.argv[1]

dvI = vaex.open_many(inputFiles)


#removes the polarity that is not of interest
if str(polarity) == "positive":
	dvI = dvI[dvI.polarity == "Positive"]
if str(polarity) == "negative":
	dvI = dvI[dvI.polarity == "Negative"]


#need to round masses to fifth decimal places to find isopairs
toRound = dvI['newMass']
toRound = np.round(toRound, decimals = 5)
dvI['Rounded'] = toRound


#Need to pull rounded masses into memory for set intersection
massesToTest = dvI.to_pandas_df(['Rounded'])
#need to add duplicate of column in order for aggregate to work properly
massesToTest['Rounded2'] = massesToTest['Rounded']

testSum = massesToTest.groupby(by=["Rounded2"]).count()

cutoff  = 10 

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
		Outputname = 'DOE_autocredential_isopairs_polarity_positive_' + str(cutoff) + '.txt' 
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
		Outputname = 'DOE_autocredential_isopairs_polarity_negative_' + str(cutoff) + '.txt'
		with open(outputPathNegative+ '/' + Outputname, 'w') as f:
			for item in isoIntersectII:
				f.write("%s\n" % item)

