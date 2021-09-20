import vaex
import numpy as np
import sys
from functools import reduce
import pandas as pd
import os
#getopt + yaml removed


path = sys.argv[2]

inputFiles = os.listdir(path)
inputFiles = [path + '/' + f for f in inputFiles]

print(inputFiles)


polarity = sys.argv[1]

dvI = vaex.open_many(inputFiles)
df_pandasI = dvI.to_pandas_df(["newMass","Intensity","polarity"], parallel=False)
if str(polarity) == "positive":
	df_pandasI = df_pandasI[df_pandasI["polarity"] == "Positive"]
if str(polarity) == "negative":
	df_pandasI = df_pandasI[df_pandasI["polarity"] == "Negative"]

df_pandasI.newMass = np.round(df_pandasI.newMass, decimals = 5)
test = df_pandasI
testSum = test.groupby(by=["newMass"]).count()

cutoff  = 10 

outputPathPositive = sys.argv[3]
outputPathNegative = sys.argv[4]

print(outputPathPositive)

if str(polarity) != "negative":
	for x in range(20):
		sd =  1
		cutoff = cutoff + 10 * sd
		print(cutoff)
		print("is cutoff")
		allSumsPassThreshold = testSum[testSum['Intensity'] > (cutoff) ]
		df_pandasIMass = allSumsPassThreshold.index.values
		isoMass = 1.00336
		df_pandasII = df_pandasIMass - isoMass
		isoIntersectII = reduce(np.intersect1d, (df_pandasIMass,df_pandasII))
		Outputname = 'DOE_autocredential_isopairs_polarity_positive_' + str(cutoff) + '.txt' 
		with open(outputPathPositive+ '/' + Outputname, 'w') as f:
			for item in isoIntersectII:
				f.write("%s\n" % item)



if str(polarity) != "positive":
	for x in range(20):
		sd =  1
		cutoff = cutoff + 10 * sd
		print(cutoff)
		print("is cutoff")
		allSumsPassThreshold = testSum[testSum['Intensity'] > (cutoff) ]
		df_pandasIMass = allSumsPassThreshold.index.values
		isoMass = 1.00336
		df_pandasII = df_pandasIMass - isoMass
		isoIntersectII = reduce(np.intersect1d, (df_pandasIMass,df_pandasII))
		Outputname = 'DOE_autocredential_isopairs_polarity_negative_' + str(cutoff) + '.txt'
		with open(outputPathNegative+ '/' + Outputname, 'w') as f:
			for item in isoIntersectII:
				f.write("%s\n" % item)

