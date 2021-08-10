import vaex
import numpy as np
import sys, getopt
from functools import reduce
import pandas as pd
import sys
import yaml
import os

#with open("/home/lconnelly/Metabolomics/Config.yaml") as file:
#  config = yaml.safe_load(file)

#config = config['samplesDirectory']
#config = config.replace("rawDataFiles", "hdf5Files")
 
#getopt.getopt(args, options, [long_options])
path = sys.argv[2]

inputFiles = os.listdir(path)
inputFiles = [path + '/' + f for f in inputFiles]

print(inputFiles)




polarity = sys.argv[1]

#fn2 = '/home/ahubbard/full_broad_dataset/text_files_problematic_files_copies/additional_files_need_correcting/*.hdf5newColumn_plate.hdf5'
#fn3 = '/home/ahubbard/full_broad_dataset/all_hdf5s/new_columns/*hdf5newColumn_plate.hdf5'
#fn1 = '/home/ahubbard/text_files_for_all_DOE_metabolomics_samples/HILIC_plate_*/*/*hdf5newColumn_plate.hdf5'
#fn2 = '/home/ahubbard/HILIC_plate*/*/*hdf5newColumn_plate.hdf5'
#fn3 = '/home/ahubbard/full_broad_dataset/all_hdf5s/new_columns/*hdf5newColumn_plate.hdf5'



dvI = vaex.open_many(inputFiles)
df_pandasI = dvI.to_pandas_df(["newMass","Intensity","polarity"], parallel=False)
if str(polarity) == "positive":
	df_pandasI = df_pandasI[df_pandasI["polarity"] == "Positive"]
if str(polarity) == "negative":
	df_pandasI = df_pandasI[df_pandasI["polarity"] == "Negative"]

df_pandasI.newMass = np.round(df_pandasI.newMass, decimals = 5)
test = df_pandasI
testSum = test.groupby(by=["newMass"]).count()





allSumsPassThreshold = testSum[testSum['Intensity'] > (50) ]
cutoff = 40
#take only those mass-intenity bins with greater than the given number of signals
allSumsPassThreshold = testSum[testSum['Intensity'] > (cutoff) ]
df_pandasIMass = allSumsPassThreshold.index.values
isoMass = 1.00336
df_pandasII = df_pandasIMass - isoMass
isoIntersectII = reduce(np.intersect1d, (df_pandasIMass,df_pandasII))
#with open('broad_autocredential_isopairs_small_bins_40_cutoff_july14th_newMassFilter_true_sorghum_set.txt', 'w') as f:
#  for item in isoIntersectII:
#    f.write("%s\n" % item)
###loop through multiple potential thresholds:
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

