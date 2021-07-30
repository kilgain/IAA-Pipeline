import vaex
import numpy as np
import sys, getopt
from functools import reduce
import pandas as pd
import sys
import yaml

with open("/home/lconnelly/Metabolomics/Config.yaml") as file:
  config = yaml.safe_load(file)

config = config['samplesDirectory']
config = config.replace("rawDataFiles", "hdf5Files")
 
#getopt.getopt(args, options, [long_options])
fn1 = config + "/*.hdf5newColumn_plate.hdf5"
print(fn1)
#fn2 = '/home/ahubbard/full_broad_dataset/text_files_problematic_files_copies/additional_files_need_correcting/*.hdf5newColumn_plate.hdf5'
#fn3 = '/home/ahubbard/full_broad_dataset/all_hdf5s/new_columns/*hdf5newColumn_plate.hdf5'
#fn1 = '/home/ahubbard/text_files_for_all_DOE_metabolomics_samples/HILIC_plate_*/*/*hdf5newColumn_plate.hdf5'
#fn2 = '/home/ahubbard/HILIC_plate*/*/*hdf5newColumn_plate.hdf5'
#fn3 = '/home/ahubbard/full_broad_dataset/all_hdf5s/new_columns/*hdf5newColumn_plate.hdf5'



dvI = vaex.open_many([fn1])
df_pandasI = dvI.to_pandas_df(["newMass","Intensity"], parallel=False)
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
  Outputname = 'DOE_autocredential_isopairs_small_bins_cutofft_july14th_newMassFilter_true_' + str(cutoff) + '.txt' 
  with open("/home/lconnelly/Metabolomics/step1Conversion/step2Isolock/step3Autocredential/autocredentialCutoffs/"+Outputname, 'w') as f:
    for item in isoIntersectII:
      f.write("%s\n" % item)
