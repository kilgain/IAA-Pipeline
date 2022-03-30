import vaex
import numpy as np
import sys, getopt
from functools import reduce
import yaml
import os.path

theCommands = sys.argv[1]
testFile = theCommands.split(',')
testFile = testFile[0]

referenceFile = sys.argv[2]
 
outputPath = sys.argv[3]

reference = vaex.open(referenceFile)

 
reference = reference.to_pandas_df(parallel=False)

#Eliminate noise by looking only at isopairs in the top 5% of signals globally    
p_95 = reference.Intensity.quantile(0.95)

referenceQuant = reference[reference.Intensity.gt(p_95)]

#round so you can more easily find isopairs
lockTest1 = referenceQuant['Mass']
lockTest1Round = round(lockTest1, 4)
  

vaexOpen = testFile

queryVaex = vaex.open(vaexOpen)
query = queryVaex.to_pandas_df(parallel=False)


isopairs = []
p_95 = query.Intensity.quantile(0.95)
#print(zFilesSubset)
queryQuant = query[query.Intensity.gt(p_95)]
#loop through different ppm shifts between reference and test file. Find the one that gives the most isopairs. Adjust masses accordingly.
allowedRange = [-25,-24,-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
for y in allowedRange:
 lockTest = queryQuant['Mass']
 lockTest = lockTest - (1.00336 + (lockTest * .000001 *y))
 lockTestRound = round(lockTest, 4)
 lockTestCommon = reduce(np.intersect1d, (lockTest1Round,lockTestRound))
 #print(y)
 #print(len(lockTestCommon))
 isopairs.append(len(lockTestCommon))


#get the max offset
offset = allowedRange[isopairs.index(max(isopairs))]

queryVaex['newMass'] = queryVaex['Mass'] - offset*.000001*queryVaex['Mass']


#export the vaex file with the added column now
#vaexOpen = /home/lconnelly/Metabolomics/outputs/hdf5Files/Analysis1/file1.txt.hdf5
#vaexOpen = vaexOpen.replace("hdf5Files", "hdf5FilesNewColumn")

outputPath = outputPath + "/" +  os.path.basename(vaexOpen) + "newColumn_plate.hdf5"

#vaexOpenWrite=vaexOpen + "newColumn_plate.hdf5"
queryVaex.export_hdf5(outputPath)
