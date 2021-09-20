import vaex
import numpy as np
import sys, getopt
from functools import reduce
import yaml

theCommands = sys.argv[1]
testFile = theCommands.split(',')
testFile = testFile[0]

referenceFile = sys.argv[2]
 

zFilesSubsetI = vaex.open(referenceFile)

 
zFilesSubsetI = zFilesSubsetI.to_pandas_df(parallel=False)

    
p_75 = zFilesSubsetI.Intensity.quantile(0.95)

zFilesSubsetQuantI = zFilesSubsetI[zFilesSubsetI.Intensity.gt(p_75)]
    
lockTest1 = zFilesSubsetQuantI['Mass']
lockTest1Round = round(lockTest1, 4)
  
print("pulled out the first file")

vaexOpen = testFile
print(vaexOpen)

zFilesSubsetVaex = vaex.open(vaexOpen)
zFilesSubset = zFilesSubsetVaex.to_pandas_df(parallel=False)


isopairs = []
p_75 = zFilesSubset.Intensity.quantile(0.95)
print(zFilesSubset)
zFilesSubsetQuant = zFilesSubset[zFilesSubset.Intensity.gt(p_75)]
allowedRange = [-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
for y in allowedRange:
 lockTest = zFilesSubsetQuant['Mass']
 lockTest = lockTest - (1.00336 + (lockTest * .000001 *y))
 lockTestRound = round(lockTest, 4)
 lockTestCommon = reduce(np.intersect1d, (lockTest1Round,lockTestRound))
 print(y)
 print(len(lockTestCommon))
 isopairs.append(len(lockTestCommon))


#get the max offset
offset = allowedRange[isopairs.index(max(isopairs))]
print(offset)
print(" is the offset")


zFilesSubsetVaex['newMass'] = zFilesSubsetVaex['Mass'] - offset*.000001*zFilesSubsetVaex['Mass']


#export the vaex file with the added column now
#vaexOpen = /home/lconnelly/Metabolomics/outputs/hdf5Files/Analysis1/file1.txt.hdf5
vaexOpen = vaexOpen.replace("hdf5Files", "hdf5FilesNewColumn")

vaexOpenWrite=vaexOpen + "newColumn_plate.hdf5"
zFilesSubsetVaex.export_hdf5(vaexOpenWrite)
