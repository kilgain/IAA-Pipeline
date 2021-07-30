import vaex
import numpy as np
import sys, getopt
from functools import reduce
import yaml

theCommands = sys.argv[1]
testFile = theCommands.split(',')
testFile = testFile[0]
    
 
#dvI = vaex.open('/home/ahubbard/all_text_files_six_plates/setaria_files_all/setaria/to_combined_all_files/combined_chromatograms_with_header.csv.hdf5')
#dvI = vaex.open('/shares/ibaxter_share/private/ahubbard/all_text_files_six_plates/setaria_files_all/sorghum_vaex/combined_chromatograms_with_header.csv.hdf5')
#dv=dvI 
    

#align everything to the first file    
#this was from it was justo ne plate#zFilesSubsetI = vaex.open('/home/ahubbard/all_text_files_six_plates/setaria_files_all/setaria/Plate1HILIC_setaria')
#this was the one plate of the broadzFilesSubsetI = vaex.open('/home/ahubbard/full_broad_dataset/all_hdf5s/0000g_rx_XAV_iHMP2_HIL_PREFB01.rawoutput.txt.hdf5')
#zFilesSubsetI = vaex.open('/home/ahubbard/HILIC_plate_2/Setaria/HILIC_Set_1704_2_4_H14_85_PH_068.rawoutput.txt.hdf5')
#zFilesSubsetI = vaex.open('/home/ahubbard/text_files_for_all_DOE_metabolomics_samples/HILIC_plate_2/Setaria/HILIC_Set_1704_2_4_H14_85_PH_068.rawoutput.txt.hdf5')


#used for the HILIC
with open("/home/lconnelly/Metabolomics/Config.yaml") as file:
  config = yaml.safe_load(file)

zFilesSubsetI = vaex.open(config['referenceFile'])
testFile = config['samplesDirectory'] + "/" + testFile



#/home/ahubbard/text_files_for_all_DOE_metabolomics_samples/RPLC_plate1/Setaria/RPLC_Set_1704_2_4_H14_85_PH_068.rawoutput.txt.hdf5 

#/home/ahubbard/text_files_for_all_DOE_metabolomics_samples/HILIC_plate_2/Setaria/HILIC_Set_1704_2_4_H14_85_PH_068.rawoutput.txt.hdf5newColumn_plate.hdf5
#/home/ahubbard/text_files_for_all_DOE_metabolomics_samples/HILIC_plate_2/Setaria

 
zFilesSubsetI = zFilesSubsetI.to_pandas_df(parallel=False)



    
p_75 = zFilesSubsetI.Intensity.quantile(0.95)
#p_60 = lockTest2.quantile(0.6)
zFilesSubsetQuantI = zFilesSubsetI[zFilesSubsetI.Intensity.gt(p_75)]
    
lockTest1 = zFilesSubsetQuantI['Mass']
lockTest1Round = round(lockTest1, 4)
  
print("pulled out the first file")

#z = subsetChromatogramII.to_pandas_df(parallel=False)

#vaexOpen = "/home/ahubbard/" + testFile
vaexOpen = testFile
vaexOpen = vaexOpen.replace("rawDataFiles", "hdf5Files")
vaexOpen = vaexOpen + ".hdf5"
print(vaexOpen)
#vaexOpen = testFile

zFilesSubsetVaex = vaex.open(vaexOpen)
zFilesSubset = zFilesSubsetVaex.to_pandas_df(parallel=False)
#zFilesSubset = z[z.File.isin(files)]
#print("on 7")



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

#filesSubset = dvI[dvI.File.isin(files1)]


#get the max offset
offset = allowedRange[isopairs.index(max(isopairs))]
print(offset)
print(" is the offset")

#update the 

zFilesSubsetVaex['newMass'] = zFilesSubsetVaex['Mass'] - offset*.000001*zFilesSubsetVaex['Mass']


#export the vaex file with the added column now
vaexOpenWrite=vaexOpen + "newColumn_plate.hdf5"
zFilesSubsetVaex.export_hdf5(vaexOpenWrite)
