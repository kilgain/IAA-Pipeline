import vaex
import numpy as np
import pandas as pd
import sys, getopt
from functools import reduce
from pathlib import Path
import os
import yaml




fileRead = sys.argv[1]
#testMasses = pd.read_csv("/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/massesFileRounded.txt")

testMasses = pd.read_csv(sys.argv[2])


#testMasses = pd.read_csv("/home/ahubbard/all_text_files_six_plates/setaria_files_all/setaria/all_vaex_files_readin/just_the_standards.txt")
#testMasses = pd.read_csv("/home/ahubbard/all_text_files_six_plates/setaria_files_all/setaria/all_vaex_files_readin/just_the_standards.txt")
#testMasses  = pd.read_csv("/home/ahubbard/full_broad_dataset/all_autocredential_scripts_for_qcing/autocredential_85_per_bin_all_600_april_13_2021.txt")

##testMasses  = pd.read_csv("/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/setaria_plate_9_standards/standards.txt")


#testMasses  = pd.read_csv("/home/ahubbard/just_citulline_and_isotopologue.txt")

#testMasses  = pd.read_csv("/home/ahubbard/just_citulline_and_isotopologue.txt")

#testMasses = pd.read_csv("/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/centroided_masses_june16h_bins_filtered_at50_counts.txt")



#testMasses = pd.read_csv("/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/just_isotopes_citrulline.txt")


#testMasses = pd.read_csv("/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/centroided_masses_july16th_bins_filtered_at100_counts.txt")

#testMasses = pd.read_csv("/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/standards_from_louis_to_pullout.txt")

###testMasses = pd.read_csv("/home/ahubbard/yeastStandards.csv")



###testMasses = pd.read_csv("/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/unknown_lipids_from_michael_just_masses.txt")
#/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files



#/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/standards_from_louis_to_pullout.txt 


#just_isotopes_citrulline.txt

#/home/ahubbard/just_citulline_and_isotopologue.txt

#myMass = 176.1034
#myIsotope = 176.1034 + 1.00336

#testMasses = [myMass, myIsotope]

#/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/setaria_plate_9_standards/standards.txt

#testMasses  = pd.read_csv("/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/not_deisotope_but_AC_setaria_HILIC.txt")


#/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/not_deisotope_but_AC_setaria_HILIC.txt




testMasses = testMasses['x'].values 
#testMasses = testMasses[0:10]
#testMasses = [176.1031,177.1001,177.1073, 182.12326, 183.1203,183.1275]
#print(testMasses)
#print(" is testMasses")


#testMass = [176.1034]
#testMasses = round(testMasses, ndigits = 5)

#testMasses = np.unique(dvI[['scan']])
fullarray = []


#fileRead = "/home/ahubbard/" + str(fileRead[0])
fileRead = str(fileRead)


dv = vaex.open(fileRead)
#loop through the array of scans
for i in testMasses:
	print(i)
  #output the glutamine info into text file
  #mass = i
	i = float(i)
	i = round(i, ndigits = 5)
	mass = i
  #actually, make it an isotopologue!
  #mass = mass + 1.00336
	ppm = 10
	myMassRange = .000001 * ppm * mass
	subsetChromatogramI = dv[dv.newMass > (mass - myMassRange)]
	subsetChromatogramII = subsetChromatogramI[subsetChromatogramI.newMass < (mass + myMassRange)]

  #z = subsetChromatogramII.to_pandas_df(parallel=False)
  #tempPath = "/home/ahubbard/ST001293/converted_data_for_autocredential/positive_data/temp_files"
  #tempPath = "/home/ahubbard/ST001293/converted_data_for_autocredential/positive_data/temp_files_xcms"

  #tempPath = "/home/ahubbard/ST001293/converted_data_for_autocredential/positive_data/autocredential_strict_nov_2020"  
#  tempPath = "/home/ahubbard/ST001293/converted_data_for_autocredential/positive_data/autocredential_lite_nov_2020"
#  tempPath = "/home/ahubbard/full_broad_dataset/temp_files_random"
#  tempPath = "/home/ahubbard/full_broad_dataset/temp_files_just_one_isopairs_jan14_2020"
#  tempPath = "/home/ahubbard/full_broad_dataset/temp_files_7900_matching_pub"
#  tempPath = "/home/ahubbard/full_broad_dataset/ALL_QI_MASSES"
  #tempPath = "/home/ahubbard/ST001293/converted_data_for_autocredential/positive_data/AC_ultralite_nov_25_2020"


  #tempPath = "/home/ahubbard/all_text_files_six_plates/setaria_files_all/setaria/to_combined_all_files/temp_files_real_masses"


#  tempPath = "/home/ahubbard/all_text_files_six_plates/setaria_files_all/setaria/to_combined_all_files/temp_files_more_stringent_filtering"
#  tempPath = "/home/ahubbard/all_text_files_six_plates/setaria_files_all/setaria/to_combined_all_files/temp_files_single_isopair_but_filtered"
  #temp_files_single_isopair_but_filtered

  #tempPath = "/shares/ibaxter_share/private/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/setaria_plate9_chromatogram_files_final/" + str(i)
  #tempPath = "/shares/ibaxter_share/private/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/setaria_plate_9_standards/" + str(i)
  #tempPath = "/shares/ibaxter_share/private/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/all_broad_autocredential/" + str(i)

  #tempPath = "/shares/ibaxter_share/private/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/HILIC_12_plates_setaria/" + str(i)
  #tempPath = "/shares/ibaxter_share/private/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/HILIC_12_plates_setaria_standards/" + str(i)
  #tempPath = "/shares/ibaxter_share/private/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/temp_files_broad_citrulline/" + str(i)
  #tempPath = "/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/june16th_bins_filtered_at_50_counts/" + str(i)


  #tempPath = "/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/broad_isotopes_multiple/" + str(i)
  #tempPath = "/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/DOE_setaria_and_HILIC/" + str(i)
  #tempPath = "/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/standards_DOE_HILIC_louis_matches/" + str(i)
  #tempPath = "/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/standards_louis_more_comprehensive/" + str(i)
  #tempPath = "/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/michael_unknown_lipids/" + str(i)

  #/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/standards_DOE_HILIC_louis_matches

	with open("/home/lconnelly/Metabolomics/Config.yaml") as file:
		config = yaml.safe_load(file)


	root = config['metabolomicsRoot']
	tempPath = root + "/step1Conversion/step2Isolock/step3Autocredential/step4Subset/tempFiles5C/" + str(i)


	if not os.path.exists(tempPath):
		os.makedirs(tempPath, exist_ok=True)

  

	#Path(tempPath).mkdir(parents=True, exist_ok=True)


  #AC_ultralite_nov_25_2020
	base=os.path.basename(fileRead)
	fullName = tempPath + "/" + str(i) + base+ ".txt"
	print(fullName)
	print("printing out now")
	print(subsetChromatogramII)
	print("is what we are printing")
	print(len(subsetChromatogramII))
	if len(subsetChromatogramII) > 0:
		subsetChromatogramII.export_csv(path=fullname)
#old path = fullname
