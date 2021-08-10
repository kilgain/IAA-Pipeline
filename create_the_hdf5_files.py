import vaex
import numpy as np
import sys, getopt
from functools import reduce
import pandas as pd
import os
import yaml

theCommands = sys.argv[1]
testFile = theCommands.split(',')
print(testFile)


#/home/lconnelly/Metabolomics/outputs/hdf5Files/Analysis1
outputPath = sys.argv[2]


#runNumber = sys.argv[3]

#with open("/home/lconnelly/Metabolomics/Config.yaml") as file:
#	config = yaml.safe_load(file)


#loop through the array of files
for i in testFile:
	file = i
	#file = /home/lconnelly/Metabolomics/rawDataFiles/Analysis1/file1.txt

	basename = os.path.basename(file)


	file2 = basename + ".hdf5"
#	file2 = file2.replace("rawDataFiles", "hdf5Files")

	outputPath = outputPath + '/' + file2
	print(file)
	df = pd.read_csv(file, sep=' ')
	print(" we read the file")
	df.columns=['File','Mass','Intensity','scan','polarity','Mstype']
	dfVaex = vaex.from_pandas(df)
	dfVaex.export_hdf5(outputPath, progress=True)
