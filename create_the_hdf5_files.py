import vaex
import sys
import pandas as pd
import os
#numpy functool getopt removed

theCommands = sys.argv[1]
testFile = theCommands.split(',')
print(testFile)


#output path example
#/home/lconnelly/Metabolomics/outputs/hdf5Files/Analysis1
outputPath = sys.argv[2]





#loop through the array of files
for i in testFile:
	file = i
	#example file
	#file = /home/lconnelly/Metabolomics/rawDataFiles/Analysis1/file1.txt
	basename = os.path.basename(file)
	file2 = basename + ".hdf5"
	outputPath = outputPath + '/' + file2
	print(file)
	df = pd.read_csv(file, sep=' ')
	print(" we read the file")
	df.columns=['File','Mass','Intensity','scan','polarity','Mstype']
	dfVaex = vaex.from_pandas(df)
	dfVaex.export_hdf5(outputPath, progress=True)
