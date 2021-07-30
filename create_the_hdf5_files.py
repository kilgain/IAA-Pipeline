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




with open("/home/lconnelly/Metabolomics/Config.yaml") as file:
	config = yaml.safe_load(file)


#loop through the array of files
for i in testFile:
	file = config['samplesDirectory'] + "/" + i
	file2 = file + ".hdf5"
	file2 = file2.replace("rawDataFiles", "hdf5Files")
	print(file)
	df = pd.read_csv(file, sep=' ')
	print(" we read the file")
	df.columns=['File','Mass','Intensity','scan','polarity','Mstype']
	dfVaex = vaex.from_pandas(df)
	dfVaex.export_hdf5(file2, progress=True)
