import vaex
import sys
import pandas as pd
import os
import time

start = time.perf_counter()

theCommands = sys.argv[1]
testFile = theCommands.split(',')
print(testFile)


#output path example
#/home/lconnelly/Metabolomics/outputs/hdf5Files/Analysis1
outputPath = str(sys.argv[2])





#loop through the array of files, convert one at a time
for i in testFile:
	file = i
	#example file
	#file = /home/lconnelly/Metabolomics/rawDataFiles/Analysis1/file1.txt
	basename = os.path.basename(file)
	file2 = basename + ".hdf5"
	file2 = outputPath + '/' + file2
	print(file)
	df = pd.read_csv(file, sep=' ')
	print(" we read the file")
	df.columns=['File','Mass','Intensity','scan','polarity','Mstype']
	dfVaex = vaex.from_pandas(df)
	dfVaex.export_hdf5(file2, progress=True)


stop = time.perf_counter()

timeOut = outputPath.replace('hdf5Files', 'runTimes')
print(timeOut, 'is where runtime will go')


with open(timeOut + '/hdf5Files/perJob/jobRuntime' + basename, 'w') as f:
	f.write(str(stop-start) + ' \n')
	f.close()

