from dask_jobqueue.htcondor import HTCondorCluster
from dask.distributed import Client, progress, wait
import os
import yaml
import glob
from contextlib import suppress
import pandas as pd
import sys
import time


def runScript(x,y,z="", xx = "", xy = ""):
	os.system("Rscript " + x + " " + y + " " + z + " " + xx + " " + xy)


def runPythonScript(x,y = "", z = "", xx = "", xy = "", xz = "", yx = ""):
	os.system("python3 " + x + " " + y + " " + z + " " + xx + " " + xy + " " + xz + " " + yx)

def runMono(x,y,z='', xx=''):
	os.system("mono " + x + " " + y + " " + z + " " + xx)


#root = directory which contains dask script
root = '/home/lconnelly/MetabolomicsNewDeployment'

#all outputs will be sorted based on this run number
runNumber = 'DOE'

#	 0 1 2 3 4 5 6 7 8 9 10
toRun = [0,0,0,0,0,0,0,0,0,0,1]

#Converts .raw files to text files (skip if raw files have already been processed)
runRawFileProcessing = False


#Determine noise cutoff automatically, or manually input a value (multiples of 10. 10 is enough for small sample sizes. Manual used if autoSelect is False.)

autoSelect = False
#value = '40' #Used in DOE for full dataset
value = '10' #used in prelim

#note, initialization should always be set to 1
#indexes for toRun
#0 - initialization (will skip if "file" column is detected)
#1 - conversion to hdf5
#2 - isolock newColumn
#3 - autocredential positive mode
#4 - autocredential negative mode
#5 - centroidFxn positive mode
#6 - centroidFxn negative mode
#7 - tempFile creation for positive
#8 - tempFile creation for negative
#9 - massFeature for positive
#10- massFeature for negative

if __name__ == "__main__":

	if not os.path.exists(root + "/outputs/runTimes/" + runNumber + '/rawFileProcessing/perJob'):
		os.makedirs(root + "/outputs/runTimes/" + runNumber + '/rawFileProcessing/perJob', exist_ok=True)


	with open("/home/lconnelly/MetabolomicsNewDeployment/Config.yaml", "r") as file:
		print(type(file))
		config = yaml.safe_load(file)
	if runRawFileProcessing:
		start = time.perf_counter()
		with open("/home/lconnelly/MetabolomicsNewDeployment/Config.yaml", "r") as file:
			config = yaml.safe_load(file)
		#directoryOfRawFiles = '/home/ahubbard/all_setaria_and_sorghum_raw_files_to_backup/Setaria_Sorghum_Metabolomics'
		#directoryOfRawFiles = '/home/ahubbard/Metabolomics/rawDataFiles/testShrikaar'		
		#directoryOfRawFiles = '/home/ahubbard/ChengZhao_Carnegie'

		directoryOfRawFiles = root + '/rawDataFiles/' + runNumber
		print(directoryOfRawFiles, 'is the raw file directory')

		#for HILIC
		#rawFiles = glob.glob(directoryOfRawFiles+'/*/*HILIC*/*.raw', recursive = True)

		#for RPLC
		#rawFiles = glob.glob(directoryOfRawFiles+'/*/*RPLC*/*.raw', recursive = True)
		
		#for default
		rawFiles = [os.path.join(directoryOfRawFiles, file) for file in os.listdir(directoryOfRawFiles)]		



                #TESTING:
		#rawFiles = rawFiles[0:8]
		#print('HERE IS LENGTH OF RAWFILES:', len(rawFiles))
		####
		print(rawFiles[0], 'is the first rawFile')
		print(rawFiles[1], 'is second rawFile')
		outputDirectory = root + '/outputs/processedRawFiles/' + runNumber
		print(outputDirectory, 'is where the processes .raw files will be located')
		
		if os.path.exists(outputDirectory):
			print('processed raw files detected, only running files not present in output directory')
			rawFilesBasename = []
			for path in rawFiles:
				rawFilesBasename.append(os.path.basename(path))
			pathData = pd.DataFrame(list(zip(rawFiles,rawFilesBasename)), columns=['fullPath', 'basename']).reset_index(drop=True)
			rawFilesBasename = set(rawFilesBasename)
			rawFilesAlreadyProcessed = os.listdir(outputDirectory)
			rawFilesAlreadyProcessed = [x[:len(x)-10] for x in rawFilesAlreadyProcessed]
			rawFilesAlreadyProcessed = set(rawFilesAlreadyProcessed)
			rawFilesToProcess = rawFilesBasename.difference(rawFilesAlreadyProcessed)
			rawFilesToProcess = list(rawFilesToProcess)
			rawFiles = []
			for path in range(len(rawFilesToProcess)):
				rawFiles.append(list(pathData['fullPath'][pathData['basename'] == rawFilesToProcess[path]])[0])


		if not os.path.exists(outputDirectory):
			os.makedirs(outputDirectory, exist_ok=True)
		print('LENGTH OF RAW FILES TO BE PROCESSED:', len(rawFiles))
		cluster = HTCondorCluster(n_workers = 1, cores=1, memory = "8 GB", disk = "4 GB",local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		cluster.scale(jobs = 300)
		client = Client(cluster)
		processed = []
		for f in rawFiles:
			processed.append(client.submit(runMono, root + "/executables/fileProcessingMono/fileProcessor.exe", f, outputDirectory))
		print('Processing Raw Files!')
		progress(processed)
		client.gather(processed)
		client.shutdown()

		stop = time.perf_counter()
		runTime = stop - start

		f = open(root + '/outputs/runTimes/' + runNumber + '/rawFileProcessingOverall.txt', 'w')
		f.writelines(str(runTime))
		f.close()

		#create fileList for anovalign later:
		fileList = os.listdir(root + '/outputs/processedRawFiles/' + runNumber)
		with open(root + '/fileListForAnovalign' + runNumber + '.txt', 'w') as f:
			for item in fileList:
				f.write("%s\n" % item)



	#Step1: Add File Column and Collapse Duplicate Masses
	if toRun[0]:
		#initiates job cluster (for Condor)
		cluster = HTCondorCluster(n_workers = 1, cores=1, memory = "8 GB", disk = "4 GB",local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		with open("/home/lconnelly/MetabolomicsNewDeployment/Config.yaml", "r") as file:
			config = yaml.safe_load(file)
		initConfig = root  + '/outputs/processedRawFiles/' +  runNumber
		print(initConfig, "is config path")
		files = os.listdir(initConfig)
		#sets job cluster size equal to number of files
		cluster.scale(jobs = 250)
		#initiates job client that assigns jobs using created cluster
		client = Client(cluster)
		#list of job futures
		processed = []
		
		filesToRun = []
		for f in files:
			toCheck = open(initConfig + '/' + f, 'r')
			line1 = toCheck.readline()
			toCheck.close()
			if not '"File"' in line1:
				filesToRun.append(initConfig + '/' + f)

		#if not '"File"' in line1:
		#submits all files for formatting
		if len(filesToRun) > 0:
			for f in filesToRun:
				processed.append(client.submit(runScript, root + "/executables/fileAdd.R", f))
		#progress bar
		print('Formatting Files!')
		progress(processed)
		#prevents progress until all jobs are done
		client.gather(processed)
		#destroys job cluster
		client.shutdown()


	#Step2: Convert .txt Files to HDF5 Format

	if toRun[1] == 1:
		
		if not os.path.exists(root + "/outputs/runTimes/" + runNumber + '/hdf5Files/perJob'):
			os.makedirs(root + "/outputs/runTimes/" + runNumber + '/hdf5Files/perJob', exist_ok=True)


		logPath = root + '/logs/hdf5Files/' + runNumber
		
		if os.path.exists(logPath):
			for x, y, z in os.walk(logPath, topdown=False):
				for name in z:
					#print('removing', os.path.join(x, name))
					os.remove(os.path.join(x, name))

		if not os.path.exists(logPath):
			os.makedirs(logPath, exist_ok=True)


		
		start = time.perf_counter()

		cluster = HTCondorCluster(cores=1, n_workers=1, memory = "4 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True}, log_directory=logPath)
		cluster.scale(jobs = 250)
		client = Client(cluster)
		processed = []


		initConfig = root  + '/outputs/processedRawFiles/' +  runNumber
		outputPath = root + "/outputs/hdf5Files/" + runNumber

		#FOR MANUAL SELECTION OF PROCESSED RAW FILES:
		initConfig = '/home/ahubbard/Metabolomics/rawDataFiles/T4_T6_T8_setariaII'


		if not os.path.exists(outputPath):
			os.makedirs(outputPath, exist_ok=True)

	
		filesBasename = os.listdir(initConfig)
		print(filesBasename[0])
		files = [initConfig + '/' + f for f in filesBasename]
		print(len(files), "is length of input for conversion")
	

		pathData = pd.DataFrame(list(zip(files,filesBasename)), columns=['fullPath', 'basename']).reset_index(drop=True)
		filesBasename = set(filesBasename)
		filesAlreadyProcessed = os.listdir(outputPath)
		filesAlreadyProcessed = [x[:len(x)-5] for x in filesAlreadyProcessed]
		print(filesAlreadyProcessed[0:3])
		filesAlreadyProcessed = set(filesAlreadyProcessed)
		filesToProcess = filesBasename.difference(filesAlreadyProcessed)
		filesToProcess = list(filesToProcess)
		files = []
		for path in range(len(filesToProcess)):
			files.append(list(pathData['fullPath'][pathData['basename'] == filesToProcess[path]])[0])

		print(len(files), 'is length of input for conversion after set operations')
		print("is output directory")
		print(outputPath)

		for f in files:
			processed.append(client.submit(runPythonScript, root + "/executables/create_the_hdf5_files.py", f, outputPath))
		print('Converting Files To Hdf5!')
		print(len(processed))
		progress(processed)
		client.gather(processed)
		client.shutdown()


		stop = time.perf_counter()

		runTime = stop - start

		f = open(root + '/outputs/runTimes/' + runNumber + '/hdf5Files/HDF5conversionOverall.txt', 'w')
		f.writelines(str(runTime))
		f.close()

		time.sleep(3)
		outputPath = root + '/outputs/runTimes/' + runNumber + '/hdf5Files'
		runPythonScript(root + '/executables/parse_logs.py', logPath, outputPath) 
		

	#Step3: Run Isolock On Each File
	if toRun[2] == 1:


		
		if not os.path.exists(root + "/outputs/runTimes/" + runNumber + '/hdf5FilesNewColumn/perJob'):
			os.makedirs(root + "/outputs/runTimes/" + runNumber + '/hdf5FilesNewColumn/perJob', exist_ok=True)

		logPath = root + '/logs/hdf5FilesNewColumn/' + runNumber

		if os.path.exists(logPath):
			for x, y, z in os.walk(logPath, topdown=False):
				for name in z:
					os.remove(os.path.join(x, name))

		if not os.path.exists(logPath):
			os.makedirs(logPath, exist_ok=True)



		start = time.perf_counter()
		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True}, log_directory = logPath)
		cluster.scale(jobs = 250)
		client = Client(cluster)
		processed = []

		files = os.listdir(root + '/outputs/hdf5Files/' + runNumber)
		files = [root + '/outputs/hdf5Files/' + runNumber + '/' + f for f in files]
		refFile = files[0]
		#refFile = '/home/lconnelly/MetabolomicsNewDeployment/outputs/hdf5Files/DOE/HILIC_Set_1342_3_3_H11_85_TB_setaria_12_0127.rawoutput.txt.hdf5'
		#refFile = '/home/lconnelly/MetabolomicsNewDeployment/outputs/hdf5Files/broadDataTest/0411_XAV_iHMP2_HIL-SM-7MCUF.rawoutput.txt.hdf5'

		outFileLocation = root + '/refFileForIsolock' + runNumber + '.txt'
		with open(outFileLocation, "w") as fileOut:
			fileOut.write(refFile)
		fileOut.close()
		
		print(refFile)
		print("is refFile")




		outputPath = root + "/outputs/hdf5FilesNewColumn/" + runNumber
		if not os.path.exists(outputPath):
                	os.makedirs(outputPath, exist_ok=True)

	

		for f in files:
			processed.append(client.submit(runPythonScript, root + "/executables/isolockChunk.py", f, refFile, outputPath))
		print('Running Isolock!')
		progress(processed)
		client.gather(processed)
		client.shutdown()
		stop = time.perf_counter()

		runTime = stop - start
		
		f = open(root + '/outputs/runTimes/' + runNumber + '/hdf5FilesNewColumn/isolockRuntimeOverall.txt', 'w')
		f.writelines(str(runTime))
		f.close()


		time.sleep(3)
		outputPath = root + '/outputs/runTimes/' + runNumber + '/hdf5FilesNewColumn'
		runPythonScript(root + '/executables/parse_logs.py', logPath, outputPath) 	


	#Step4: Run Autocredential
	#Positive Mode Autocredential

	if toRun[3] == 1:


		if not os.path.exists(root + "/outputs/runTimes/" + runNumber + '/autocredentialPositive/perJob'):
			os.makedirs(root + "/outputs/runTimes/" + runNumber + '/autocredentialPositive/perJob', exist_ok=True)

		logPath = root + '/logs/autocredentialPositive/' + runNumber


		#clear old logs before adding new
		if os.path.exists(logPath):
			for x, y, z in os.walk(logPath, topdown=False):
				for name in z:
					os.remove(os.path.join(x, name))

				
		if not os.path.exists(logPath):
			os.makedirs(logPath, exist_ok=True)


		start = time.perf_counter()
		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "10 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv":True}, log_directory = logPath)
		cluster.scale(jobs = 250)
		client = Client(cluster)
		processed = []


		path = root + '/outputs/hdf5FilesNewColumn/' + runNumber
		#inputFiles = os.listdir(path)
		#inputFiles = [path + '/' + f for f in inputFiles]

		#print(inputFiles[0])
		#print('is first input for autocred')


		outputPathPositive = root + "/outputs/autocredentialPositive/" + runNumber
		outputPathNegative = root + "/outputs/autocredentialNegative/" + runNumber

		if not os.path.exists(outputPathPositive):
			os.makedirs(outputPathPositive, exist_ok=True)


		if not os.path.exists(outputPathNegative):
			os.makedirs(outputPathNegative, exist_ok=True)

		for i in range(65, 1205, 5):
			processed.append(client.submit(runPythonScript, root + "/executables/autocredentialVaexChunkUpdated.py", "positive", path, outputPathPositive, outputPathNegative, str(i)))

		print('Running Autocredential Positive Mode!')
		progress(processed)
		client.gather(processed)
		client.shutdown()
		

		#stitch together chunks
		allOutputs = os.listdir(outputPathPositive)
		for i in range(10,201,10):
			output = root + '/outputs/autocredentialPositive/' + runNumber + '/DOE_autocredential_isopairs_polarity_positive_' + str(i) + '.txt'
			w = open(output, 'w')
			search = '_' + str(i) + '_chunk'
			toSelect = [search in f for f in allOutputs]


			toStitch = [i for (i,v) in zip(allOutputs, toSelect) if v]

			for entry in toStitch:
				f = open(root + '/outputs/autocredentialPositive/' + runNumber + '/' + entry, 'r')
				lines = f.readlines()
				f.close()
				if len(lines) == 0:
					continue
				w.writelines(lines)
			w.close()
			toSort = pd.read_table(output, header = None)
			toSort = toSort.sort_values(0, axis=0, ascending=True)
			toSort = pd.unique(toSort[0])
			toSort = pd.DataFrame(toSort)
			toSort.to_string(output, header = False, index = False, index_names = False)
			for entry in toStitch:
				os.remove(root + '/outputs/autocredentialPositive/' + runNumber + '/' + entry)

		stop = time.perf_counter()

		runTime = stop - start
		
		f = open(root + '/outputs/runTimes/' + runNumber + '/autocredentialPositive/autocredentialPositiveOverall.txt', 'w')
		f.writelines(str(runTime))
		f.close()
		
		time.sleep(3)
		outputPath = root + '/outputs/runTimes/' + runNumber + '/autocredentialPositive'
		runPythonScript(root + '/executables/parse_logs.py', logPath, outputPath)



	#Negative Mode Autocredential

	if toRun[4] == 1:


		if not os.path.exists(root + "/outputs/runTimes/" + runNumber + '/autocredentialNegative/perJob'):
			os.makedirs(root + "/outputs/runTimes/" + runNumber + '/autocredentialNegative/perJob', exist_ok=True)


		logPath = root + '/logs/autocredentialNegative/' + runNumber

		#clear old logs before adding new
		if os.path.exists(logPath):
			for x, y, z in os.walk(logPath, topdown=False):
				for name in z:
					os.remove(os.path.join(x, name))


		if not os.path.exists(logPath):
			os.makedirs(logPath, exist_ok=True)


		start = time.perf_counter()
		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "10 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv":True}, log_directory = logPath)
		cluster.scale(jobs = 250)
		client = Client(cluster)
		processed = []

		path = root + '/outputs/hdf5FilesNewColumn/' + runNumber
		#inputFiles = os.listdir(path)
		#inputFiles = [path + '/' + f for f in inputFiles]


		outputPathPositive = root + "/outputs/autocredentialPositive/" + runNumber
		outputPathNegative = root + "/outputs/autocredentialNegative/" + runNumber


		if not os.path.exists(outputPathPositive):
			os.makedirs(outputPathPositive, exist_ok=True)


		if not os.path.exists(outputPathNegative):
			os.makedirs(outputPathNegative, exist_ok=True)


		for i in range(65,1205,5):
			processed.append(client.submit(runPythonScript, root + "/executables/autocredentialVaexChunkUpdated.py", "negative", path, outputPathPositive, outputPathNegative, str(i)))
		print('Running Autocredential Negative Mode!')
		progress(processed)
		client.gather(processed)
		client.shutdown()


		#stitch together outputs
		allOutputs = os.listdir(outputPathNegative)
		for i in range(10,201,10):
			output = root + '/outputs/autocredentialNegative/' + runNumber + '/DOE_autocredential_isopairs_polarity_negative_' + str(i) + '.txt'
			w = open(output, 'w')
			search = '_' + str(i) + '_chunk'
			toSelect = [search in f for f in allOutputs]
			toStitch = [i for (i,v) in zip(allOutputs, toSelect) if v]
			for entry in toStitch:
				f = open(root + '/outputs/autocredentialNegative/' + runNumber + '/' + entry, 'r')
				lines = f.readlines()
				f.close()
				if len(lines) == 0:
					continue
				w.writelines(lines)
			w.close()
			toSort = pd.read_table(output, header = None)
			toSort = toSort.sort_values(0, axis=0, ascending=True)
			toSort = pd.unique(toSort[0])
			toSort = pd.DataFrame(toSort)
			toSort.to_string(output, header = False, index = False, index_names = False)
			for entry in toStitch:
				os.remove(root + '/outputs/autocredentialNegative/' + runNumber + '/' + entry)

		stop = time.perf_counter()
		
		runTime = stop - start

		f = open(root + '/outputs/runTimes/' + runNumber + '/autocredentialNegative/autocredentialNegativeOverall.txt', 'w')
		f.writelines(str(runTime))
		f.close()

		time.sleep(3)
		outputPath = root + '/outputs/runTimes/' + runNumber + '/autocredentialNegative'
		runPythonScript(root + '/executables/parse_logs.py', logPath, outputPath)



	#Step5: Centroid the Masses

	#Positive Mode Centroid

	if toRun[5] == 1:

		if not os.path.exists(root + "/outputs/runTimes/" + runNumber + '/centroidFxnPositive'):
			os.makedirs(root + "/outputs/runTimes/" + runNumber + '/centroidFxnPositive', exist_ok=True)

		start = time.perf_counter()
		cluster = HTCondorCluster(cores=1, n_workers=1, memory="4 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		cluster.scale(jobs = 10)
		client = Client(cluster)
		processed = []

		path = config['metabolomicsRoot'] + "/outputs/autocredentialPositive/" + runNumber
		files = os.listdir(path)
		files = [path + '/' + f for f in files]

		outputPath = config['metabolomicsRoot'] + '/outputs/centroidFxnPositive/' + runNumber
	

		print("here is outputPath for centroidFxn:")
		print(outputPath)

		print("here is first file input for centroid:")
		print(files[0])

	
		if not os.path.exists(outputPath):
			os.makedirs(outputPath, exist_ok=True)


		print("here is script path:")
		print(root + "/executables/centroidFxn.R")

		for f in files:
			processed.append(client.submit(runScript, root + "/executables/centroidFxn.R", f, outputPath))
		print('Centroiding Positive Mode Masses!')
		progress(processed)
		client.gather(processed)
		client.shutdown()

		stop = time.perf_counter()

		runTime = stop - start

		f = open(root + '/outputs/runTimes/' + runNumber + '/centroidFxnPositive/centroidPositiveOverall.txt', 'w')
		f.writelines(str(runTime))
		f.close()

	#Negative Mode Centroid:
	
	if toRun[6] == 1:
	
		if not os.path.exists(root + "/outputs/runTimes/" + runNumber + '/centroidFxnNegative'):
			os.makedirs(root + "/outputs/runTimes/" + runNumber + '/centroidFxnNegative', exist_ok=True)


		start = time.perf_counter()
		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		cluster.scale(jobs = 10)
		client = Client(cluster)
		processed = []


		path = config['metabolomicsRoot'] + "/outputs/autocredentialNegative/" + runNumber
		files = os.listdir(path)
		files = [path + '/' + f for f in files]
		outputPath = config['metabolomicsRoot'] + '/outputs/centroidFxnNegative/' + runNumber

		print(outputPath)
		print(files[0])

		if not os.path.exists(outputPath):
			os.makedirs(outputPath, exist_ok=True)


		for f in files:
			processed.append(client.submit(runScript, root + "/executables/centroidFxn.R", f, outputPath))
		print('Centroiding Negative Mode Masses!')
		progress(processed)
		client.gather(processed)
		client.shutdown()
		stop = time.perf_counter()

		runTime = stop - start

		f = open(root + '/outputs/runTimes/' + runNumber + '/centroidFxnNegative/centroidNegativeOverall.txt', 'w')
		f.writelines(str(runTime))
		f.close()



	#Step6: Subset the Masses

	#Positive Mode TempFile Generation

	if toRun[7] == 1:
		
		if not os.path.exists(root + "/outputs/runTimes/" + runNumber + '/tempFilesPositive/perJob'):
			os.makedirs(root + "/outputs/runTimes/" + runNumber + '/tempFilesPositive/perJob', exist_ok=True)


		logPath = root + '/logs/tempFilesPositive/' + runNumber


		#clear old logs before adding new
		if os.path.exists(logPath):
			for x, y, z in os.walk(logPath, topdown=False):
				for name in z:
					os.remove(os.path.join(x, name))

		if not os.path.exists(logPath):
			os.makedirs(logPath, exist_ok=True)


		start = time.perf_counter()
		root = config['metabolomicsRoot']
		
		#these 2 lines are deprecated
		#files = os.listdir(root + "/outputs/hdf5FilesNewColumn/" + runNumber)
		#files = [root + "/outputs/hdf5FilesNewColumn/" + runNumber + "/"  + f for f in files]		

		#these 2 lines should be in final version
		#rawFiles = os.listdir(root + "/outputs/processedRawFiles/" + runNumber)
		#fileToUse = [root + "/outputs/processedRawFiles/" + runNumber + '/' + f for f in rawFiles]
		

		#for testing:
		#rawFiles = os.listdir(root + "/outputs/processedRawFiles/" + DOE)
		#fileToUse = [root + "/outputs/processedRawFiles/DOE" + '/' + f for f in rawFiles]

		#in final
		#fileToUse = fileToUse[0]


		
		#can manually change this path to choose cutoff manually
		#maxMass = '/home/lconnelly/MetabolomicsNewDeployment/outputs/centroidFxnPositive/broadData/DOE_autocredential_isopairs_polarity_positive_50.txtcentroidMasses.txt'
		#maxMass = '/home/lconnelly/MetabolomicsNewDeployment/outputs/centroidFxnPositive/broadData/DOE_autocredential_isopairs_polarity_positive_40.txtcentroidMasses.txt'
		
		maxMass = root + '/outputs/centroidFxnPositive/' + runNumber + '/DOE_autocredential_isopairs_polarity_positive_' + value + '.txtcentroidMasses.txt'

		#autoselect
		if autoSelect == True:
			runScript(root + "/executables/binSelector.R", fileToUse, str(len(files)), root + '/outputs/centroidFxnPositive/' + runNumber + '/closest')
			path = root + '/outputs/centroidFxnPositive/' + runNumber + '/closest'
			file = open(path, 'r')
			value = file.readline()
			print(value)
			print('VALUE ^^^')
			maxMass = root + '/outputs/centroidFxnPositive/' + runNumber + '/DOE_autocredential_isopairs_polarity_positive_' + value + '.txtcentroidMasses.txt' 
			print(maxMass)
			print("is the centroidMasses file used")


		outputPath = config['metabolomicsRoot'] + '/outputs/tempFilesPositive/' + runNumber

		if not os.path.exists(outputPath):
			os.makedirs(outputPath, exist_ok=True)

		if not os.path.exists(root + '/tempTextFiles'):
			os.makedirs(root + '/tempTextFiles', exist_ok=True)

		
		#maxMass = '/home/lconnelly/MetabolomicsNewDeployment/outputs/centroidFxnNegative/DOE/DOE_autocredential_isopairs_polarity_negative_10.txtcentroidMasses.txt'
		
		#Assign masses in maxmass to chunks
		masses = pd.read_csv(maxMass)


		cluster = HTCondorCluster(cores=1, n_workers=1, memory="8 GB", disk = "6 GB", local_directory="$_CONDOR_SCRATCH_DIR", log_directory = logPath, job_extra={"getenv": True})
		cluster.scale(jobs = 250)
		client = Client(cluster)
		processed = []
		#n = maximum number of masses to give to a worker
		n = 40
		for i in range(65, 1205, 5):
			try:
				slice = masses.x[masses.x >= i-5]
				slice = slice[slice < i]
				list = [slice[i:i+n] for i in range(0, slice.shape[0],n)]
				files = os.listdir(root + '/outputs/hdf5FilesNewColumn/' + runNumber + '/chunk' + str(i) + 'iso')
				files = [root + '/outputs/hdf5FilesNewColumn/' + runNumber + '/chunk' + str(i) + 'iso/' + f for f in files]
				files = ','.join(files)
				fileList = open(root + "/tempTextFiles/fToRun" + str(i) + ".txt", "w")
				fileList.write(files)
				fileList.close()

				for element in list:
					element = element.reset_index(drop = True)
					element = element.to_string(index = False, header = False)
					element = element.replace('\n', ',')
					massesToRun = open(root + "/tempTextFiles/massesToRun" + element[0:10] + ".txt", "w")
					massesToRun.write(element)
					massesToRun.close()				
					processed.append(client.submit(runPythonScript, root + "/executables/create_all_the_temp_files_chunks.py", root + "/tempTextFiles/fToRun" + str(i) + ".txt", root + "/tempTextFiles/massesToRun" + element[0:10] + ".txt", outputPath, "positive"))
				#print('try succeeded for', str(i))
			except Exception as e:
				#if i == 95:
					#print(str(e))
				#print('failed for', str(i))
				continue
		print(len(processed))
		print('Making Positive Mode Masses Temp Files!')
		progress(processed)
		print(processed)
		client.gather(processed)
		client.shutdown()
		
		stop = time.perf_counter()
		runTime = stop - start
		
		f = open(root + '/outputs/runTimes/' + runNumber + '/tempFilesPositive/tempFilesPositiveOverall.txt', 'w')
		f.writelines(str(runTime))
		f.close() 	

		time.sleep(3)
		outputPath = root + '/outputs/runTimes/' + runNumber + '/tempFilesPositive'
		runPythonScript(root + '/executables/parse_logs.py', logPath, outputPath)		

		#for i in range(95, 1205, 5):
			#try:
				#os.remove(root + "/tempTextFiles/fToRun" + str(i) + ".txt")
			#except:
				#continue

	#Negative Mode TempFile Generation
	if toRun[8] == 1:


		if not os.path.exists(root + "/outputs/runTimes/" + runNumber + '/tempFilesNegative/perJob'):
			os.makedirs(root + "/outputs/runTimes/" + runNumber + '/tempFilesNegative/perJob', exist_ok=True)

		logPath = root + '/logs/tempFilesNegative/' + runNumber

		#clear old logs before adding new
		if os.path.exists(logPath):
			for x, y, z in os.walk(logPath, topdown=False):
				for name in z:
					os.remove(os.path.join(x, name))
		
		if not os.path.exists(logPath):
			os.makedirs(logPath, exist_ok=True)


		start = time.perf_counter()
		root = config['metabolomicsRoot']

		#these 2 lines are deprecated
		#files = os.listdir(root + "/outputs/hdf5FilesNewColumn/" + runNumber)
		#files = [root + "/outputs/hdf5FilesNewColumn/"+ runNumber + '/' + f for f in files]
		

		#should be in final version
		#get a raw file to use for bin selection
		rawFiles = os.listdir(root + "/outputs/processedRawFiles/" + runNumber)
		fileToUse = [root + "/outputs/processedRawFiles/" + runNumber + '/' + f for f in rawFiles]
		fileToUse = fileToUse[0]


	        #select which cutoff value to use from autocredential

                #MANUAL SELECTION
		#maxMass = '/home/lconnelly/MetabolomicsNewDeployment/outputs/centroidFxnNegative/ChengZhao_Carnegie/DOE_autocredential_isopairs_polarity_negative_5.txtcentroidMasses.txt'

		maxMass = root + '/outputs/centroidFxnNegative/' + runNumber + '/DOE_autocredential_isopairs_polarity_negative_' + value + '.txtcentroidMasses.txt'

		#auto select
		if autoSelect == True:
			runScript(root + "/executables/binSelector.R", fileToUse, str(len(files)), root + '/outputs/centroidFxnNegative/' + runNumber + '/closest')
			path = root + '/outputs/centroidFxnNegative/' + runNumber + '/closest'
			file = open(path, 'r')
			value = file.readline()
			maxMass = root + '/outputs/centroidFxnNegative/' + runNumber + '/DOE_autocredential_isopairs_polarity_negative_' + value + '.txtcentroidMasses.txt' 
			print(maxMass)
			print("is the centroidMasses file used")


		outputPath = config['metabolomicsRoot'] + '/outputs/tempFilesNegative/' + runNumber

		if not os.path.exists(outputPath):
			os.makedirs(outputPath, exist_ok=True)

		if not os.path.exists(root + '/tempTextFiles'):
			os.makedirs(root + '/tempTextFiles', exist_ok=True)


		masses = pd.read_csv(maxMass)

		cluster = HTCondorCluster(cores=1, n_workers=1, memory="8 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		cluster.scale(jobs = 250)
		client = Client(cluster)
		processed = []

		#n = max # of masses to give to a worker
		n = 40
		for i in range(65, 1205, 5):
			try:
				slice = masses.x[masses.x >= i-5]
				slice = slice[slice < i]
				list = [slice[i:i+n] for i in range(0, slice.shape[0],n)]
				files = os.listdir(root + '/outputs/hdf5FilesNewColumn/' + runNumber + '/chunk' + str(i) + 'iso')
				files = [root + '/outputs/hdf5FilesNewColumn/' + runNumber + '/chunk' + str(i) + 'iso/' + f for f in files]
				files = ','.join(files)
				fileList = open(root + "/tempTextFiles/fToRun" + str(i) + ".txt", "w")
				fileList.write(files)
				fileList.close()
			
				for element in list:
					element = element.reset_index(drop = True)
					element = element.to_string(index = False, header = False)
					element = element.replace('\n', ',')
					massesToRun = open(root + "/tempTextFiles/massesToRun" + element[0:10] + ".txt", "w")
					massesToRun.write(element)
					massesToRun.close()
					processed.append(client.submit(runPythonScript, root + "/executables/create_all_the_temp_files_chunks.py", root + "/tempTextFiles/fToRun" + str(i) + ".txt", root + "/tempTextFiles/massesToRun" + element[0:10] + ".txt", outputPath, "negative"))
			except Exception as e:
				continue

		print('Making Negative Mode Masses Temp Files!')
		progress(processed)
		client.gather(processed)
		client.shutdown()

		stop = time.perf_counter()

		runTime = stop - start

		f = open(root + '/outputs/runTimes/' + runNumber + '/tempFilesNegative/tempFilesNegativeOverall.txt', 'w')
		f.writelines(str(runTime))
		f.close()

		time.sleep(3)
		outputPath = root + '/outputs/runTimes/' + runNumber + '/tempFilesNegative'
		runPythonScript(root + '/executables/parse_logs.py', logPath, outputPath)


	#Step7: Merge files for each mass and generate scan/mass v intensity plots for each


	#Positive Mode Mass Feature
	if toRun[9] == 1:

		if not os.path.exists(root + "/outputs/runTimes/" + runNumber + '/vecOfIntensitiesForEachMassPositive/perJob'):
			os.makedirs(root + "/outputs/runTimes/" + runNumber + '/vecOfIntensitiesForEachMassPositive/perJob', exist_ok=True)

		logPath = root + '/logs/vecOfIntensitiesForEachMassPositive/' + runNumber

		#clear old logs before adding new
		if os.path.exists(logPath):
			for x, y, z in os.walk(logPath, topdown=False):
				for name in z:
					try:
						os.remove(os.path.join(x, name))
					except:
						continue

		if not os.path.exists(logPath):
			os.makedirs(logPath, exist_ok=True)

		




		start = time.perf_counter()
		if not os.path.exists(config['metabolomicsRoot']+"/outputs/outputPdfsPositive/"+ runNumber):
			os.makedirs(config['metabolomicsRoot']+"/outputs/outputPdfsPositive/"+runNumber, exist_ok=True)

		plotPath = config['metabolomicsRoot']+"/outputs/outputPdfsPositive/"+ runNumber

		if not os.path.exists(config['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassPositive/"+runNumber):
			os.makedirs(config['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassPositive/"+runNumber, exist_ok=True)

		vecOfIntensitiesPath = config['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassPositive/"+runNumber

		massDirectories = os.listdir(root + '/outputs/tempFilesPositive/' + runNumber)
		massDirectories = [root + '/outputs/tempFilesPositive/' + runNumber + '/' + m for m in massDirectories]



		print("here is first mass directory:", massDirectories[0])

		cluster = HTCondorCluster(cores=1, n_workers=1, memory="8 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True}, log_directory = logPath)
		cluster.scale(jobs = 250)
		client = Client(cluster)
		processed = []

		for m in massDirectories:
			#NOTE!!!! LINE CHANGED FOR TESTING
			#processed.append(client.submit(runScript, root + "/executables/newestAnovalignChunks.R", m, root + "/outputs/processedRawFiles/DOE", vecOfIntensitiesPath))
			processed.append(client.submit(runScript, root + "/executables/newestAnovalignChunks.R", m, root + "/fileListForAnovalign" + runNumber + ".txt", vecOfIntensitiesPath, 'positive'))
		print('Finding Positive Mode Mass Feature Profiles!')
		progress(processed)
		client.gather(processed)
		client.shutdown()	
		stop = time.perf_counter()
		runTime = stop - start
		
		
		f = open(root + '/outputs/runTimes/' + runNumber + '/vecOfIntensitiesForEachMassPositive/anovalignPositiveOverall.txt', 'w')
		f.writelines(str(runTime))
		f.close()

		time.sleep(3)
		outputPath = root + '/outputs/runTimes/' + runNumber + '/vecOfIntensitiesForEachMassPositive'
		runPythonScript(root + '/executables/parse_logs.py', logPath, outputPath)


	#Negative Mode Mass Feature

	if toRun[10] == 1:

		if not os.path.exists(root + "/outputs/runTimes/" + runNumber + '/vecOfIntensitiesForEachMassNegative/perJob'):
			os.makedirs(root + "/outputs/runTimes/" + runNumber + '/vecOfIntensitiesForEachMassNegative/perJob', exist_ok=True)


		logPath = root + '/logs/vecOfIntensitiesForEachMassNegative/' + runNumber

		#clear old logs before adding new
		if os.path.exists(logPath):
			for x, y, z in os.walk(logPath, topdown=False):
				for name in z:
					try:
						os.remove(os.path.join(x, name))
					except:
						continue

		if not os.path.exists(logPath):
			os.makedirs(logPath, exist_ok=True)


		start = time.perf_counter()
		if not os.path.exists(config['metabolomicsRoot']+"/outputs/outputPdfsNegative/"+runNumber):
			os.makedirs(config['metabolomicsRoot']+"/outputs/outputPdfsNegative/"+runNumber, exist_ok=True)
		
		plotPath = config['metabolomicsRoot']+"/outputs/outputPdfsNegative/"+runNumber

		if not os.path.exists(config['metabolomicsRoot'] + "/outputs/vecOfIntensitiesForEachMassNegative/"+runNumber):
			os.makedirs(config['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassNegative/"+runNumber, exist_ok=True)

		vecOfIntensitiesPath = config['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassNegative/"+runNumber

		massDirectories = os.listdir(root + "/outputs/tempFilesNegative/" + runNumber)
		massDirectories = [root + '/outputs/tempFilesNegative/' + runNumber + '/' +  m for m in massDirectories]

		cluster = HTCondorCluster(cores=1, n_workers=1, memory="8 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", log_directory=logPath,  job_extra={"getenv": True})		
		cluster.scale(jobs = 250)
		client = Client(cluster)
		processed = []
		for m in massDirectories:
			processed.append(client.submit(runScript, root + "/executables/newestAnovalignTempFilesDeveloping.R", m, root + "/fileListForAnovalign" + runNumber + ".txt", vecOfIntensitiesPath, 'negative'))
		print('Finding Negative Mode Mass Feature Profiles!')
		progress(processed)
		client.gather(processed)
		client.shutdown()

		stop = time.perf_counter()
		
		runTime = stop - start

		f = open(root + '/outputs/runTimes/' + runNumber + '/vecOfIntensitiesForEachMassNegative/anovalignNegativeOverall.txt', 'w')
		f.writelines(str(runTime))
		f.close()


		time.sleep(3)
		outputPath = root + '/outputs/runTimes/' + runNumber + '/vecOfIntensitiesForEachMassNegative'
		runPythonScript(root + '/executables/parse_logs.py', logPath, outputPath)
