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
#example:
root = '/home/lconnelly/MetabolomicsNewDeployment'

#all outputs will be sorted based on this run number
runNumber = 'DOE'

#select which steps to run
#	 0 1 2 3 4 5 6 7 8 9 10
toRun = [1,0,1,0,0,0,0,0,0,0,0]


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

#Converts .raw files to text files (skip if raw files have already been processed)
runRawFileProcessing = False

#Determine noise cutoff automatically, or manually input a value (multiples of 10. 10 is enough for small sample sizes. Manual used if autoSelect is False.)
autoSelect = False
value = '40'





if __name__ == "__main__":

	if not os.path.exists(root + "/outputs/runTimes/" + runNumber):
		os.makedirs(root + "/outputs/runTimes/" + runNumber, exist_ok=True)


	with open("/home/lconnelly/MetabolomicsNewDeployment/Config.yaml", "r") as file:
		config = yaml.safe_load(file)
	if runRawFileProcessing:
		start = time.perf_counter()
		with open("/home/lconnelly/MetabolomicsNewDeployment/Config.yaml", "r") as file:
			config = yaml.safe_load(file)

		directoryOfRawFiles = root + '/rawDataFiles/' + runNumber
		print(directoryOfRawFiles, 'is the raw file directory')

		#should be a list  of paths to each individual .raw file
		rawFiles = [os.path.join(directoryOfRawFiles, file) for file in os.listdir(directoryOfRawFiles)]		


		print(rawFiles[0], 'is the first rawFile')
		outputDirectory = root + '/outputs/processedRawFiles/' + runNumber
		print(outputDirectory, 'is where the processed .raw files will be located')
		

		#check for any previous outputs, dont rerun them
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

		#create output directory
		if not os.path.exists(outputDirectory):
			os.makedirs(outputDirectory, exist_ok=True)
		print('number of raw files to be processed:', len(rawFiles))
		cluster = HTCondorCluster(n_workers = 1, cores=1, memory = "8 GB", disk = "4 GB",local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		cluster.scale(jobs = 250)
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

		f = open(root + '/outputs/runTimes/' + runNumber + '/rawFileProcessing.txt', 'w')
		f.writelines(str(runTime))
		f.close()



	#Step1: Add File Column and Collapse Duplicate Masses
	if toRun[0]:
		#initiates job cluster (for Condor), these are the resources per worker
		cluster = HTCondorCluster(n_workers = 1, cores=1, memory = "8 GB", disk = "4 GB",local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		with open("/home/lconnelly/MetabolomicsNewDeployment/Config.yaml", "r") as file:
			config = yaml.safe_load(file)
		initConfig = root  + '/outputs/processedRawFiles/' +  runNumber
		print(initConfig, "is config path")
		files = os.listdir(initConfig)
		#sets number of workers
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
		
		start = time.perf_counter()

		cluster = HTCondorCluster(cores=1, n_workers=1, memory = "4 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		cluster.scale(jobs = 300)
		client = Client(cluster)
		processed = []


		initConfig = root  + '/outputs/processedRawFiles/' +  runNumber
		outputPath = root + "/outputs/hdf5Files/" + runNumber


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

		print(len(files), 'is length of input for conversion after checking for old outputs')

		for f in files:
			processed.append(client.submit(runPythonScript, root + "/executables/create_the_hdf5_files.py", f, outputPath))
		print('Converting Files To Hdf5!')
		print(len(processed))
		progress(processed)
		client.gather(processed)
		client.shutdown()


		stop = time.perf_counter()

		runTime = stop - start

		f = open(root + '/outputs/runTimes/' + runNumber + '/HDF5conversion.txt', 'w')
		f.writelines(str(runTime))
		f.close()

		

	#Step3: Run Isolock On Each File
	if toRun[2] == 1:

		start = time.perf_counter()
		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		cluster.scale(jobs = 300)
		client = Client(cluster)
		processed = []

		files = os.listdir(root + '/outputs/hdf5Files/' + runNumber)
		files = [root + '/outputs/hdf5Files/' + runNumber + '/' + f for f in files]
		refFile = files[0]
		#can manually changed the reference file used here:
		#refFile = '/path'


		#the reference file used will be stored in a text file in the root directory
		#refFileForIsolockRunNumber.txt
		outFileLocation = root + '/refFileForIsolock' + runNumber + '.txt'
		with open(outFileLocation, "w") as fileOut:
			fileOut.write(refFile)
		fileOut.close()
		
		print(refFile)
		print("^ is refFile ^")




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
		
		f = open(root + '/outputs/runTimes/' + runNumber + '/isolockRuntime.txt', 'w')
		f.writelines(str(runTime))
		f.close()

	#Step4: Run Autocredential

	#Positive Mode Autocredential

	if toRun[3] == 1:




		start = time.perf_counter()
		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "10 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv":True})
		cluster.scale(jobs = 250)
		client = Client(cluster)
		processed = []


		path = root + '/outputs/hdf5FilesNewColumn/' + runNumber


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
		

		#output is in chunks, need to stitch together chunks
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
		
		f = open(root + '/outputs/runTimes/' + runNumber + '/autocredentialPositive.txt', 'w')
		f.writelines(str(runTime))
		f.close()

	#Negative Mode Autocredential

	if toRun[4] == 1:

		start = time.perf_counter()
		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "10 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv":True})
		cluster.scale(jobs = 250)
		client = Client(cluster)
		processed = []

		path = root + '/outputs/hdf5FilesNewColumn/' + runNumber
		inputFiles = os.listdir(path)
		inputFiles = [path + '/' + f for f in inputFiles]


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

		f = open(root + '/outputs/runTimes/' + runNumber + '/autocredentialNegative.txt', 'w')
		f.writelines(str(runTime))
		f.close()


	#Step5: Centroid the Masses

	#Positive Mode Centroid

	if toRun[5] == 1:

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

		f = open(root + '/outputs/runTimes/' + runNumber + '/centroidPositive.txt', 'w')
		f.writelines(str(runTime))
		f.close()

	#Negative Mode Centroid:
	
	if toRun[6] == 1:


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

		f = open(root + '/outputs/runTimes/' + runNumber + '/centroidNegative.txt', 'w')
		f.writelines(str(runTime))
		f.close()



	#Step6: Subset the Masses

	#Positive Mode TempFile Generation

	if toRun[7] == 1:
		
		start = time.perf_counter()
		root = config['metabolomicsRoot']

		
		#used for auto selection of noise cutoff
		rawFiles = os.listdir(root + "/outputs/processedRawFiles/" + runNumber)
		fileToUse = [root + "/outputs/processedRawFiles/" + runNumber + '/' + f for f in rawFiles]
		fileToUse = fileToUse[0]


		#can manually change this path to choose cutoff manually
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

		
		#Assign masses in maxmass to chunks
		masses = pd.read_csv(maxMass)


		cluster = HTCondorCluster(cores=1, n_workers=1, memory="8 GB", disk = "6 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		cluster.scale(jobs = 300)
		client = Client(cluster)
		processed = []
		#n = maximum number of masses to give to a worker
		n = 50
		for i in range(65, 1205, 5):
			try:
				slice = masses.x[masses.x >= i-6.1]
				slice = slice[slice < i+1.1]
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
					processed.append(client.submit(runPythonScript, root + "/executables/create_all_the_temp_files_chunks.py", root + "/tempTextFiles/fToRun" + str(i) + ".txt", root + "/tempTextFiles/massesToRun" + element[0:10] + ".txt", outputPath))
				#print('try succeeded for', str(i))
			except Exception as e:
				continue
		print('Making Positive Mode Masses Temp Files!')
		progress(processed)
		print(processed)
		client.gather(processed)
		client.shutdown()
		
		stop = time.perf_counter()
		runTime = stop - start
		
		f = open(root + '/outputs/runTimes/' + runNumber + '/tempFilesPositive.txt', 'w')
		f.writelines(str(runTime))
		f.close() 	


	#Negative Mode TempFile Generation
	

	if toRun[8] == 1:

		start = time.perf_counter()
		root = config['metabolomicsRoot']
		files = os.listdir(root + "/outputs/hdf5FilesNewColumn/" + runNumber)
		files = [root + "/outputs/hdf5FilesNewColumn/"+ runNumber + '/' + f for f in files]
		
		#get a raw file to use for bin selection
		rawFiles = os.listdir(root + "/outputs/processedRawFiles/" + runNumber)
		fileToUse = [root + "/outputs/processedRawFiles/" + runNumber + '/' + f for f in rawFiles]
		fileToUse = fileToUse[0]


	        #select which cutoff value to use from autocredential

                #MANUAL SELECTION
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
		cluster.scale(jobs = 300)
		client = Client(cluster)
		processed = []

		#n = max # of masses to give to a worker
		n = 50
		for i in range(65, 1205, 5):
			try:
				slice = masses.x[masses.x >= i-6.1]
				slice = slice[slice < i+1.1]
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
					processed.append(client.submit(runPythonScript, root + "/executables/create_all_the_temp_files_chunks.py", root + "/tempTextFiles/fToRun" + str(i) + ".txt", root + "/tempTextFiles/massesToRun" + element[0:10] + ".txt", outputPath))
			except Exception as e:
				continue

		print('Making Negative Mode Masses Temp Files!')
		progress(processed)
		client.gather(processed)
		client.shutdown()

		stop = time.perf_counter()

		runTime = stop - start

		f = open(root + '/outputs/runTimes/' + runNumber + '/tempFilesNegative.txt', 'w')
		f.writelines(str(runTime))
		f.close()


	#Step7: Merge files for each mass and generate scan/mass v intensity plots for each


	#Positive Mode Mass Feature

	if toRun[9] == 1:

		start = time.perf_counter()

		if not os.path.exists(config['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassPositive/"+runNumber):
			os.makedirs(config['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassPositive/"+runNumber, exist_ok=True)

		vecOfIntensitiesPath = config['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassPositive/"+runNumber

		massDirectories = os.listdir(root + '/outputs/tempFilesPositive/' + runNumber)
		massDirectories = [root + '/outputs/tempFilesPositive/' + runNumber + '/' + m for m in massDirectories]



		print("here is first mass directory:", massDirectories[0])

		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		cluster.scale(jobs = 300)
		client = Client(cluster)
		processed = []

		for m in massDirectories:
			processed.append(client.submit(runScript, root + "/executables/newestAnovalignChunks.R", m, root + "/outputs/processedRawFiles/" + runNumber, vecOfIntensitiesPath))
		print('Finding Positive Mode Mass Feature Profiles!')
		progress(processed)
		client.gather(processed)
		client.shutdown()	
		stop = time.perf_counter()
		runTime = stop - start
		
		
		f = open(root + '/outputs/runTimes/' + runNumber + '/anovalignPositive.txt', 'w')
		f.writelines(str(runTime))
		f.close()



	#Negative Mode Mass Feature

	if toRun[10] == 1:

		start = time.perf_counter()

		if not os.path.exists(config['metabolomicsRoot'] + "/outputs/vecOfIntensitiesForEachMassNegative/"+runNumber):
			os.makedirs(config['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassNegative/"+runNumber, exist_ok=True)

		vecOfIntensitiesPath = config['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassNegative/"+runNumber

		massDirectories = os.listdir(root + "/outputs/tempFilesNegative/" + runNumber)
		massDirectories = [root + '/outputs/tempFilesNegative/' + runNumber + '/' +  m for m in massDirectories]

		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})		
		cluster.scale(jobs = 250)
		client = Client(cluster)
		processed = []
		for m in massDirectories:
			processed.append(client.submit(runScript, root + "/executables/newestAnovalignChunks.R", m, root + "/outputs/processedRawFiles/" + runNumber, vecOfIntensitiesPath))
		print('Finding Negative Mode Mass Feature Profiles!')
		progress(processed)
		client.gather(processed)
		client.shutdown()

		stop = time.perf_counter()
		
		runTime = stop - start

		f = open(root + '/outputs/runTimes/' + runNumber + '/anovalignNegative.txt', 'w')
		f.writelines(str(runTime))
		f.close()
