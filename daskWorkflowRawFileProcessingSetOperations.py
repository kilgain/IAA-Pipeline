from dask_jobqueue.htcondor import HTCondorCluster
from dask.distributed import Client, progress, wait
import os
import yaml
import glob
from contextlib import suppress
import pandas as pd

def runScript(x,y,z="", xx = "", xy = ""):
	os.system("Rscript " + x + " " + y + " " + z + " " + xx + " " + xy)


def runPythonScript(x,y = "", z = "", xx = "", xy = ""):
	os.system("python3 " + x + " " + y + " " + z + " " + xx + " " + xy)

def runMono(x,y,z='', xx=''):
	os.system("mono " + x + " " + y + " " + z + " " + xx)


#root = directory which contains dask script
root = '/home/lconnelly/MetabolomicsNewDeployment'

#all outputs will be sorted based on this run number
runNumber = 'testFullPipelineIsoNew'

#	 0 1 2 3 4 5 6 7 8 9 10
toRun = [1,1,1,1,0,1,0,1,0,1,0]
runRawFileProcessing = True

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


	with open("/home/lconnelly/MetabolomicsNewDeployment/Config.yaml", "r") as file:
		print(type(file))
		config = yaml.safe_load(file)
	if runRawFileProcessing:
		with open("/home/lconnelly/MetabolomicsNewDeployment/Config.yaml", "r") as file:
			config = yaml.safe_load(file)
		directoryOfRawFiles = '/home/ahubbard/all_setaria_and_sorghum_raw_files_to_backup/Setaria_Sorghum_Metabolomics'
		
		#directoryOfRawFiles = root + '/rawDataFiles/' + runNumber
		print(directoryOfRawFiles, 'is the raw file directory')
		rawFiles = glob.glob(directoryOfRawFiles+'/*/*HILIC*/*.raw', recursive = True)
		#TESTING:
		rawFiles = rawFiles[0:6]
		print('HERE IS LENGTH OF RAWFILES:', len(rawFiles))
		####
		print(rawFiles[0], 'is the first rawFile')
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
		cluster.scale(jobs = 10)
		client = Client(cluster)
		processed = []
		for f in rawFiles:
			processed.append(client.submit(runMono, root + "/executables/fileProcessingMono/fileProcessor.exe", f, outputDirectory))
		print('Processing Raw Files!')
		progress(processed)
		client.gather(processed)
		client.shutdown()



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
		cluster.scale(jobs = 10)
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

		cluster = HTCondorCluster(cores=1, n_workers=1, memory = "4 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		cluster.scale(jobs = 10)
		client = Client(cluster)
		processed = []


		initConfig = root  + '/outputs/processedRawFiles/' +  runNumber
		outputPath = root + "/outputs/hdf5Files/" + runNumber


		if not os.path.exists(outputPath):
			os.makedirs(outputPath, exist_ok=True)

	
		#initConfig = root  + '/outputs/processedRawFiles/' +  runNumber
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
		progress(processed)
		client.gather(processed)
		client.shutdown()



	#Step3: Run Isolock On Each File
	if toRun[2] == 1:
		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True}, log_directory='/home/lconnelly/MetabolomicsNewDeployment/logs')
		cluster.scale(jobs = 10)
		client = Client(cluster)
		processed = []

		files = os.listdir(root + '/outputs/hdf5Files/' + runNumber)
		files = [root + '/outputs/hdf5Files/' + runNumber + '/' + f for f in files]
		refFile = files[0]
		#refFile = '/home/lconnelly/MetabolomicsNewDeployment/outputs/hdf5Files/DOE/HILIC_Set_1342_3_3_H11_85_TB_setaria_12_0127.rawoutput.txt.hdf5'


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
			processed.append(client.submit(runPythonScript, root + "/executables/isolock.py", f, refFile))
		print('Running Isolock!')
		progress(processed)
		client.gather(processed)
		client.shutdown()

	#Step4: Run Autocredential

	#Positive Mode Autocredential

	if toRun[3] == 1:

		cluster = HTCondorCluster(cores=1, n_workers=1, memory="100 GB", disk = "10 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv":True})
		cluster.scale(jobs = 1)
		client = Client(cluster)
		processed = []


		path = root + '/outputs/hdf5FilesNewColumn/' + runNumber
		inputFiles = os.listdir(path)
		inputFiles = [path + '/' + f for f in inputFiles]

		print(inputFiles[0])
		print('is first input for autocred')


		outputPathPositive = root + "/outputs/autocredentialPositive/" + runNumber
		outputPathNegative = root + "/outputs/autocredentialNegative/" + runNumber

		if not os.path.exists(outputPathPositive):
			os.makedirs(outputPathPositive, exist_ok=True)


		if not os.path.exists(outputPathNegative):
			os.makedirs(outputPathNegative, exist_ok=True)

		processed.append(client.submit(runPythonScript, root + "/executables/autocredentialVaex.py", "positive", path, outputPathPositive, outputPathNegative))
		print('Running Autocredential Positive Mode!')
		progress(processed)
		client.gather(processed)
		client.shutdown()

	#Negative Mode Autocredential

	if toRun[4] == 1:


		cluster = HTCondorCluster(cores=1, n_workers=1, memory="50 GB", disk = "10 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv":True})
		cluster.scale(jobs = 1)
		client = Client(cluster)
		processed = []


		outputPathPositive = root + "/outputs/autocredentialPositive/" + runNumber
		outputPathNegative = root + "/outputs/autocredentialNegative/" + runNumber


		if not os.path.exists(outputPathPositive):
			os.makedirs(outputPathPositive, exist_ok=True)


		if not os.path.exists(outputPathNegative):
			os.makedirs(outputPathNegative, exist_ok=True)


		path = root + '/outputs/hdf5FilesNewColumn/' + runNumber
		inputFiles = os.listdir(path)
		inputFiles = [path + '/' + f for f in inputFiles]

		processed.append(client.submit(runPythonScript, root + "/executables/autocredentialVaex.py", "negative", path, outputPathPositive, outputPathNegative))
		print('Running Autocredential Negative Mode!')
		progress(processed)
		client.gather(processed)
		client.shutdown()

	#Step5: Centroid the Masses

	#Positive Mode Centroid

	if toRun[5] == 1:
		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
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

	#Negative Mode Centroid:
	
	if toRun[6] == 1:

		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		cluster.scale(jobs = len(files))
		client = Client(cluster)
		processed = []


		path = config['metabolomicsRoot'] + "/outputs/autocredentialNegative/" + runNumber
		files = os.listdir(path)
		files = [path + '/' + f for f in files]
		outputPath = config['metabolomicsRoot'] + '/outputs/centroidFxnNegative/' + runNumber

		if not os.path.exists(outputPath):
			os.makedirs(outputPath, exist_ok=True)


		for f in files:
			processed.append(client.submit(runScript, root + "/executables/centroidFxn.R", f, outputPath))
		print('Centroiding Negative Mode Masses!')
		progress(processed)
		client.gather(processed)
		client.shutdown()


	#Step6: Subset the Masses

	#Positive Mode TempFile Generation

	if toRun[7] == 1:

		root = config['metabolomicsRoot']
		files = os.listdir(root + "/outputs/hdf5FilesNewColumn/" + runNumber)
		files = [root + "/outputs/hdf5FilesNewColumn/DOE/" + f for f in files]
		

		print(files[0])

		#select which cutoff value to use from autocredential
		centroidMasses = config['metabolomicsRoot']
		centroidMasses = centroidMasses + "/outputs/centroidFxnPositive/" + runNumber
		massFiles = os.listdir(centroidMasses)
		massFiles = [centroidMasses + "/" + s for s in massFiles]	
		maxMass = max(massFiles, key = lambda x: os.stat(x).st_size)
	
		#can manually change this path to choose cutoff manually
		#maxMass = '/home/lconnelly/MetabolomicsNewDeployment/outputs/centroidFxnPositive/broadData/DOE_autocredential_isopairs_polarity_positive_50.txtcentroidMasses.txt'
		
		print(maxMass)
		print("is the centroidMasses file used")


		outputPath = config['metabolomicsRoot'] + '/outputs/tempFilesPositive/' + runNumber

		if not os.path.exists(outputPath):
			os.makedirs(outputPath, exist_ok=True)


		cluster = HTCondorCluster(cores=1, n_workers=1, memory="30 GB", disk = "6 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		cluster.scale(jobs = 250)
		client = Client(cluster)
		processed = []
		for f in files:
			processed.append(client.submit(runPythonScript, root + "/executables/create_all_the_temp_files.py", f, maxMass, outputPath, "Positive"))
		print('Making Positive Mode Masses Temp Files!')
		progress(processed)
		client.gather(processed)
		client.shutdown()
	

	#Negative Mode TempFile Generation
	

	if toRun[8] == 1:
		root = config['metabolomicsRoot']
		files = os.listdir(root + "/outputs/hdf5FilesNewColumn/" + runNumber)
		files = [root + "/outputs/hdf5FilesNewColumn/"+ runNumber + '/' + f for f in files]


	        #select which cutoff value to use from autocredential
		centroidMasses = config['metabolomicsRoot']
		centroidMasses = centroidMasses + "/outputs/centroidFxnNegative/" + runNumber
		massFiles = os.listdir(centroidMasses)
		massFiles = [centroidMasses + "/" + s for s in massFiles]
		maxMass = max(massFiles, key = lambda x: os.stat(x).st_size)

		print(maxMass)
		print("is the centroidMasses file used")



		outputPath = config['metabolomicsRoot'] + '/outputs/tempFilesNegative/' + runNumber

		if not os.path.exists(outputPath):
			os.makedirs(outputPath, exist_ok=True)



		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		cluster.scale(jobs = len(files))
		client = Client(cluster)
		processed = []
		for f in files:
			processed.append(client.submit(runPythonScript, root + "/executables/create_all_the_temp_files.py", f, maxMass, outputPath, "Negative"))
		print('Making Negative Mode Masses Temp Files!')
		progress(processed)
		client.gather(processed)
		client.shutdown()


	#Step7: Merge files for each mass and generate scan/mass v intensity plots for each


	#Positive Mode Mass Feature

	if toRun[9] == 1:
		if not os.path.exists(config['metabolomicsRoot']+"/outputs/outputPdfsPositive/"+ runNumber):
			os.makedirs(config['metabolomicsRoot']+"/outputs/outputPdfsPositive/"+runNumber, exist_ok=True)

		plotPath = config['metabolomicsRoot']+"/outputs/outputPdfsPositive/"+ runNumber

		if not os.path.exists(config['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassPositive/"+runNumber):
			os.makedirs(config['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassPositive/"+runNumber, exist_ok=True)

		vecOfIntensitiesPath = config['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassPositive/"+runNumber

		massDirectories = os.listdir(root + '/outputs/tempFilesPositive/' + runNumber)
		massDirectories = [root + '/outputs/tempFilesPositive/' + runNumber + '/' + m for m in massDirectories]



		print("here is first mass directory:", massDirectories[0])

		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		cluster.scale(jobs = 250)
		client = Client(cluster)
		processed = []

		for m in massDirectories:
			processed.append(client.submit(runScript, root + "/executables/newestAnovalign.R", m, config['fileList'], plotPath, vecOfIntensitiesPath))
		print('Finding Positive Mode Mass Feature Profiles!')
		progress(processed)
		client.gather(processed)
		client.shutdown()	


	#Negative Mode Mass Feature

	if toRun[10] == 1:
		if not os.path.exists(config['metabolomicsRoot']+"/outputs/outputPdfsNegative/"+runNumber):
			os.makedirs(config['metabolomicsRoot']+"/outputs/outputPdfsNegative/"+runNumber, exist_ok=True)
		
		plotPath = config['metabolomisRoot']+"/outputs/outputPdfsNegative/"+runNumber

		if not os.path.exists(config['metabolomicsRoot'] + "/outputs/vecOfIntensitiesForEachMassNegative/"+runNumber):
			os.makedirs(config['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassNegative/"+runNumber, exist_ok=True)

		vecOfIntensitiesPath = config['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassNegative/"+runNumber

		massDirectories = os.listdir(root + "/outputs/tempFilesNegative/" + runNumber)
		massDirectories = [root + '/outputs/tempFilesNegative/' + runNumber + m for m in massDirectories]

		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})		
		cluster.scale(jobs = len(massDirectories))
		client = Client(cluster)
		processed = []
		for m in massDirectories:
			processed.append(client.submit(runScript, root + "/executables/newestAnovalign.R", m, config['fileList'], plotPath, vecOfIntensitiesPath))
		print('Finding Negative Mode Mass Feature Profiles!')
		progress(processed)
		client.gather(processed)
		client.shutdown()
