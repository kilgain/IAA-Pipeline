from dask_jobqueue.htcondor import HTCondorCluster
from dask.distributed import Client, progress, wait
import os
import yaml
import glob

def runScript(x,y,z="", xx = "", xy = ""):
	os.system("Rscript " + x + " " + y + " " + z + " " + xx + " " + xy)


def runPythonScript(x,y = "", z = "", xx = "", xy = ""):
	os.system("python3 " + x + " " + y + " " + z + " " + xx + " " + xy)


#root = directory which contains dask script
root = '/home/lconnelly/MetabolomicsNewDeployment'

#all outputs will be sorted based on this run number
runNumber = 'setariaTest'


toRun = [1,0,0,0,0,0,0,1,0,0]


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


if __name__ == "__main__":


	with open("/home/lconnelly/MetabolomicsNewDeployment/Config.yaml", "r") as file:
		print(type(file))
		config = yaml.safe_load(file)



	#Step1: Add File Column and Collapse Duplicate Masses
	if toRun[0] == 1:
		#initiates job cluster (for Condor)
		cluster = HTCondorCluster(n_workers = 1, cores=1, memory = "4 GB", disk = "4 GB",local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		with open("/home/lconnelly/MetabolomicsNewDeployment/Config.yaml", "r") as file:
			config = yaml.safe_load(file)
		initConfig = config['samplesDirectory']  + '/' +  runNumber
		print(initConfig, "is config path")
		files = os.listdir(initConfig)
		#sets job cluster size equal to number of files
		cluster.scale(jobs = len(files))
		#initiates job client that assigns jobs using created cluster
		client = Client(cluster)
		#list of job futures
		processed = []
		f = open(initConfig + '/' + files[0], 'r')
		line1 = f.readline()
		f.close()		
		files = [initConfig + '/' + f for f in files]

		if not '"File"' in line1:
		#submits all files for formatting
			for f in files:
				processed.append(client.submit(runScript, root + "/executables/fileAdd.R", f))
		#progress bar
		progress(processed)
		#prevents progress until all jobs are done
		client.gather(processed)
		#destroys job cluster
		client.shutdown()


	#Step2: Convert .txt Files to HDF5 Format

	if toRun[1] == 1:

		cluster = HTCondorCluster(cores=1, n_workers=1, memory = "4 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		cluster.scale(jobs = len(files))
		client = Client(cluster)
		processed = []


		outputPath = root + "/outputs/hdf5Files/" + runNumber

		if not os.path.exists(outputPath):
			os.makedirs(outputPath, exist_ok=True)

	
		print("is input for conversion")
		print(files[0])
	
		print("is output directory")
		print(outputPath)

		for f in files:
			processed.append(client.submit(runPythonScript, root + "/executables/create_the_hdf5_files.py", f, outputPath))
		progress(processed)
		client.gather(processed)
		client.shutdown()



	#Step3: Run Isolock On Each File
	if toRun[2] == 1:
		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		cluster.scale(jobs = len(files))
		client = Client(cluster)
		processed = []

		files = os.listdir(root + '/outputs/hdf5Files/' + runNumber)
		files = [root + '/outputs/hdf5Files/' + runNumber + '/' + f for f in files]
		refFile = files[0]

		print(refFile)
		print("is refFile")




		outputPath = root + "/outputs/hdf5FilesNewColumn/" + runNumber
		if not os.path.exists(outputPath):
                	os.makedirs(outputPath, exist_ok=True)

	

		for f in files:
			processed.append(client.submit(runPythonScript, root + "/executables/python_script_determine_offset.py", f, refFile))
		progress(processed)
		client.gather(processed)
		client.shutdown()

	#Step4: Run Autocredential



	#Positive:

	if toRun[3] == 1:

		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv":True})
		cluster.scale(jobs = 1)
		client = Client(cluster)
		processed = []


		#path = /home/lconnelly/Metabolomics/rawDataFiles/Analysis1
		path = config['samplesDirectory'] + '/' +  runNumber
		#path = /home/lconnelly/Metabolomics/outputs/hdf5FilesNewColumn/Analysis1
		path = path.replace("rawDataFiles", "outputs/hdf5FilesNewColumn")
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

		processed.append(client.submit(runPythonScript, root + "/executables/autocredential.py", "positive", path, outputPathPositive, outputPathNegative))

		progress(processed)
		client.gather(processed)
		client.shutdown()

	#Negative:

	if toRun[4] == 1:


		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv":True})
		cluster.scale(jobs = 1)
		client = Client(cluster)
		processed = []


		outputPathPositive = root + "/outputs/autocredentialPositive/" + runNumber
		outputPathNegative = root + "/outputs/autocredentialNegative/" + runNumber


		if not os.path.exists(outputPathPositive):
			os.makedirs(outputPathPositive, exist_ok=True)


		if not os.path.exists(outputPathNegative):
			os.makedirs(outputPathNegative, exist_ok=True)


		#path = /home/lconnelly/Metabolomics/rawDataFiles/Analysis1
		path = config['samplesDirectory'] + '/' + runNumber
                #path = /home/lconnelly/Metabolomics/outputs/hdf5FilesNewColumn/Analysis1
		path = path.replace("rawDataFiles", "outputs/hdf5FilesNewColumn")
		inputFiles = os.listdir(path)
		inputFiles = [path + '/' + f for f in inputFiles]

		processed.append(client.submit(runPythonScript, root + "/executables/autocredential.py", "negative", path, outputPathPositive, outputPathNegative))

		progress(processed)
		client.gather(processed)
		client.shutdown()

	#Step5: Centroid the Masses
	if toRun[5] == 1:
		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		cluster.scale(jobs = 10)
		client = Client(cluster)
		processed = []

	#Positive:

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

		progress(processed)
		client.gather(processed)
		client.shutdown()

	#Negative:
	
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

		progress(processed)
		client.gather(processed)
		client.shutdown()


	#Step6: Subset the Masses

	if toRun[7] == 1:
	#positive:

		root = config['metabolomicsRoot']
		files = os.listdir(root + "/outputs/hdf5FilesNewColumn/" + runNumber)
		files = [root + "/outputs/hdf5FilesNewColumn/"+ runNumber + '/' + f for f in files]

		print(files[0])

		#select which cutoff value to use from autocredential
		centroidMasses = config['metabolomicsRoot']
		centroidMasses = centroidMasses + "/outputs/centroidFxnPositive/" + runNumber
		massFiles = os.listdir(centroidMasses)
		massFiles = [centroidMasses + "/" + s for s in massFiles]	
		maxMass = max(massFiles, key = lambda x: os.stat(x).st_size)
	
		print(maxMass)
		print("is the centroidMasses file used")


		outputPath = config['metabolomicsRoot'] + '/outputs/tempFilesPositive/' + runNumber

		if not os.path.exists(outputPath):
			os.makedirs(outputPath, exist_ok=True)


		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
		cluster.scale(jobs = 250)
		client = Client(cluster)
		processed = []
		for f in files:
			processed.append(client.submit(runPythonScript, root + "/executables/create_all_the_temp_files.py", f, maxMass, outputPath, "Positive"))
	
		progress(processed)
		client.gather(processed)
		client.shutdown()
	

	#Negative:
	

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

		progress(processed)
		client.gather(processed)
		client.shutdown()


	#Step7: Merge files for each mass and generate scan/mass v intensity plots for each


	#Positive
	if not os.path.exists(config['metabolomicsRoot']+"/outputs/outputPdfsPositive/"+ runNumber):
		os.makedirs(config['metabolomicsRoot']+"/outputs/outputPdfsPositive/"+runNumber, exist_ok=True)

	plotPath = config['metabolomicsRoot']+"/outputs/outputPdfsPositive/"+ runNumber

	if not os.path.exists(config['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassPositive/"+runNumber):
		os.makedirs(config['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassPositive/"+runNumber, exist_ok=True)

	vecOfIntensitiesPath = config['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassPositive/"+runNumber

	massDirectories = os.listdir(root + '/outputs/tempFilesPositive/' + runNumber)
	massDirectories = [root + '/outputs/tempFilesPositive/' + runNumber + '/' + m for m in massDirectories]




	cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "4 GB", local_directory="$_CONDOR_SCRATCH_DIR", job_extra={"getenv": True})
	cluster.scale(jobs = len(massDirectories))
	client = Client(cluster)
	processed = []

	if toRun[9] == 1:
		for m in massDirectories:
			processed.append(client.submit(runScript, root + "/executables/read_in_the_masses_tempfiles.R", m, config['fileList'], plotPath, vecOfIntensitiesPath))

	progress(processed)
	client.gather(processed)
	client.shutdown()	

