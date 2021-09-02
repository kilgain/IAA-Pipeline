from dask_jobqueue.htcondor import HTCondorCluster
from dask.distributed import Client, progress
import os
import yaml
import glob

def runScript(x,y,z="", xx = "", xy = ""):
	os.system("Rscript " + x + " " + y + " " + z + " " + xx + " " + xy)


def runPythonScript(x,y = "", z = "", xx = "", xy = ""):
	os.system("python3 " + x + " " + y + " " + z + " " + xx + " " + xy)

root = '/home/lconnelly/Metabolomics'
runNumber = '300Standards'
toRun = [0,0,0,0,0,0,0,0,0,1]
toRunInit = 0


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


	with open("/home/lconnelly/Metabolomics/Config.yaml") as file:
		yaml = yaml.safe_load(file)



	#Step1: Add File Column and Collapse Duplicate Masses
	if toRunInit == 1:
		cluster = HTCondorCluster(n_workers = 1, cores=1, memory = "4 GB", disk = "4 GB", local_directory = root)
		with open("/home/lconnelly/Metabolomics/Config.yaml") as file:
			config = yaml.safe_load(file)
	#config = /home/lconnelly/Metabolomics/rawDataFiles/Analysis1
		config = config['samplesDirectory']  + '/' +  runNumber
		files = os.listdir(config)
		cluster.scale(jobs = len(files))
		client = Client(cluster)
	#list of job futures
		processed = []
	#for loop of submissions
		f = open(config + '/' + files[0], 'r')
		line1 = f.readline()
		f.close()		
		files = [config + '/' + f for f in files]

		if not '"File"' in line1:
			for f in files:
				processed.append(client.submit(runScript, root + "/executables/fileAdd.R", f))
		progress(processed)
		client.gather(processed)
		client.shutdown()


	#Step2: Convert .txt Files to HDF5 Format

	if toRun[1] == 1:

		cluster = HTCondorCluster(cores=1, n_workers=1, memory = "4 GB", disk = "4 GB", local_directory = root)
		cluster.scale(jobs = len(files))
		client = Client(cluster)
		processed = []



	#config = /home/lconnelly/Metabolomics/rawDataFiles/Analysis1
		config = config.replace("rawDataFiles", "outputs/hdf5Files")

	#config = /home/lconnelly/Metabolomics/outputs/hdf5Files/Analysis1
		if not os.path.exists(config):
			os.makedirs(config, exist_ok=True)


		outputPath = config

		with open("/home/lconnelly/Metabolomics/Config.yaml") as file:
			yaml = yaml.safe_load(file)



	
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
		cluster = HTCondorCluster(cores=1, n_workers = 1, memory = "4 GB", disk = "4 GB", local_directory = root)
		cluster.scale(jobs = len(files))
		client = Client(cluster)
		processed = []

	#gets basename of paths

	#config = /home/lconnelly/Metabolomics/outputs/hdf5Files/Analysis1
		files = os.listdir(config)
	#need full path of files here
	#/home/lconnelly/Metabolomics/outputs/hdf5Files/Analysis1/file1.txt.hdf5
		files = [config + '/' + f for f in files]
		refFile = files[0]

		print(refFile)
		print("is refFile")



	#config = /home/lconnelly/Metabolomics/outputs/hdf5Files/Analysis1
	
		config = config.replace("hdf5Files", "hdf5FilesNewColumn")



		print(config)
		print("is output path for isolock")
	#config = /home/lconnelly/Metabolomics/outputs/hdf5FilesNewColumn/Analysis1

		if not os.path.exists(config):
                	os.makedirs(config, exist_ok=True)

	

		for f in files:
			processed.append(client.submit(runPythonScript, root + "/executables/python_script_determine_offset.py", f, refFile))
		progress(processed)
		client.gather(processed)
		client.shutdown()

	#Step4: Run Autocredential


	cluster = HTCondorCluster(cores=1, n_workers = 1, memory = "20 GB", disk = "4 GB", local_directory = root)
	cluster.scale(jobs = 1)
	client = Client(cluster)
	processed = []



	if yaml['polarity'] != "negative":
	#Positive:

		if toRun[3] == 1:



		#path = /home/lconnelly/Metabolomics/rawDataFiles/Analysis1
			path = yaml['samplesDirectory'] + '/' +  runNumber
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


	if yaml['polarity'] != "positive":
	#Negative:

		if toRun[4] == 1:
			outputPathPositive = root + "/outputs/autocredentialPositive/" + runNumber
			outputPathNegative = root + "/outputs/autocredentialNegative/" + runNumber


			if not os.path.exists(outputPathPositive):
				os.makedirs(outputPathPositive, exist_ok=True)


			if not os.path.exists(outputPathNegative):
				os.makedirs(outputPathNegative, exist_ok=True)




		#path = /home/lconnelly/Metabolomics/rawDataFiles/Analysis1
			path = yaml['samplesDirectory'] + '/' + runNumber
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
		cluster = HTCondorCluster(cores=1, n_workers = 1, memory = "4 GB", disk = "4 GB", local_directory = root)
		cluster.scale(jobs = 10)
		client = Client(cluster)
		processed = []

	#Positive:

		path = yaml['metabolomicsRoot'] + "/outputs/autocredentialPositive/" + runNumber
		files = os.listdir(path)
		files = [path + '/' + f for f in files]

		outputPath = yaml['metabolomicsRoot'] + '/outputs/centroidFxnPositive/' + runNumber
	

		print("here is outputPath for centroidFxn:")
		print(outputPath)

		print("here is first file input for centroid:")
		print(files[0])

	
		if not os.path.exists(outputPath):
			os.makedirs(outputPath, exist_ok=True)


		print("here is script path:")
		print(root + "/executables/centroidFxn.R")

		for f in files:
			print(f)
			print(outputPath)
			processed.append(client.submit(runScript, root + "/executables/centroidFxn.R", f, outputPath))

		progress(processed)
		client.gather(processed)
		client.shutdown()

	#Negative:
	
	if toRun[6] == 1:

		cluster = HTCondorCluster(cores=1, n_workers = 1, memory = "4 GB", disk = "4 GB", local_directory = root)
		cluster.scale(jobs = len(files))
		client = Client(cluster)
		processed = []



		path = yaml['metabolomicsRoot'] + "/outputs/autocredentialNegative/" + runNumber
		files = os.listdir(path)
		files = [path + '/' + f for f in files]
		outputPath = yaml['metabolomicsRoot'] + '/outputs/centroidFxnNegative/' + runNumber

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

		root = yaml['metabolomicsRoot']
#FOR BASE CODE:
#	files = os.listdir(root + "/outputs/hdf5FilesNewColumn/" + runNumber)

	#For DOE
		files1 = glob.glob('/home/ahubbard/text_files_for_all_DOE_metabolomics_samples/HILIC_plate_*/*/*hdf5newColumn_plate.hdf5')
		files2 = glob.glob('/home/ahubbard/HILIC_plate*/*/*hdf5newColumn_plate.hdf5')
	
		files = files1+files2

#	files = [root + "/outputs/hdf5FilesNewColumn/"+ runNumber + '/' + f for f in files]

		print(files[0])

	#select which cutoff value to use from autocredential
		centroidMasses = yaml['metabolomicsRoot']
		centroidMasses = centroidMasses + "/outputs/centroidFxnPositive/" + runNumber
		massFiles = os.listdir(centroidMasses)
		massFiles = [centroidMasses + "/" + s for s in massFiles]	
		maxMass = max(massFiles, key = lambda x: os.stat(x).st_size)
	
		print(maxMass)
		print("is the centroidMasses file used")



		outputPath = yaml['metabolomicsRoot'] + '/outputs/tempFilesPositive/' + runNumber

		if not os.path.exists(outputPath):
			os.makedirs(outputPath, exist_ok=True)



		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "4 GB", local_directory = root)
		cluster.scale(jobs = len(files))
		client = Client(cluster)
		processed = []
		for f in files:
			processed.append(client.submit(runPythonScript, root + "/executables/create_all_the_temp_files.py", f, maxMass, outputPath))
	
		progress(processed)
		client.gather(processed)
		client.shutdown()
	

	#Negative:
	

	if toRun[8] == 1:
		root = yaml['metabolomicsRoot']
		files = os.listdir(root + "/outputs/hdf5FilesNewColumn/" + runNumber)
		files = [root + "/outputs/hdf5FilesNewColumn/"+ runNumber + '/' + f for f in files]


        #select which cutoff value to use from autocredential
		centroidMasses = yaml['metabolomicsRoot']
		centroidMasses = centroidMasses + "/outputs/centroidFxnNegative/" + runNumber
		massFiles = os.listdir(centroidMasses)
		massFiles = [centroidMasses + "/" + s for s in massFiles]
		maxMass = max(massFiles, key = lambda x: os.stat(x).st_size)

		print(maxMass)
		print("is the centroidMasses file used")



		outputPath = yaml['metabolomicsRoot'] + '/outputs/tempFilesNegative/' + runNumber

		if not os.path.exists(outputPath):
			os.makedirs(outputPath, exist_ok=True)



		cluster = HTCondorCluster(cores=1, n_workers=1, memory="10 GB", disk = "4 GB", local_directory = root)
		cluster.scale(jobs = len(files))
		client = Client(cluster)
		processed = []
		for f in files:
			processed.append(client.submit(runPythonScript, root + "/executables/create_all_the_temp_files.py", f, maxMass, outputPath))

		progress(processed)
		client.gather(processed)
		client.shutdown()


	#Step7: Merge files for each mass and generate scan/mass v intensity plots for each


	#Positive
	if not os.path.exists(yaml['metabolomicsRoot']+"/outputs/outputPdfsPositive/"+ runNumber):
		os.makedirs(yaml['metabolomicsRoot']+"/outputs/outputPdfsPositive/"+runNumber, exist_ok=True)

	plotPath = yaml['metabolomicsRoot']+"/outputs/outputPdfsPositive/"+ runNumber

	if not os.path.exists(yaml['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassPositive/"+runNumber):
		os.makedirs(yaml['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassPositive/"+runNumber, exist_ok=True)

	vecOfIntensitiesPath = yaml['metabolomicsRoot']+"/outputs/vecOfIntensitiesForEachMassPositive/"+runNumber

	massDirectories = os.listdir(root + '/outputs/tempFilesPositive/' + runNumber)
	massDirectories = [root + '/outputs/tempFilesPositive/' + runNumber + '/' + m for m in massDirectories]


#	massDirectories = os.listdir('/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/DOE_setaria_and_HILIC/')
#	massDirectories = ['/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/DOE_setaria_and_HILIC/' + m for m in massDirectories]
	


	cluster = HTCondorCluster(n_workers=1, cores=1, memory="10 GB", disk = "4 GB", local_directory = root)
	cluster.scale(jobs = len(massDirectories))
	client = Client(cluster)
	processed = []

	if toRun[9] == 1:
		for m in massDirectories:
			processed.append(client.submit(runScript, root + "/executables/read_in_the_masses_tempfiles.R", m, yaml['fileList'], plotPath, vecOfIntensitiesPath))

	progress(processed)
	client.gather(processed)
	client.shutdown()	

