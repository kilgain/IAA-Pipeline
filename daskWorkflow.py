from dask_jobqueue.htcondor import HTCondorCluster
from dask.distributed import Client, progress
import os
import yaml


def runScript(x,y):
	os.system("Rscript " + x + " " + y)


def runPythonScript(x,y, z = ""):
	os.system("python3 " + x + " " + y + " " + z)


if __name__ == "__main__":
	#Step1: Add File Column and Collapse Duplicate Masses
	cluster = HTCondorCluster(cores=10, memory = "4 GB", disk = "4 GB", local_directory = "/home/lconnelly/Metabolomics")
	with open("/home/lconnelly/Metabolomics/Config.yaml") as file:
		config = yaml.safe_load(file)
	config = config['samplesDirectory']
	files = os.listdir(config)
	cluster.scale(jobs = len(files))
	client = Client(cluster)
	#list of job futures
	processed = []
	#for loop of submissions
	f = open(config + '/' + files[0], 'r')
	line1 = f.readline()
	f.close()		


	if not '"File"' in line1:
		for f in files:
			processed.append(client.submit(runScript, "/home/lconnelly/Metabolomics/step1Conversion/fileAdd.R", f))
	progress(processed)
	client.gather(processed)
	client.shutdown()


	#Step2: Convert .txt Files to HDF5 Format
	cluster = HTCondorCluster(cores=10, memory = "4 GB", disk = "4 GB", local_directory = "/home/lconnelly/Metabolomics")
	cluster.scale(jobs = len(files))
	client = Client(cluster)
	processed = []

	config = config.replace("rawDataFiles", "hdf5Files")
	x = os.path.isdir(config)
	if x == False:
		os.mkdir(config)
		for f in files:
			processed.append(client.submit(runPythonScript, "/home/lconnelly/Metabolomics/step1Conversion/create_the_hdf5_files.py", f))
	progress(processed)
	client.gather(processed)
	client.shutdown()



	#Step3: Run Isolock On Each File
	cluster = HTCondorCluster(cores=10, memory = "4 GB", disk = "4 GB", local_directory = "/home/lconnelly/Metabolomics")
	cluster.scale(jobs = len(files))
	client = Client(cluster)
	processed = []
	x = os.listdir(config)
	x = [i for i in x if ".hdf5newColumn" in i]
	x = x[0]

	
	if not os.path.isfile(config + '/' + x):
		for f in files:
			processed.append(client.submit(runPythonScript, "/home/lconnelly/Metabolomics/step1Conversion/step2Isolock/python_script_determine_offset.py", f))
	progress(processed)
	client.gather(processed)
	client.shutdown()

	#Step4: Run Autocredential
	cluster = HTCondorCluster(cores=10, memory = "4 GB", disk = "4 GB", local_directory = "/home/lconnelly/Metabolomics")
	cluster.scale(jobs = 1)
	client = Client(cluster)
	processed = []

	with open("/home/lconnelly/Metabolomics/Config.yaml") as file:
        	config = yaml.safe_load(file)

	config = config['metabolomicsRoot']
	config = config + "/step1Conversion/step2Isolock/step3Autocredential/autocredentialCutoffs"
	x = os.path.isdir(config)
	if x == False:
		os.mkdir(config)
		processed.append(client.submit(runPythonScript, "/home/lconnelly/Metabolomics/step1Conversion/step2Isolock/step3Autocredential/autocredential.py"))
	progress(processed)
	client.gather(processed)
	client.shutdown()

	#Step5: Centroid the Masses
	files = os.listdir(config)
	cluster = HTCondorCluster(cores=10, memory = "4 GB", disk = "4 GB", local_directory = "/home/lconnelly/Metabolomics")
	cluster.scale(jobs = len(files))
	client = Client(cluster)
	processed = []
	with open("/home/lconnelly/Metabolomics/Config.yaml") as file:
		config = yaml.safe_load(file)
	config = config['metabolomicsRoot']
	config = config + "/step1Conversion/step2Isolock/step3Autocredential/centroidMasses"
	x = os.path.isdir(config)
	if x == False:
		os.mkdir(config)
		for f in files:
			processed.append(client.submit(runScript, "/home/lconnelly/Metabolomics/step1Conversion/step2Isolock/step3Autocredential/centroidFxn.R", f))
	progress(processed)
	client.gather(processed)
	client.shutdown()

	#Step6: Subset the Masses
	with open("/home/lconnelly/Metabolomics/Config.yaml") as file:
		config = yaml.safe_load(file)

	root = config['metabolomicsRoot']
	files = os.listdir(root + "/hdf5Files")
	subset = "newColumn"
	files = [i for i in files if subset in i]
	files = [root + "/hdf5Files/" + f for f in files]


	#select which cutoff value to use from autocredential
	centroidMasses = config['metabolomicsRoot']
	centroidMasses = centroidMasses + "/step1Conversion/step2Isolock/step3Autocredential/centroidMasses"
	massFiles = os.listdir(centroidMasses)
	massFiles = [centroidMasses + "/" + s for s in massFiles]	
	maxMass = max(massFiles, key = lambda x: os.stat(x).st_size)
	


	cluster = HTCondorCluster(cores=10, memory="10 GB", disk = "4 GB", local_directory = "home/lconnelly/Metabolomics")
	cluster.scale(jobs = len(files))
	client = Client(cluster)
	processed = []
	if not os.path.isdir(config['tempFiles']):
		for f in files:
			processed.append(client.submit(runPythonScript, "/home/lconnelly/Metabolomics/step1Conversion/step2Isolock/step3Autocredential/step4Subset/create_all_the_temp_files.py", f, maxMass))
	progress(processed)
	client.gather(processed)
	client.shutdown()
	

	#Step7: Merge files for each mass and generate scan/mass v intensity plots for each
	if not os.path.exists(config['metabolomicsRoot']+"/outputPdfs/"):
		os.makedirs(config['metabolomicsRoot']+"/outputPdfs/", exist_ok=True)

	if not os.path.exists(config['metabolomicsRoot']+"/outputCorrelations/"):
                os.makedirs(config['metabolomicsRoot']+"/outputCorrelations/", exist_ok=True)

	massDirectories = os.listdir(config['tempFiles'])
	massDirectories = [config['tempFiles'] + "/" + m for m in massDirectories]
	cluster = HTCondorCluster(cores=10, memory="10 GB", disk = "4 GB", local_directory = "home/lconnelly/Metabolomics")
	cluster.scale(jobs = len(massDirectories))
	client = Client(cluster)
	processed = []
	#if not os.path.exists(massDirectories[0] + '/' + os.path.basename(massDirectories[0] + 'merged.txt')):
	for m in massDirectories:
		processed.append(client.submit(runScript, "/home/lconnelly/Metabolomics/step1Conversion/step2Isolock/step3Autocredential/step4Subset/step5Merge/read_in_the_masses_tempfiles.R", m))
	progress(processed)
	client.shutdown()	
