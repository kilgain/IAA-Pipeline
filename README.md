# IAA-Pipeline
Scripts to run Isolock/Autocredential/Anovalign

# NOTE! This is a work in progress, further instructions and details will be added soon!

- For questions, please reach out to AHubbard@danforthcenter.org and LConnelly@danforthcenter.org

- Make sure all output directories exist before running the scripts!

# Step 1: Activate conda environment

- Create the environment from the environment.yml file:
     `conda env create -f environment.yml`

- The first line of the yml file sets the new environment's name. 

- Activate the new environment: conda activate IAA

- Verify that the new environment was installed correctly:
    `conda env list`



# Step 2: Convert .raw Files To .txt Containing Relevant Information.

  How to run on command line using fileProcessor.exe:
  
  `mono /path/to/fileProcessor.exe /path/to/input.raw /path/to/outputDirectory`
  
  Example output: soon™

# Step 3: Add File Column To .txt Output:

  How to run on command line using fileAdd.R:
  
  `Rscript /path/to/fileAdd.R /path/to/rawoutput.txt`
  
  Example output: soon™

# Step 4: Convert .txt Files to HDF5 Format:

  How to run on command line using create_the_hdf5_files.py:
  
  `python /path/to/create_the_hdf_files.py /path/to/input.txt /path/to/outputDirectory`
  
  Example output: soon™

# Step 5: Run Isolock on Each File:

  How to run on command line using isolock.py:
  
  `python /path/to/isolock.py /path/to/input.hdf5 /path/to/HDF5FileToUseForReference.hdf5 /path/to/outputDirectory`
  
  Example output: soon™

# Step 6: Run Autocredential

  Note: Autocredential filters the mass spectra by positive or negative polarity. To run positive mode, set polarity = "positive". To run negative mode, set polarity =   "negative". To get output for both, set polarity = "both". Make sure the output paths exist before running. Provide a dummy path for the polarity mode you do not wish to run (this will be addressed in a future release).

  How to run on command line using autocredentialVaex.py:
  
  `python /path/to/autocredentialVaex.py "polarity" /path/to/directoryOfInputFiles /path/for/positiveModeOutputDirectory /path/for/negativeModeOutputDirectory`
  
  Example output: soon™
  
# Step 7: Centroid Masses From Autocredential

  Autocredential makes a set of outputs based on how stringently the mass spectra was filtered. For most applications, the filter value of 20 will be adequate. For large-scale applications (dozens to hundreds of files), this filter value should be higher (30-50). Future releases will auto-determine filter to use.

  How to run on command line using centroidFxn.R:
  
  `Rscript /path/to/centroidFxn.R /path/to/autocredentialOutput.txt /path/to/outputDirectory`
  
  Example output: soon™
  
# Step 8: Create TempFiles For Each Mass From CentroidFxn (all datapoints within +-10 ppm of mass)

  How to run on command line using create_all_the_temp_files.py:
  
  `python /path/to/create_all_the_temp_files.py /path/to/isolockOutputFile.hdf5 /path/to/centroidFxnOutput.txt /path/to/outputDirectory`
  
  Example output: soon™
  
# Step 9: Find Mass Features:

  Note: This step requires a text file which contains the names of all your files in the following format:
  "x"
  file1.rawoutput.txt
  file2.rawoutput.txt
  ...

  How to run on command line using newestAnovalign.R:
  
  `Rscript /path/to/newestAnovalign.R /path/to/tempFileDirectory /path/to/textFileOfAllFileNames.txt /path/to/outputDirectory`
  
  Example output: soon™
  
  

