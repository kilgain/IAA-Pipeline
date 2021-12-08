# IAA-Pipeline
Scripts to run Isolock/Autocredential/Anovalign

# NOTE! This is a work in progress, further instructions and details will be added soon!

#Step 1: Activate conda environment



# Step 2: Convert .raw Files To .txt Containing Relevant Information.

  How to run on command line using fileProcessor.exe:
  
  Example output:

# Step 3: Add File Column To .txt Output:

  How to run on command line using fileAdd.R:
  
  Example output:

# Step 4: Convert .txt Files to HDF5 Format:

  How to run on command line using create_the_hdf5_files.py:
  
  Example output:

# Step 5: Run Isolock on Each File:

  How to run on command line using isolock.py:
  
  Example output:

# Step 6: Run Autocredential

  How to run on command line using autocredentialVaex.py:
  
  Example output:
  
# Step 7: Centroid Masses From Autocredential

  Autocredential makes a set of outputs based on how stringently the mass spectra was filtered. For most applications, the filter value of 20 will be adequate. For large-scale applications (dozens to hundreds of files), this filter value should be higher (30-50). Future releases will auto-determine filter to use.

  How to run on command line using centroidFxn.R:
  
  Example output:
  
# Step 8: Create TempFiles For Each Mass From CentroidFxn (all datapoints within +-10 ppm of mass)

  How to run on command line using create_all_the_temp_files.py:
  
  Example output:
  
# Step 9: Find Mass Features:

  How to run on command line using newestAnovalign.R:
  
  Example output:
  
  

