import vaex
import pandas as pd
import sys
import os
#remove yaml, getopt, functools, pathlib, numpy



fileRead = sys.argv[1]

testMasses = pd.read_csv(sys.argv[2])


#for citrulline test
#myMass = 176.1034
#myIsotope = 176.1034 + 1.00336
#testMasses = [myMass, myIsotope]


testMasses = testMasses['x'].values
#for running subset of masses/isotopologues of masses by hand 
#testMasses = testMasses[0:10]
#testMasses = [176.1031,177.1001,177.1073, 182.12326, 183.1203,183.1275]

fileRead = str(fileRead)
dv = vaex.open(fileRead)

#loop through the array of scans
for i in testMasses:
	print(i)
	i = float(i)
	i = round(i, ndigits = 5)
	mass = i
  #can use this line to make it an isotopologue!
	#mass = mass + 1.00336
	ppm = 10
	myMassRange = .000001 * ppm * mass
	subsetChromatogramI = dv[dv.newMass > (mass - myMassRange)]
	subsetChromatogramII = subsetChromatogramI[subsetChromatogramI.newMass < (mass + myMassRange)]

	tempPath = sys.argv[3]

	if not os.path.exists(tempPath+'/'+str(i)):
		os.makedirs(tempPath+'/'+str(i), exist_ok=True)

  
	base=os.path.basename(fileRead)
	fullName = tempPath + "/" + str(i) + '/' + str(i) + base+ ".txt"
	print(fullName)
	print("printing out now")
	print(subsetChromatogramII)
	print("is what we are printing")
	print(len(subsetChromatogramII))
	if len(subsetChromatogramII) > 0:
		subsetChromatogramII.export_csv(path=fullName)
