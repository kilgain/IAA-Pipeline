library(data.table)
library(dplyr)
library(stringr)


callAllOtherFxns = function(segment1)
{
  
  table = segment1
  print(length(table$scan))
  print("^^^ is length of data loaded into callALlOtherFxns")
  #read in the full temp file into a single table
  megaTableII = as.data.table(table)
  print(length(megaTableII$scan))
  print(" ^^^ is length of data loaded into data.table")
  megaTableII$File <- str_remove(megaTableII$File, ".rawoutput.txt")
  print("made megatable")
  megaTableII$Mass = NULL
  megaTableII$polarity = NULL
  megaTableII$Mstype = NULL
  megaTableII$file = NULL
  print(head(megaTableII))

  colnames(megaTableII) =  c("file", "intensity","scan", "compound")
  
  print(head(megaTableII))
  megaTableII$compound = as.numeric(megaTableII$compound)
  megaTableII$compound[is.na(megaTableII$compound)] = 0 
  
  myMass = mean(megaTableII$compound)
  myRT = median(megaTableII$scan)
  print("adding cols")

  myCompound = myMass
  
  
  #print(megaTableII)
  
  
  #note that this may cause a crash if we don't have a situation where
  #every file has the mass ...
  #get the files present in the file
  myFilesInSubset = unique(megaTableII$file)
  
  print(i)
  print("is i")
  #myMass = vecOfMasses[i]
  myMass = as.numeric(myMass)
  
  #starts at the more tolerant ppm threshold, note well
  ppmAllowed = 5
  rangeAllowed = ppmAllowed * .000001 * myMass
  #rangeAllowed = .0001
  #get the information just for that one mass
  #but, in the FULL set of files
  #actually, we've already done the subsetting
  mySubsetFull = megaTableII
 
  #randomly choose one of the samples to align against as the central one we'll align against
  #actually, don't do it randomly - choose the sample with the most entries

  print("going into callPeaks")
  
  #get the start and the stop info for the regions that we want to sum up
  #using the complete chromatogram

  #first run everything through at the file level
  
  #table = read.table(file = "/Users/ahubbard/Downloads/justin_test_data/temp_files_300_samples_setaria_hilic/temp_files_more_stringent_filtering/365.107659166667.txt", header = TRUE, sep = ",")
  #table1SubPlot = twoPeaks[twoPeaks$File == "HILIC_Set_94_1_1_H1_85_TB_setaria_12_0511.rawoutput.txt",]
  
  realignedFilesII = megaTableII
  
  #realignedFilesII = generalizeWarping(megaTableII)
  #realignedFilesII = realignedFilesII[realignedFilesII$scan > 0,]
  #print("made it past warping")
  
  #focus only on the scans that are outliers in BOTH mass and intensity!
  
  
  #peakSearch = realignedFilesII
  #actually, don't rely on the aligned files now!
  
  #print( peakSearch)
  #print("is  peakSearch")
  
  
  
  myWindowStartVec = c(min(megaTableII$scan))
  myWindowStopVec = c(max(megaTableII$scan))
  
  #myWindowStartVec = myStartWindows
  #myWindowStopVec = myEndWindows
  
  #have a table that's going to include all of the summed signals for each
  #mass across each file
  
  #keep track of all of the information for that given mass-feature across all of the files
  filesPresent = c(1:length(myFilesInSubset))
  
  
  #Now we finish the RT window pipeline by summing up the signals
  #across all scans for a given mass from the aligned chromatograms
  #focusing on all rhe regions id's by shatterFind
  
  print(myWindowStartVec)
  for(coordinates in 1:length(myWindowStartVec))
  {
    myWindowStart =  myWindowStartVec[coordinates]
    myWindowStop = myWindowStopVec[coordinates]
    
    
    #have a table to store all of the information
    allInfoTableFile = matrix(nrow= 0, ncol = length(c("compound","rt", fileList)))
    allInfoTableFile = as.data.frame(allInfoTableFile)
    
    #note that this could be challenging if we don't have 
    #the files that we need in the subset
    colnames(allInfoTableFile) =c("compound","rt", as.character(fileList))
    
    
    #include the mass info
    allInfoTableFile[1,colnames(allInfoTableFile) == "compound"] = myCompound
    
    #make sure to inclue the RT info as well
    #myRT = mean(myWindowStart, myWindowStop)
    
    
    allInfoTableFile[1,colnames(allInfoTableFile) == "rt"] = myRT
    
    
    # myCompound = unique(mySubset$compound)[1]
    myCompound = myMass
    for(k in 1:length(myFilesInSubset))
    {
      intensity = 0
      myAdjustedTable = realignedFilesII[realignedFilesII$file == myFilesInSubset[k],]
      print(nrow(myAdjustedTable ))
      print("is the nrow of my table subsetted by the file")
      
      if(nrow(myAdjustedTable) > 5)
      {
        
        toSmoothAndSum = myAdjustedTable
        toSmoothAndSum$file = NULL
        print(head(toSmoothAndSum))
        colnames(toSmoothAndSum) = c("intensity", "scan", "compound")
        #smoothedAddUp = toSmoothAndSum$intensity
        
        
        #myAdjustedTable = as.data.table(myAdjustedTable)
        subsetToAdd = toSmoothAndSum[toSmoothAndSum$scan > myWindowStart & toSmoothAndSum$scan < myWindowStop, ]$intensity
        subsetToAdd = subsetToAdd[subsetToAdd != 0]
        smoothedAddUp = subsetToAdd
        
        #quantileTest = quantile(subsetToAdd, c(.25, .85), na.rm = TRUE)
        
        
        if(length(subsetToAdd) > 0)
        {
          
          #print(quantileTest)
          #print(" is quantileTest")
          #print("greater than 10 in length")
          
          #Old Quantile Filter = [smoothedAddUp > quantileTest[1] & smoothedAddUp < quantileTest[2]]
          
          intensity = sum(smoothedAddUp)
        }
        
        #if we have almost all zeroes, taking the quantile is not going to be useful
        if(length(subsetToAdd) == 0)
        {
          

          intensity = 0
        }
      }
      
      #if(nrow(myAdjustedTable) < 25)
      #{
      #  intensity = 0
      #}
      #intensity = mean(subsetToAdd$intensity)
      print(intensity)
      print("is the intensity")
      allInfoTableFile[1,colnames(allInfoTableFile) == myFilesInSubset[k]] = intensity
      print(myFilesInSubset[k])
      print(allInfoTableFile[1, colnames(allInfoTableFile) == myFilesInSubset[k]])
    }
    
    #if(sigRatioFirst > 3)
    #{
    
    allInfoTableFile[is.na(allInfoTableFile)] = 0
    #allInfoTableAllFiles = rbind(allInfoTableAllFiles, allInfoTableFile)
    #}	
  }
  
  
  #make sure to turn all NA's to 0
  allInfoTableFile[is.na(allInfoTableFile)] = 0
  
  return(allInfoTableFile)
}







args = commandArgs(trailingOnly=TRUE)
#vecOfMasses = basename(args[1])
stringOfPaths = args[1]

#arg1 = arrayOfPathsToTempFiles
#/home/lconnelly/Metabolomics/outputs/tempFilesPositive/mass1/, /home/lconnelly/Metabolomics/outputs/tempFilesPositive/mass2/, ...
#basename(arg1) = arrayOfMasses 
#mass1, mass2...

#In order of: PlantH+, PlantNa+, StandardH+, StandardNa+, PlantH+Iso, PlantNa+Iso, StandardH+Iso, StandardNa+Iso
#vecOfPaths = c("C:/Users/Louis/OneDrive/Desktop/TestThreeTempFilesNewBoundaryFinder/5CPlant/148.0609",
#               "C:/Users/Louis/OneDrive/Desktop/TestThreeTempFilesNewBoundaryFinder/5CPlantNa/170.04284",
#               "C:/Users/Louis/OneDrive/Desktop/TestThreeTempFilesNewBoundaryFinder/5CYeast/153.0779",
#               "C:/Users/Louis/OneDrive/Desktop/TestThreeTempFilesNewBoundaryFinder/5CYeastNa/175.05984",
#               "C:/Users/Louis/OneDrive/Desktop/TestThreeTempFilesNewBoundaryFinder/5CPlantIso/149.06426",
#               "C:/Users/Louis/OneDrive/Desktop/TestThreeTempFilesNewBoundaryFinder/5CPlantNaIso/171.0462",
#               "C:/Users/Louis/OneDrive/Desktop/TestThreeTempFilesNewBoundaryFinder/5CYeastIso/154.08126",
#               "C:/Users/Louis/OneDrive/Desktop/TestThreeTempFilesNewBoundaryFinder/5CYeastNaIso/176.0632")

vecOfPaths = unlist(strsplit(stringOfPaths, ","))
vecOfMasses = basename(vecOfPaths)
#mergedfileEachMass = list()
#for each mass in vecOfMAsses we will have main species 
#and adducts
#{

#if(i == 1){(filePath = plant5C + vecOfMasses 
#if(i == 2)((filePath = plant5CIso + vecOfMasses
#for each filePath, is merghed file here? if not - make it

#outputPdfPath = args[2]
#pdf(outputPdfPath)
mergedFilesEachMass = vector()
for(i in 1:length(vecOfPaths)) {
  massSubDirectory = vecOfPaths[i]
  file_list <- list.files(path=massSubDirectory, full.names = TRUE)
  file_list = file_list[sapply(file_list, file.size) > 0]
  #print(head(file_list))
  fullFileName = paste(massSubDirectory, "/",  vecOfMasses[i],"merged.txt", sep = "")
  mergedFilesEachMass = c(mergedFilesEachMass, fullFileName)
  #print("merged file name is:")
  #print(fullFileName)
  if(!(fullFileName %in% file_list))
  {
    print("before read in")
    allFiles = read.table(file_list[1], header = TRUE, sep = ",")
    print("after read in")
    print(head(allFiles))
    print(" is allFiles before joining")
    file_list = file_list[-1]
    if(length(file_list) > 0)
    {
      for (j in 1:length(file_list))
      {
        print(j)
        print("is j")
        if(file_list[j] != fullFileName)
        {
          oneFile = read.table(file_list[j], header = TRUE, sep = ",")
          print(head(oneFile))
          allFiles = rbind(allFiles, oneFile)
        }
      }
    }
    print(paste("mergedTempFile nrow =  ", length(allFiles$newMass)))
    print("we're writing merged tempfile output now..")
    #Creating merged temp file

    write.table(allFiles, file=fullFileName, row.names = FALSE, sep = "\t")

  }else {
#    allFiles = read.table(fullFileName, sep = "\t")
    print(paste("found merged file for", vecOfMasses[i]))
  }
  
}



#once we have a mergedFile add the table to the list of mergedFiles
print(mergedFilesEachMass)
#mergedfilesEachMass[[1]] = mergedFileMainIon
#mergedfileEachMass[[2]] = mergedFileIsotopologue



#}

#For manual Examination:
# Yeast5C = read.table(mergedFilesEachMass[3], header = TRUE, sep = "\t")
# Yeast5CNa = read.table(mergedFilesEachMass[4], header = TRUE, sep = "\t")
# Yeast5CIso = read.table(mergedFilesEachMass[7], header = TRUE, sep = "\t")
# Yeast5CNaIso = read.table(mergedFilesEachMass[8], header = TRUE, sep = "\t")
# 
# 
# Plant5C = read.table(mergedFilesEachMass[1], header = TRUE, sep = "\t")
# Plant5CNa = read.table(mergedFilesEachMass[2], header = TRUE, sep = "\t")
# Plant5CIso = read.table(mergedFilesEachMass[5], header = TRUE, sep = "\t")
# Plant5CNaIso = read.table(mergedFilesEachMass[6], header = TRUE, sep = "\t")


#Once we;ve create the list of mergedfesEachMAss, we're going to loop through
#and look for the scans that are common to all after polishing


fileOrder = c("Plant", "PlantNa", "Yeast", "YeastNa", "PlantIso", "PlantNaIso")
plotPath = args[3]
plotPath = paste(plotPath, "/", vecOfMasses[1], ".png", sep = "")
print(plotPath)
png(plotPath, width = 1440, height = 7200)
par(mfrow = c(8,2))
for(i in 1:length(vecOfPaths)) {
  if(i == 1) {
    #reference is first file in inputs- with multiple adducts this will be plant endogenous
    reference = read.table(mergedFilesEachMass[i], header = TRUE, sep = "\t")
    reference = reference[reference$polarity == "Positive",]
    print("made it past reference read-in")    
    
    #keep only those points in the top 50% of intensities
    referencePolished = reference %>% filter(Intensity > quantile(reference$Intensity, p = .5))
    #referenceFirstPolish is the dataset that will be used for boundaryExpansion
    referenceFirstPolish = referencePolished
    plot(referencePolished$newMass, referencePolished$Intensity, main = "reference polished top 50% intensity")
    plot(referencePolished$scan, referencePolished$Intensity, main = "reference polished top 50% intensity")
    #calculates numpoints at each scan and the scans that are in the top 20% of numpoints
    scanCount = aggregate(referencePolished$Intensity, by=list(Category=referencePolished$scan), FUN=length)
    importantScans = scanCount[scanCount$x > 30,]
    #importantScans = scanCount[scanCount$x > quantile(scanCount$x, p = .8),]
    #keep only those scans in the top 20% of numpoints(Scan)
    referencePolished = referencePolished[referencePolished$scan %in% importantScans$Category,]
    maxSignal = importantScans[which.max(importantScans$x),]$Category
    bounds = c((maxSignal-10), (maxSignal+10))
    #this is old set of datapoints used in expansion
    #referenceFirstPolish = referencePolished
    plot(referencePolished$newMass, referencePolished$Intensity, main = "reference polished intensity + top20% scan#")
    plot(referencePolished$scan, referencePolished$Intensity, main = "reference polished intensity + top20% scan#", xlim = c(0,3000))
    print("made it past reference polishing")
  }
}
# else {
#    #read in adduct- filter top 50% intensity and top 20% scan buckets
#    file2 = read.table(mergedFilesEachMass[i], header = TRUE, sep = "\t")
#    file2 = file2[file2$polarity == "Positive",]
#    file2 = file2 %>% filter(Intensity > quantile(file2$Intensity, p = .5))
#    scanCount = aggregate(file2$Intensity, by=list(Category=file2$scan), FUN=length)
#    importantScans = scanCount[scanCount$x > quantile(scanCount$x, p = .8),]
#    file2 = file2[file2$scan %in% importantScans$Category,]
#    plot(file2$scan, file2$Intensity, main = paste("file", fileOrder[i], "polished intensity + top20% scan#"), xlim = c(0,3000))
#    referencePolishedBeforeSubset = referencePolished
#    referencePolished = referencePolished[referencePolished$scan %in% importantScans$Category,]
#    #if using this adduct would result in the number of points in the reference going to 0, do not use this adduct for polishing
#    #NOTE!!!! This is not a good approach! I just did this to prevent the algorithm from throwing an error!
#    if(length(referencePolished$Intensity) == 0) {
#	referencePolished = referencePolishedBeforeSubset
#    }
#    print(length(referencePolished$Intensity))
#    plot(referencePolished$scan, referencePolished$Intensity, main = paste("reference polished intensity + top20% scan# after file",fileOrder[i], sep = " "), xlim = c(0,3000))
#    print(paste("we made it past file...", i))
#  }
#}

#scans that are high frequency/intensity across all adducts/isotopologues minus yeast isotopologues
#bounds = range(referencePolished$scan)
oldStart = bounds[1]
oldStop = bounds[2]

plot(referencePolished$scan, referencePolished$Intensity, main = "initial bounds identified in referencePolished")
abline(v = bounds[1], col = 'red')
abline(v = bounds[2], col = 'red')

#BoundaryFinder



start = bounds[1]
stop = bounds[2]
print(paste(start, "is initial start"))
print(paste(stop, "is initial stop"))
expanding = TRUE
while(expanding) {
  
  maxScan = stop
  scanAhead = maxScan + 1
  x = 0
  scans = vector()
  correlationsFilesGreaterThan10 = vector()

  # look ahead in a window of 5 scans
  while(x < 10){
    
    #get the data for the reference and query
    sub1 = referenceFirstPolish[referenceFirstPolish$scan == maxScan,]
    sub1 <- sub1 %>% 
      group_by(File) %>% 
      filter(Intensity == suppressWarnings(max(Intensity))) %>% 
      distinct
    sub2 = referenceFirstPolish[referenceFirstPolish$scan == scanAhead,]
    sub2 <- sub2 %>% 
      group_by(File) %>% 
      filter(Intensity ==suppressWarnings(max(Intensity))) %>% 
      distinct
    
    #keep only the datapoints from files that show up in both our reference and query
    forCor1 = sub1[sub1$File %in% intersect(sub1$File, sub2$File),]
    forCor2 = sub2[sub2$File %in% intersect(sub1$File, sub2$File),]
    
    forCor1 = forCor1[order(forCor1$File),]
    forCor2 = forCor2[order(forCor2$File),]
    
    #determine the correlation
    correlation = cor(forCor1$Intensity, forCor2$Intensity)
    
    #if we have more than 10 files, giving us good statistical power
    # record the correlation and the scan 
    if(length(forCor1$Intensity) > 10) {
      #print(paste(correlation, "is correlation of scan", maxScanPlant, "with scan", scanAhead))
      correlationsFilesGreaterThan10 = c(correlationsFilesGreaterThan10, correlation)
      scans = c(scans, scanAhead)
    }
    scanAhead = scanAhead + 1
    x = x + 1
  }
  #print(max(correlationsFilesGreaterThan10))
  #plot(x = scans, y =correlationsFilesGreaterThan10)

  #if, in our window of 5 scans, one of them had a cor > .75 w/ our reference
  #(and also 10 files or more were shared, expand outwards)
  if(max(correlationsFilesGreaterThan10) > .75) {
    print(max(correlationsFilesGreaterThan10))
    index = which.max(correlationsFilesGreaterThan10)
    stopOld = stop
    stop = scans[index] 
    shift = stop - stopOld
    start = start - shift
    print(paste(start, "is updated start"))
    print(paste(stop, "is updated stop"))
  }else{
    expanding = FALSE
    print(paste(start, "is updated start"))
    print(paste(stop, "is updated stop"))
    print("expansion has terminated")
  }
}

plot(reference$scan, reference$Intensity, main = "old bounds (red), new bounds (green)")
abline(v = oldStart, col = 'red')
abline(v = oldStop, col = 'red')
abline(v = start, col = 'green')
abline(v = stop, col = 'green')

#We are going to take c("Plant", "PlantNa", "Yeast", "YeastNa", "PlantIso", "PlantNaIso"), subset them by the bounds identified (start/stop), and then run CallAllOtherFunctions on them to generate a VecOfIntensities

#get all file names
fileList = read.table(args[2], header = T)$x
fileList = str_remove(fileList, ".rawoutput.txt")
fileList = unique(fileList)

allInfoTableAllAdducts = matrix(nrow= 0, ncol = length(c("compound","rt", fileList)))
allInfoTableAllAdducts = as.data.frame(allInfoTableAllAdducts)
colnames(allInfoTableAllAdducts) =c("compound","rt", as.character(fileList))


for(i in 1:length(vecOfPaths)) {
 adduct = read.table(mergedFilesEachMass[i], header = TRUE, sep = "\t")
 adduct = adduct %>% filter(scan > start & scan < stop & polarity == "Positive")
 #run anovalign here
 #take anovaligned points and call boundary finder again here
 #subset by new boundary finder bounds and run callAllOtherFxns
 tableForAdduct = callAllOtherFxns(adduct)
 allInfoTableAllAdducts = rbind(allInfoTableAllAdducts,tableForAdduct)
}
#rownames(allInfoTableAllAdducts) = fileOrder
write.csv(allInfoTableAllAdducts, file =paste(args[4], "/", vecOfMasses[1], ".txt", sep = ""), row.names=F)

dev.off()
