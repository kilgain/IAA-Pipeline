

args = commandArgs(trailingOnly=TRUE)
#vecOfMasses = basename(args[1])
vecOfPaths = args[1]

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
        if(file_list[j] != mergedExclude)
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
for(i in 1:6) {
  if(i == 1) {
    reference = read.table(mergedFilesEachMass[i], header = TRUE, sep = "\t")
    reference = reference[reference$polarity == "Positive",]
    
    
    referencePolished = reference %>% filter(Intensity > quantile(reference$Intensity, p = .5))
    plot(referencePolished$newMass, referencePolished$Intensity, main = "reference polished top 50% intensity")
    plot(referencePolished$scan, referencePolished$Intensity, main = "reference polished top 50% intensity")
    scanCount = aggregate(referencePolished$Intensity, by=list(Category=referencePolished$scan), FUN=length)
    importantScans = scanCount[scanCount$x > quantile(scanCount$x, p = .8),]
    
    referencePolished = referencePolished[referencePolished$scan %in% importantScans$Category,]
    referenceFirstPolish = referencePolished
    plot(referencePolished$newMass, referencePolished$Intensity, main = "reference polished intensity + top20% scan#")
    plot(referencePolished$scan, referencePolished$Intensity, main = "reference polished intensity + top20% scan#", xlim = c(0,3000))
  } else {
    file2 = read.table(mergedFilesEachMass[i], header = TRUE, sep = "\t")
    file2 = file2[file2$polarity == "Positive",]
    file2 = file2 %>% filter(Intensity > quantile(file2$Intensity, p = .5))
    scanCount = aggregate(file2$Intensity, by=list(Category=file2$scan), FUN=length)
    importantScans = scanCount[scanCount$x > quantile(scanCount$x, p = .8),]
    file2 = file2[file2$scan %in% importantScans$Category,]
    plot(file2$scan, file2$Intensity, main = "file 2 polished intensity + top20% scan#", xlim = c(0,3000))
    referencePolished = referencePolished[referencePolished$scan %in% importantScans$Category,]
    plot(referencePolished$scan, referencePolished$Intensity, main = "reference polished intensity + top20% scan# after file2", xlim = c(0,3000))
  }
}

#scans that are high frequency/intensity across all adducts/isotopologues minus yeast isotopologues
bounds = range(referencePolished$scan)
oldStart = bounds[1]
oldStop = bounds[2]

plot(reference$scan, reference$Intensity, main = "initial bounds identified")
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
  while(x < 5){
    
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
    
    
    forCor1 = sub1[sub1$File %in% intersect(sub1$File, sub2$File),]
    forCor2 = sub2[sub2$File %in% intersect(sub1$File, sub2$File),]
    
    forCor1 = forCor1[order(forCor1$File),]
    forCor2 = forCor2[order(forCor2$File),]
    
    
    correlation = cor(forCor1$Intensity, forCor2$Intensity)
    
    
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

#dev.off()
