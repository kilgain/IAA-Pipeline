library(data.table)
library(dplyr)
library(yaml)
library(stringr)
#library(gtools)
#library(minpack.lm)
#library(splines)



args = commandArgs(trailingOnly=TRUE)

#config = read_yaml("/home/ahubbard/Metabolomics/Config.yaml")

#args = c("/Users/ahubbard/151.06013", "/Users/ahubbard/files_for_config.txt", "/home/ahubbard/Metabolomics/outputs/outputPdfsPositive/T4_T6_T8_setariaII","/home/ahubbard/Metabolomics/outputs/vecOfIntensitiesForEachMassPositive/T4_T6_T8_setariaII")

vecOfMasses = basename(args[1])
currentDir = args[1]


print(vecOfMasses)
myMassTest = as.numeric(vecOfMasses)

#old adjustment for 27/36 samples v full DOE= - (.00000367)*myMassTest


myMassTest = myMassTest
myMassTest = round(myMassTest, digits = 5)
print(myMassTest)
myMassTestQuery = myMassTest


MassSubDirectory = currentDir
print(MassSubDirectory)
print(" is where we're going to read from")

File_list <- list.files(path=MassSubDirectory, full.names = TRUE)
File_list = File_list[sapply(File_list, file.size) > 0]

fullFileName = paste(currentDir, "/",  vecOfMasses,"merged.txt", sep = "")
mergedExclude = fullFileName

print(fullFileName)

print(File_list)
print(!(fullFileName %in% File_list))


if(!(fullFileName %in% File_list))
{
  
  
  allFiles = read.table(File_list[1], header = TRUE, sep = ",")
  #print(head(allFiles))
  #colnames(allFiles) = c("row","File","Mass","Instensity","Scan", "Mode","MsType")
  
  print(allFiles)
  print(" is allFiles")
  
  
  #print(mergedExclude)
  File_list = File_list[-1]
  if(length(File_list) > 0)
  {
    for (j in 1:length(File_list))
    {
      print(j)
      print("is j")
      if(File_list[j] != mergedExclude)
      {   
        
        #oneFile = tryCatch(read.table(File_list[j], header = TRUE, sep = '|'), error=function(e) NULL)
        oneFile = read.table(File_list[j], header = TRUE, sep = ",")
        print(head(oneFile))
        #colnames(oneFile) = c("File","Mass","Instensity","Scan", "Mode","MsType", "newMass")
        allFiles = rbind(allFiles, oneFile)
        #print(j)
        #print(myMassTest)
      }
    }
  }
  print("we're writing the output now")
  #write to the temp table
  
  write.table(allFiles, file=fullFileName, row.names = FALSE, sep = ",")
  
}else {
  
  allFiles = read.table(fullFileName, header = TRUE, sep = "\t")
  
}	

##TEMPORARY POLARITY FILTER FOR TESTING:

allFiles = allFiles[allFiles$polarity == "Positive",]

#fxn to look for correlation btwn scans
corFxn = function(table, scan1, scan2)
{
  
  table = table
  
  
  scan1 = scan1
  scan2 = scan2
  
  sub1 = table[table$scan == scan1,]
  sub1 <- sub1 %>% 
    group_by(File) %>% 
    filter(Intensity == suppressWarnings(max(Intensity))) %>% 
    distinct
  sub2 = table[table$scan == scan2,]
  sub2 <- sub2 %>% 
    group_by(File) %>% 
    filter(Intensity ==suppressWarnings(max(Intensity))) %>% 
    distinct
  
  #keep only the datapoints from Files that show up in both our reference and query
  forCor1 = sub1[sub1$File %in% intersect(sub1$File, sub2$File),]
  forCor2 = sub2[sub2$File %in% intersect(sub1$File, sub2$File),]
  
  forCor1 = forCor1[order(forCor1$File),]
  forCor2 = forCor2[order(forCor2$File),]
  
  
  
  #if enough Files shared, return the correlation
  if(length(intersect(sub1$File, sub2$File)) > 10) 
  {
    #determine the correlation
    correlation = cor(forCor1$Intensity, forCor2$Intensity)
  }
  
  
  #if enough Files shared, return the correlation
  if(length(intersect(sub1$File, sub2$File)) <= 10) 
  {
    #determine the correlation
    correlation = 0
  }
  
  
  return(correlation)
}

#make boundary finder a function:
expandBounds = function(referenceFirstPolish)
{
  #expandBounds will accept a slice of a chromatogram
  #and use the correlation decay approach to figure out if it contains a region of likely signal.

  #start = bounds[1]
  #stop = bounds[2]
  #print(paste(start, "is initial start"))
  #print(paste(stop, "is initial stop"))
  expanding = TRUE
  referenceFirstPolish = referenceFirstPolish
  #print(head(referenceFirstPolish))
  #take only the high priority scans (those with more than 30 signals)
  
  
  #print(referenceFirstPolish)
  #print("we are going to aggregate now")

   #THIS IS HARDCODED ... but only for now

  #first step is to take the merged chromatogram, and count the numer of signals in each scan
  #this will allow simple denoising by subtracting scans that don't have a lot og signal
  fullTableInfo = aggregate(referenceFirstPolish$Intensity, by=list(Category= referenceFirstPolish$scan), FUN=length)
  fullTableInfo = fullTableInfo[ fullTableInfo$x > 5,]
  
  #maxSignal is going to be the reference scan
  maxSignal = fullTableInfo[which.max(fullTableInfo$x),]$Category
  
  print(paste(maxSignal, "is maxSignal"))

  #this will be the 'meat' of the signal
  maximumSignalBounds = c((maxSignal - 10), (maxSignal + 10))
  
  #look to either side of this main signal to see if we can extend
  mainSlice = fullTableInfo$Category[fullTableInfo$Category >  maximumSignalBounds[1] & fullTableInfo$Category <  maximumSignalBounds[2] ]
  
  #at this point, I have the center of a peak and the the +/- 10 scan windows  
  #I'll also have an array of 'high priority' scans to the right 
  #AND to the left of this region.  The funamentaly question is:  can I extend the bounds of my central signal
  #to the left OR right


  
  #order all of the scans
  majorSignalScans = fullTableInfo$Category[order(fullTableInfo$Category)]
  
  #possible scans to the right
  toIterateOutRight = majorSignalScans[majorSignalScans  > max(mainSlice)]
  
  #possible scans to the left
  toIterateOutLeft = majorSignalScans[majorSignalScans  < min(mainSlice)]
  toIterateOutLeft = rev(toIterateOutLeft)
  
  expanding = TRUE
  
  #move outwards from the main slice
  while(expanding) 
  {
    
    #maxScan = stop
    #scanAhead = maxScan + 1
    x = 1
    scans = vector()
    correlationsFilesGreaterThan10 = vector()
    
    rightWardScans = vector()
    leftWardScans = vector()
    #look to either side of the core of the signal in 5 scan increments
    while(x < 6)
    {
      
      #print(x)
      #print(" is x")
      
      #get the data for the reference and query
      maxScanRight = toIterateOutRight[x] 
      maxScanLeft = toIterateOutLeft[x] 
      
      #call the function to expand out via the correlation now
      rightwardCor = corFxn(referenceFirstPolish, max(mainSlice), maxScanRight)
      rightWardScans = c(rightWardScans , abs(rightwardCor))
      
      leftwardCor = corFxn(referenceFirstPolish, min(mainSlice), maxScanLeft)
      leftWardScans = c(leftWardScans, abs(leftwardCor))
      x = x + 1
    }
    
    
    
    leftWardScans = abs(leftWardScans )
    rightWardScans = abs(rightWardScans )
    myExtension = 0
    
    
    #we've looked five scans on either side.  Now, can we successfully expand or not?
    
    #we're extending rightward
    if(max(rightWardScans) > max(leftWardScans))
    {
      
      #print("better to expand out rightward")
      #the first five elements of the rightward bounds will be our possible extended scans
      myExtendedScan = toIterateOutRight[which.max(rightWardScans)]
      
      #what's the difference between the main slice and the extended scan
      myExtension = abs(max(mainSlice) - myExtendedScan)
      
    }
    
    #we're extending leftward
    if(max(rightWardScans) < max(leftWardScans))
    {
      
      #print("better to expand out leftward")
      #the first five elements of the leftward bounds will be our possible extended scans
      myExtendedScan = toIterateOutLeft[which.max(leftWardScans)]
      
      #what's the difference between the main slice and the extended scan
      myExtension = abs(min(mainSlice) - myExtendedScan)
      
    }
    
    
    #we've discovered that we can extend the bound.  So, let's implement the extension
    #for the right and the left bounds
    #print(paste("here is myExtension:", myExtension))
    #print(paste("here is max(mainSlice)", max(mainSlice)))
    #print(paste("here is min(mainSlice)", min(mainSlice)))
    rightBoundNew = max(mainSlice) + myExtension
    leftBoundNew = min(mainSlice) - myExtension 
    
    
    #update the main slice to include our expanded region now
    #taking scan from both the left and the right depending on which
    #region yielded the best expansion
    mainSlice = c(mainSlice, toIterateOutRight[toIterateOutRight < rightBoundNew])
    mainSlice = c(mainSlice, toIterateOutLeft[toIterateOutLeft > leftBoundNew ])
    
    #now we've added the potential scans 
    mainSlice = mainSlice[order(mainSlice)]
    
    #print(mainSlice)
    #print(" is main slice")
    
    #if we've extended bounds, remove the scans from the potential 
    toIterateOutRight = toIterateOutRight[toIterateOutRight > rightBoundNew]
    toIterateOutLeft = toIterateOutLeft[toIterateOutLeft < leftBoundNew ]
    
    
    #stop expanding if no good candidates
    
    #get the output from BOTH expansion processes
    #determine if we're in a situation that NEITHER are good
    
    
    corsFromBothExpansions = c(rightWardScans, leftWardScans)
    
    if(max(corsFromBothExpansions) < .75 )
    {
      expanding = FALSE
      #print("oops, can't expand now")
    }
    
    #print("expanding")
    
  }
  
  
  

  
  #return the bounds here, now
  bounds = c(min(mainSlice), max(mainSlice))
  
  
  #if our extesnion is less than 10, it's negligble - don't count it
  if(range(mainSlice)[2] - range(mainSlice)[1] < 10) {
    bounds = c(Inf,Inf)
  }
  
  print(paste("bounds:", bounds[1], "-", bounds[2]))
  return(bounds)
}


#Pipeline for Mass features:
vecOfStarts = vector()

boundaryFinder = function(mySubset, mySubsetOG)
{
  mySubset = mySubset

  
  #copy of data to be used for graphing at the end of the boundary-finding. signals are iteratively removed from MySubset
  #In order to avoid recapturing the same signal again and again.
  mySubsetOG = mySubsetOG
  #exit condition for the while loop.
  #turns to false if signal/noise of new peak is low, or if the highest signal in the peak is 15% of the original signal
  lastPeakGood = TRUE
  
  #storage vectors for bounds identified
  vecOfStarts = vector()
  vecOfStops = vector()
  

  #While loops that iteratively isolates peaks until an exit condition is hit (lastPeakGood)
  while (lastPeakGood) {
  
      #print("iterating!")
      #print(head(mySubsetOG))
    
    
      #expandBounds is going to to use the correlation decay approach
      #to find signal regions
      #each time we find a signal region, we're going to subtract it from the 
      #main table, re-run exapnadBounds to see if it can find any additional signal regions
      #and we're going to be keeping track of each of these signal regions (which will be passed to boundaryCollapser)
      boundsPrevious = expandBounds(mySubset)
    #get the area within the bounds
      withInBounds = mySubset[mySubset$scan >= boundsPrevious[1] & mySubset$scan <= boundsPrevious[2],]
    
    #only add a slice if its signal/noise is high enough!
      
      #to claculate the signal for the S/N we are going to add up all of the signal in the region between
      #the bounds 
      signal = sum(withInBounds$Intensity)
      
      
      
      
      #we'll make sure to have removed the signal region from the chromatogram
      signalRemoved = rbind(mySubset[mySubset$scan < boundsPrevious[1],],  mySubset[mySubset$scan > boundsPrevious[2],])
    

      #although we're removing signal regions from the chromatogram, we're all storing the original original chromatogram
      #in order to calculate a representative noise region
      #IF there are 1,000 unqiue scans in my signal region
      #I also want my noise region I am comparing it to to have 1000 scans
      #SO ...if 1,000 scans in the signal region, I will randomly sample 1,000 scans from the orignal
      set.seed(12345)
      noise = mySubsetOG[mySubsetOG$scan %in% sample(unique(mySubsetOG$scan),length(unique(withInBounds$scan))), ]
      
      
      #noise, in the case that out sginal is 1000 scans, will now be the summed signal of 1000 representative 
      #scans cgosen from the original chromatogram
      #because the signal AND noise regions are of the same size, we don't need to normalize by scan  
      noise = sum(noise$Intensity)
      
      #noise = sum(forNoise$Intensity) / (max(forNoise$scan) - min(forNoise$scan))
      
    
    #calculate the signal to noise
      sigToNoise = signal/noise
    
      #only do if we have enough points
      #otherwise, make the bounds ahead initialized as null 
      boundsAhead = vector()
      
      #are we pushing our luck with calling another potential peak?
      #If we've pulled out one peak, can we find another one?
      #because if we can, we're defintely not goign to want to stop the algorithm
      if(nrow(signalRemoved) >  200)
      {
        print("bounds ahead =")
        boundsAhead = expandBounds(signalRemoved)
      }
      
      #if the remaining signal from our peak is less than 200 datapoints
      #we're not going to bother looking for another peak
      if(nrow(signalRemoved) <  200)
      {
        boundsAhead[1] = Inf
        boundsAhead[2] = Inf
      }
      
      
      #rese the chromatogram to the remaining signal once the likely peak is removed
      mySubset = signalRemoved
    
    print(paste(sigToNoise, "is sig to noise"))
    
    #Now, we're going to figure out if we want to stop the algorithm AND keep our bounds
    
    #boundaAhead[[1]] OR boundsAhead[[2]] is Inf we assume that we can't push ahead with more bounds
    #so we'll want to stop looking
     if((is.finite(boundsAhead[1]) == FALSE) || (is.finite(boundsAhead[2]) == FALSE))
    {
      
      #we want to stop looking for bounds
       #but, we we already have a good peak?
      if(sigToNoise > 2)
      {
        #if so, get its bounds 
        vecOfStarts = c(vecOfStarts, boundsPrevious[1])
        vecOfStops = c(vecOfStops, boundsPrevious[2])
      }
       
       #and quit looking for peaks!
      lastPeakGood = FALSE
      
     }
    
    #if  we have more chances to find more peaks, still keep track of the bounds that we have 
    if((is.finite(boundsAhead[1]) == TRUE) && (is.finite(boundsAhead[2]) == TRUE))
    {
      if(sigToNoise > 2)
      {
        vecOfStarts = c(vecOfStarts, boundsPrevious[1])
        vecOfStops = c(vecOfStops, boundsPrevious[2])
      }
    }
    
    }
  return(list(vecOfStarts, vecOfStops))
}


boundCollapser = function(vecOfStarts, vecOfStops) {
  vecOfStarts <- vecOfStarts
  vecOfStops <- vecOfStops
  # print(vecOfStarts)
  # print(vecOfStops)
  if (length(vecOfStarts > 0)) {
    ranges <- cbind(vecOfStarts, vecOfStops + 1)
    row.names(ranges) <- NULL
    colnames(ranges) <- c("start", "stop")
    ranges <- data.table(ranges)
    x <- ranges %>%
      arrange(start) %>%
      group_by(g = cumsum(cummax(lag(
        stop, default = first(stop)
      )) < start - 200)) %>%
      summarise(start = first(start), stop = max(stop))
    
    vecOfStartsFinal <- pull(x[, 2])
    vecOfStopsFinal <- pull(x[, 3])
    return(list(vecOfStartsFinal, vecOfStopsFinal))
  }
}





#anovalign
generalizeWarping <- function(metabolitePlotFull)
{
  metabolitePlotFull = metabolitePlotFull
  #get the max and min of the original zone submitted to the fxn
  #we'll want to make sure that the adjusted zone isn't going to exceed 
  #these limits
  maxOfSubmittedZone = max(metabolitePlotFull$scan)
  minOfSubmittedZone = min(metabolitePlotFull$scan)
  #for the reference for the t-test, use the sample with the 3rd most signal
  #why?  sometimes the highest one can be an outlier/poor chromatography that
  #'smooshes' together signals
  #maxReference = dplyr::count(metabolitePlotFull, File)[order(dplyr::count(metabolitePlotFull, File)$n, decreasing = TRUE), ][3, ]$File
  
  print(head(metabolitePlotFull[order(metabolitePlotFull$Intensity, decreasing = TRUE),]))
  
  maxReference = metabolitePlotFull[order(metabolitePlotFull$Intensity, decreasing = TRUE),][3,]$File
  #this table will contain the alignment reference AND the aligned Files for all
  #Files in each tertile
  bothFull = data.frame(matrix(nrow = 0, ncol = length(colnames(
    metabolitePlotFull
  )) + 1))
  #add a column with the scan adjusted RT's
  # colnames(bothFull) = c(colnames(metabolitePlotFull),"newRT")
  #go through each tertile of the potential signal
  #and subset by that tertile in order to run anovAlign
  for (i in 1:3)
  {
    #print(i)
    #print(" is the zone")
    if (i == 1)
    {
      metabolitePlot = metabolitePlotFull[metabolitePlotFull$Intensity < quantile(metabolitePlotFull$Intensity, .25, na.rm = T), ]
      print(i)
    }
    if (i == 2)
    {
      #metabolitePlot = metabolitePlotFull[metabolitePlotFull$Intensity < quantile(metabolitePlotFull$Intensity, .75),]
      metabolitePlot = metabolitePlotFull[metabolitePlotFull$Intensity < quantile(metabolitePlotFull$Intensity, .75, na.rm =
                                                                                    T) &
                                            metabolitePlotFull$Intensity > quantile(metabolitePlotFull$Intensity, .25, na.rm = T)  , ]
      print(i)
    }
    if (i == 3)
    {
      metabolitePlot = metabolitePlotFull[metabolitePlotFull$Intensity > quantile(metabolitePlotFull$Intensity, .75, na.rm = T), ]
      print(i)
    }
    #for the tertile, pull out the signal belonging to the reference sample to be the reference sample
    metabolitePlotSlice1 =  metabolitePlot[metabolitePlot$File == maxReference, ]
    #initialize the column that will contain the adjusted scan info
    metabolitePlotSlice1$newRT = metabolitePlotSlice1$scan
    #initialize the scan adjusted table, which will contain the original reference data
    #as well as the scan adjust info for each File
    both = metabolitePlotSlice1
    #take out the top of the reference peak
    metabolitePlotSlice1 = metabolitePlotSlice1[metabolitePlotSlice1$Intensity > quantile(metabolitePlotSlice1$Intensity, .75, na.rm = T), ]
    print(paste(maxReference, "is MaxReference"))
    for (i in 1:length(unique(metabolitePlot$File)))
    {
      #print(i)
      if (i != maxReference)
      {
        metabolitePlotSlice2 = metabolitePlot[metabolitePlot$File == unique(metabolitePlot$File)[i], ]
        metabolitePlotSlice2BeforeFilter = metabolitePlotSlice2
        #take the top of the peak to compare to the reference sample as well.  Thus, we are taking the top of each 'zone' in each File
        metabolitePlotSlice2 = metabolitePlotSlice2[metabolitePlotSlice2$Intensity > quantile(metabolitePlotSlice2$Intensity, .75, na.rm = T), ]
        if (length(metabolitePlotSlice1$scan) > 2 &
            length(metabolitePlotSlice2$scan) > 2)
        {
          #print("doing ttest at least")
          #implement the t.test() fxn to determine if there if a shift in RT's
          #between the reference sample and the query
          tTestShift = t.test(metabolitePlotSlice1$scan,
                              metabolitePlotSlice2$scan)
          #subtract the adjusted shift from the RT data for that sample
          metabolitePlotSlice2BeforeFilter$newRT = metabolitePlotSlice2BeforeFilter$scan + (tTestShift$estimate[1] - tTestShift$estimate[2])
          #print("print the estimated RT shift")
          #print((tTestShift$estimate[1] - tTestShift$estimate[2]))
          #print("is the estimated RT shift")
          #create a table recording the adjusted RT shifts for each File
          #in that given zone.
          both = rbind(both,  metabolitePlotSlice2BeforeFilter)
        }
      }
    }
    if (nrow(both) > 0)
    {
      print(nrow(both ))
      print("is nrow of both")
      bothFull = rbind(bothFull, both)
      print("binding both")
      #plot(both$scan, both$Intensity, main = "after adjusting")
      #plot(both$scan, both$Intensity, main = "before adjusting")
    }
  }
  if (nrow(bothFull) < 1)
  {
    bothFull =  metabolitePlotFull
  }
  bothFull$scan = bothFull$newRT
  bothFull$newRT = NULL
  #print(bothFull)
  #make sure that the reajuste scans don't exceed the boundaries of the original data
  bothFull = bothFull[bothFull$scan > minOfSubmittedZone & bothFull$scan < maxOfSubmittedZone,]
  print(minOfSubmittedZone)
  print("is minOfSubmittedZone")
  print(maxOfSubmittedZone)
  print("is maxOfSubmittedZone")
  #return the adjusted table
  return(bothFull)
}




#this is going to be the routine to call all of the intervals
callAllOtherFxns = function(segment1)
{
  
  table = segment1
  #read in the full temp File into a single table
  megaTableII = as.data.table(table)
  megaTableII$File <- str_remove(megaTableII$File, ".rawoutput.txt")
  print("made megatable")
  print(head(megaTableII)) 
  
  
  colnames(megaTableII) =  c("File", "compound","Intensity", "scan")
  megaTableII$compound = as.numeric(megaTableII$compound)
  megaTableII$compound[is.na(megaTableII$compound)] = 0 
  
  myMass = vecOfMasses
  myRT = median(megaTableII$scan)
  print("adding cols")
  
  myCompound = myMass
  
  
  #print(megaTableII)
  
  
  #note that this may cause a crash if we don't have a situation where
  #every File has the Mass ...
  #get the Files present in the File
  myFilesInSubset = unique(megaTableII$File)
  
  #print(i)
  #print("is i")
  #myMass = vecOfMasses[i]
  myMass = as.numeric(myMass)
  
  #starts at the more tolerant ppm threshold, note well
  ppmAllowed = 5
  rangeAllowed = ppmAllowed * .000001 * myMass
  #rangeAllowed = .0001
  #get the information just for that one Mass
  #but, in the FULL set of Files
  #actually, we've already done the subsetting
  mySubsetFull = megaTableII
  
  #randomly choose one of the samples to align against as the central one we'll align against
  #actually, don't do it randomly - choose the sample with the most entries
  
  print("going into callPeaks")
  
  #get the start and the stop info for the regions that we want to sum up
  #using the complete chromatogram
  
  #first run everything through at the File level
  
  #table = read.table(File = "/Users/ahubbard/Downloads/justin_test_data/temp_Files_300_samples_setaria_hilic/temp_Files_more_stringent_filtering/365.107659166667.txt", header = TRUE, sep = ",")
  #table1SubPlot = twoPeaks[twoPeaks$File == "HILIC_Set_94_1_1_H1_85_TB_setaria_12_0511.rawoutput.txt",]
  
  realignedFilesII = megaTableII
  
  #realignedFilesII = generalizeWarping(megaTableII)
  #realignedFilesII = realignedFilesII[realignedFilesII$scan > 0,]
  #print("made it past warping")
  
  #focus only on the scans that are outliers in BOTH Mass and Intensity!
  
  
  #peakSearch = realignedFilesII
  #actually, don't rely on the aligned Files now!
  
  #print( peakSearch)
  #print("is  peakSearch")
  
  
  
  myWindowStartVec = c(min(megaTableII$scan))
  myWindowStopVec = c(max(megaTableII$scan))
  
  #myWindowStartVec = myStartWindows
  #myWindowStopVec = myEndWindows
  
  #have a table that's going to include all of the summed signals for each
  #Mass across each File
  
  #keep track of all of the information for that given Mass-feature across all of the Files
  FilesPresent = c(1:length(myFilesInSubset))
  
  
  #Now we finish the RT window pipeline by summing up the signals
  #across all scans for a given Mass from the aligned chromatograms
  #focusing on all rhe regions id's by shatterFind
  
  print(myWindowStartVec)
  for(coordinates in 1:length(myWindowStartVec))
  {
    myWindowStart =  myWindowStartVec[coordinates]
    myWindowStop = myWindowStopVec[coordinates]
    
    
    #have a table to store all of the information
    allInfoTableFile = matrix(nrow= 0, ncol = length(c("compound","rt", numberOfSampleFiles)))
    allInfoTableFile = as.data.frame(allInfoTableFile)
    
    #note that this could be challenging if we don't have 
    #the Files that we need in the subset
    colnames(allInfoTableFile) =c("compound","rt", as.character(numberOfSampleFiles))
    
    
    #include the Mass info
    allInfoTableFile[1,colnames(allInfoTableFile) == "compound"] = myCompound
    
    #make sure to inclue the RT info as well
    #myRT = mean(myWindowStart, myWindowStop)
    
    
    allInfoTableFile[1,colnames(allInfoTableFile) == "rt"] = myRT
    
    
    # myCompound = unique(mySubset$compound)[1]
    myCompound = myMass
    for(k in 1:length(myFilesInSubset))
    {
      Intensity = 0
      myAdjustedTable = realignedFilesII[realignedFilesII$File == myFilesInSubset[k],]
      print(nrow(myAdjustedTable ))
      print("is the nrow of my table subsetted by the File")
      
      if(nrow(myAdjustedTable) > 5)
      {
        
        toSmoothAndSum = myAdjustedTable
        toSmoothAndSum$File = NULL
        colnames(toSmoothAndSum) = c("compound", "Intensity", "scan")
        #smoothedAddUp = toSmoothAndSum$Intensity
        
        
        #myAdjustedTable = as.data.table(myAdjustedTable)
        subsetToAdd = toSmoothAndSum[toSmoothAndSum$scan > myWindowStart & toSmoothAndSum$scan < myWindowStop, ]$Intensity
        subsetToAdd = subsetToAdd[subsetToAdd != 0]
        smoothedAddUp = subsetToAdd
        
        #quantileTest = quantile(subsetToAdd, c(.25, .85), na.rm = TRUE)
        
        
        if(length(subsetToAdd) > 9)
        {
          
          #print(quantileTest)
          #print(" is quantileTest")
          #print("greater than 10 in length")
          
          #Old Quantile Filter = [smoothedAddUp > quantileTest[1] & smoothedAddUp < quantileTest[2]]
          
          Intensity = sum(smoothedAddUp)
        }
        
        #if we have almost all zeroes, taking the quantile is not going to be useful
        if(length(subsetToAdd) < 10)
        {
          
          Intensity = sum(subsetToAdd)
          Intensity = 0
        }
      }
      
      #if(nrow(myAdjustedTable) < 25)
      #{
      #  Intensity = 0
      #}
      #Intensity = mean(subsetToAdd$Intensity)
      print(Intensity)
      print("is the Intensity")
      allInfoTableFile[1,colnames(allInfoTableFile) == myFilesInSubset[k]] = Intensity
      
      print("going into")
      #print(myFilesInSubset[k])
    }
    
    #if(sigRatioFirst > 3)
    #{
    
    #allInfoTableFile[is.na(allInfoTableFile)] = 0
    #allInfoTableallFiles = rbind(allInfoTableallFiles, allInfoTableFile)
    #}	
  }
  
  
  #make sure to turn all NA's to 0
  allInfoTableFile[is.na(allInfoTableFile)] = 0
  #print( allInfoTableFile)
  print(head(allInfoTableFile))
  print("is what we are outputting as a table")
  
  
  return(allInfoTableFile)
}


#mySubset = read.table(File = "/Users/ahubbard/Downloads/justin_test_data/temp_Files_300_samples_setaria_hilic/temp_Files_more_stringent_filtering/365.107659166667.txt", header = TRUE, sep = ",")
#allFiles = read.table(File = "/Users/ahubbard/Downloads/justin_test_data/temp_Files_300_samples_setaria_hilic/temp_Files_more_stringent_filtering/365.107659166667.txt", header = TRUE, sep = ",")
#allFiles = read.table(file = "/Users/ahubbard/Downloads/176.10273merged.txt", header = TRUE, sep = "\t")
#allFiles = read.table(file = "/Users/ahubbard/Downloads/176.10273merged.txt", header = TRUE, sep = "\t")

#pipeline for bounds goes here: outputs table where colnames are Files and first row is Intensity (first three columns are compound(Mass), rt, x
#Before running, need to initialize output table:
FileList = read.table(args[2], header = T)$x
FileList = str_remove(FileList, ".rawoutput.txt")

nameOfallFiles = FileList
numberOfSampleFiles = nameOfallFiles

numberOfSampleFiles = unique(numberOfSampleFiles)

allInfoTableallFiles = matrix(nrow= 0, ncol = length(c("compound","rt", numberOfSampleFiles)))
allInfoTableallFiles = as.data.frame(allInfoTableallFiles)

#note that this could be challenging if we don't have 
#the Files that we need in the subset
colnames(allInfoTableallFiles) =c("compound","rt", as.character(numberOfSampleFiles))




##finish changing the allFiles!
#allFiles = read.delim("/Users/ahubbard/Downloads/176.10273merged.txt")
#allFiles = read.delim("/Users/ahubbard/176.10269merged.txt")




#make sure that the Intensity column is numeric
allFiles$Intensity = as.numeric(allFiles$Intensity)
#replace any NA's with 0's
allFiles$Intensity[is.na(allFiles$Intensity)] = 0



#for now, we won't need the Mstype or the polarity info
#if it is in the table, remove it!
allFiles$polarity = NULL
allFiles$Massspectype = NULL
allFiles$newMass = NULL
allFiles = allFiles[allFiles$Mstype != 'Ms2',]
allFiles$Mstype = NULL
#set colnames to make them consistent with other code
colnames(allFiles) = c("File", "Mass", "Intensity", "scan")
#currently unused, will store final results after running anovalign
allFiles = allFiles[order(allFiles$File),]




allFiles$scan = as.numeric(allFiles$scan)
print(head(allFiles$scan))
print(" is allFiles$scan")

#remove the region of the retention time data referred to as "the void" - typically only low quality signals exist in this region
allFiles = allFiles %>% filter(scan < 100)
length(allFiles$scan)

#Convert Scans From Minutes to milliseconds and round them:
allFiles$scan = round(allFiles$scan * 60 * 1000, digits = 0)

#Bin the "scans" into 3000 buckets. This way each "scan bucket" contains as much information as one "scan" from before.




#Bin the "scans" into 3000 buckets. This way each "scan bucket" contains as much information as one "scan" from before.
scanMax = max(allFiles$scan)
cutValue = scanMax / 3000
toCut = ((scanMax)/600)



allFiles$scan = cut(allFiles$scan, breaks = (max(allFiles$scan)/cutValue))
levels(allFiles$scan) = seq(1,scanMax) 
allFiles$scan = as.numeric(allFiles$scan)


allFiles = allFiles[allFiles$scan > 200,]


#call the boundaryFinder fxn now
bounds = boundaryFinder(allFiles, allFiles)

vecOfStarts = bounds[[1]]
vecOfStops = bounds[[2]]


#Ordering the starts and stops to make the subsequent merging work
vecOfStarts <- vecOfStarts[order(vecOfStarts, decreasing = F)]
vecOfStops <- vecOfStops[order(vecOfStops, decreasing = F)]

#initialize vectors for final bounds to be plotted
vecOfStartsFinal <- vector()
vecOfStopsFinal <- vector()
listOfFinalBounds <- list()

#collapses those bounds that are within 200 scans of one another
#boundCollapser - if scans are w/in 200 of each other, collapse

#returns a list of two vectors with any starts/stops within 200 scans merged.
listOfFinalBounds <- boundCollapser(vecOfStarts, vecOfStops)
vecOfStartsFinal <- listOfFinalBounds[[1]]
vecOfStopsFinal <- listOfFinalBounds[[2]]


table1Original = allFiles


print(head(table1Original))
table1Original$polarity = NULL
table1Original$Massspectype = NULL
table1Original$Mstype = NULL
colnames(table1Original) = c("File", "Mass", "Intensity", "scan")




print(head(table1Original))

if(!(is.numeric(vecOfStartsFinal))){
  vecOfStartsFinal = 0
}
if(!(is.numeric(vecOfStopsFinal))){
  vecOfStopsFinal = 3000
}

#table = tableOriginal
#now, go through each of the "zones" and determine if they need to be merged or not
for(i in 1:length(vecOfStartsFinal))
{
  
  start1 = vecOfStartsFinal[i]
  stop1 =  vecOfStopsFinal[i]
  print(start1)
  print(" is start1")  
  print(stop1)
  print(" is stop1")
  print(" is now looking up the bounds !!!")
  
  
  segment1 = table1Original[table1Original$scan > start1 & table1Original$scan < stop1,]
  
  #make sure that we have more signal in our region than noise
  outsideSegment1 = table1Original[table1Original$scan < start1,]
  outsideSegment2 = table1Original[ table1Original$scan > stop1,]
  
  
  outsideSegment = rbind(outsideSegment1 , outsideSegment2)

  #I think this is causing the issue
  
  sigToNoise = (sum(segment1$Intensity) / (max(segment1$scan) - min(segment1$scan))) / (sum(outsideSegment$Intensity) / (max(outsideSegment$scan) - min(outsideSegment$scan)))
  if(sigToNoise > 2)
  {
    
    #segment1 = segment1[segment1$Intensity > quantile(segment1$Intensity, p = .9),]
    print(head(segment1))
    segment1 = generalizeWarping(segment1)
    
    #run anovAlign on the segment now! We want to re-run boundary finder here to get tighter bounds
    
    #colnames(segment1) = c("File", "Mass", "Intensity", "scan")
    segment1$scan = round(segment1$scan)
    
    bounds = boundaryFinder(segment1, allFiles)
   
    #if don't improve bounds, keep the old ones
    if(length(bounds[[1]]) == 0)
    {
      bounds[[1]] = start1
    }
    
    if(length(bounds[[2]]) == 0)
    {
      bounds[[2]] = stop1
    }
    
    
    bounds = boundCollapser(bounds[[1]], bounds[[2]])
    
    
    if( abs(bounds[[2]] - (bounds[[1]]) < 100 ))
    {
      bounds[[1]] = start1  
      bounds[[2]] = stop1   
    }
    
    segment1 = segment1[segment1$scan > bounds[[1]] & segment1$scan < bounds[[2]],]
    anovalignLeftBound = bounds[[1]][i]
    anovalignRightBound = bounds[[2]][i]
    
    
    tableForRegion = callAllOtherFxns(segment1)
    allInfoTableallFiles = rbind(allInfoTableallFiles,tableForRegion)
  }
}



write.csv(allInfoTableallFiles, file =paste(args[4], "/", vecOfMasses[1], ".txt", sep = ""), row.names=F)


name = paste(args[3], '/', vecOfMasses, ".pdf", sep = "")

print("here is segment1:")
print(head(segment1))


pdf(name)


if(exists('anovalignLeftBound')) {
print("test:")
print(anovalignLeftBound)
print(anovalignRightBound)



plot(segment1$scan, segment1$Intensity, main = "ANOVALIGN", xlim=c(0,3000))
abline(v = anovalignLeftBound, col = "green")
abline(v = anovalignRightBound, col = "green")  
}



#write.csv(allInfoTableallFiles, file =paste(args[4], "/", vecOfMasses, ".txt", sep = ""), row.names=F)



#Now we have the allInfoTableallFiles
#cols = compound, rt, [Files in order of injection]
#row1 = sum of Intensity within RT bounds +- 2.5ppm

print(head(allFiles))
print(bounds)

plot(allFiles$Mass, allFiles$Intensity, main = "10PPM Window")
abline(v = myMassTest, col = "green")

plot(allFiles$scan, allFiles$Intensity)
abline(v = bounds[1], col = "green")
abline(v = bounds[2], col = "green")

dev.off()


