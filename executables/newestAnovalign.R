library(data.table)
library(dplyr)
library(stringr)


args = commandArgs(trailingOnly=TRUE)

#all setup within this brace: 
{

vecOfMasses = basename(args[1])
currentDir = args[1]


print(vecOfMasses)
myMassTest = as.numeric(vecOfMasses)



myMassTest = myMassTest
myMassTest = round(myMassTest, digits = 5)
#print(myMassTest)
myMassTestQuery = myMassTest


MassSubDirectory = currentDir
#print(MassSubDirectory)
#print(" is where we're going to read from")

File_list <- list.files(path=MassSubDirectory, full.names = TRUE)
File_list = File_list[sapply(File_list, file.size) > 0]

fullFileName = paste(currentDir, "/",  vecOfMasses,"merged.txt", sep = "")
mergedExclude = fullFileName

#print(fullFileName)

#print(File_list)
#print(!(fullFileName %in% File_list))


#simply create the merged file!
if(!(fullFileName %in% File_list))
{
  
  
  allFiles = read.table(File_list[1], header = TRUE, sep = ",")
  

  #print(allFiles)
  #print(" is allFiles")
  
  File_list = File_list[-1]
  if(length(File_list) > 0)
  {
    for (j in 1:length(File_list))
    {
      print(j)
      print("is j")
      if(File_list[j] != mergedExclude)
      {   
        
        oneFile = read.table(File_list[j], header = TRUE, sep = ",")
        #print(head(oneFile))
        allFiles = rbind(allFiles, oneFile)
        #print(j)
        #print(myMassTest)
      }
    }
  }
  #print("we're writing the output now")
  #write to the temp table
  
  write.table(allFiles, file=fullFileName, row.names = FALSE, sep = ",")
  
}else {
  
  allFiles = read.table(fullFileName, header = TRUE, sep = ",")
  
}	

}





#print(head(allFiles))
#all functions within this brace:
{



#Pipeline for Mass features:


boundaryFinder = function(mySubset)
{
  
  allFiles = mySubset
  allInfoScans <- allFiles %>% count(scan, sort = TRUE) 
  
  #get all the information for the full plausible mass feature region
  #and then filter out the lower 50 % of scans by number of signals
  allInfoScans = allInfoScans[allInfoScans$n > quantile(allInfoScans$n, p = .5),]$scan
  allFiles = allFiles[allFiles$scan %in% allInfoScans,]
  
  #partition the EIC into division of 5 scans
  divider = 5
  numerator = 0
  
  #figure out how many times through we'll iterate
  iteration = max(allFiles$scan) /divider
  
  #keep track of the correlations for each window
  vecOfCorrelations = vector()
  vecOfRegions = vector()
  vecOfIndexes = vector()
  sigVNoise = vector()
  
  rowNums = 0
  
  #iterate through each partition
  for(i in 1:iteration)
  {
    
    start = i * divider
    stop = start + divider
    

    #table one is going to be one window
    table1BoundA = start
    table1BoundB = stop 
    
    
    #table two is going to another window
    table2BoundA = start - 1
    table2BoundB = table2BoundA - divider
    
    #create the table for the first window
    table1half = allFiles[allFiles$scan <= stop,]
    table1half = table1half[table1half$scan > start,]
    
    #create the table for the second window
    table2half = allFiles[allFiles$scan <= table2BoundA,]
    table2half = table2half[table2half$scan > table2BoundB,]
    
  
    
    
    #when calculating signal to noise and other qualities of our
    #EIC regions, we're going to make sure to use the table
    #that has the most signal, by default
    
    #if window2's table has more signal, we'll use it for the statistical information
    if(sum(table2half$Intensity) > sum(table1half$Intensity))
    {
      
      numerator = sum(table2half$Intensity)
      rowNums = nrow(table2half)
      maxSub = max(table2half$Intensity)
      
    }
    
    
    #if window1's table has more signal, we'll use it for the statistical information
    if(sum(table1half$Intensity) > sum(table2half$Intensity))
    {
      
      numerator = sum(table1half$Intensity)
      rowNums = nrow(table1half)
      maxSub = max(table1half$Intensity)
      
    }
    
    #calculate the signal/noise for the section of the EIC that we're on
    #take the most 'optimistic' perspective, by using the one of the 
    #two windows with the most signal
    sigToNoise =  numerator / (mean(allFiles$Intensity)*rowNums)
    
    
    #if we've got more signal than noise
    #and each of the tables for my window have signal, this region of the EIC
    #could contain a metabolite
    if(nrow(table2half) > 0 & nrow(table1half) > 0 & (sigToNoise > 1))
    {
      
      #keep track of the sig/noise for comparison purposes
      sigVNoise = c(sigVNoise, sigToNoise)
      
      # an summarize the window
      #by keepign the maximum signal for each file
      Grouped1B = table1half %>% group_by(File) %>% summarise(Value = max(Intensity))
      Grouped2B = table2half  %>% group_by(File) %>% summarise(Value = max(Intensity))
      
      #keep only the information for files that have data in each window
      Grouped1B = Grouped1B[Grouped1B$File %in% intersect(Grouped1B$File, Grouped2B$File),]
      Grouped2B = Grouped2B[Grouped2B$File %in% intersect(Grouped1B$File, Grouped2B$File),]
      
      #order the tables for each scan window by file
      Grouped1B = Grouped1B[order(Grouped1B$File),]
      Grouped2B = Grouped2B[order(Grouped2B$File),]
      
      #now, each window has the same number of signals (only kept signals associated with scans in each file)
      #so they are the same length, and each vector has been ordered by file
      #thus, we can calculate the correlation to see if the signal
      #across each window correlates
      theCor = cor(Grouped2B$Value, Grouped1B$Value)
      
      
      #if signal in both scan windows correlates, this section of the EIC
      #could maintain a metabolite
      if((is.na(theCor) == FALSE) & theCor > .5)
      {
        #print(theCor)
        vecOfCorrelations = c(vecOfCorrelations, theCor)
        vecOfRegions = c(vecOfRegions, start)
        vecOfIndexes = c(vecOfIndexes, i)
      }
    }
  }
  
  return(vecOfRegions)
  
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
  
  #print(head(metabolitePlotFull[order(metabolitePlotFull$Intensity, decreasing = TRUE),]))
  
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

    #get the first tertile
    if (i == 1)
    {
      metabolitePlot = metabolitePlotFull[metabolitePlotFull$Intensity <= quantile(metabolitePlotFull$Intensity, .25, na.rm = T), ]
      #print(i)
    }
    #get the second tertile
    if (i == 2)
    {
      #metabolitePlot = metabolitePlotFull[metabolitePlotFull$Intensity < quantile(metabolitePlotFull$Intensity, .75),]
      metabolitePlot = metabolitePlotFull[metabolitePlotFull$Intensity <= quantile(metabolitePlotFull$Intensity, .75, na.rm =
                                                                                    T) &
                                            metabolitePlotFull$Intensity >= quantile(metabolitePlotFull$Intensity, .25, na.rm = T)  , ]
      #print(i)
    }
    #get the third tertile
    if (i == 3)
    {
      metabolitePlot = metabolitePlotFull[metabolitePlotFull$Intensity >= quantile(metabolitePlotFull$Intensity, .75, na.rm = T), ]
      #print(i)
    }
    
    #for the tertile, pull out the signal belonging to the reference sample to be the reference sample
    metabolitePlotSlice1 =  metabolitePlot[metabolitePlot$File == maxReference, ]
    #initialize the column that will contain the adjusted scan info
    metabolitePlotSlice1$newRT = metabolitePlotSlice1$scan
    #initialize the scan adjusted table, which will contain the original reference data
    #as well as the scan adjust info for each File
    both = metabolitePlotSlice1
    #take out the top of the reference peak
    
    #here we do the filtering of the reference sample, but maybe we won't want to do this
    #metabolitePlotSlice1 = metabolitePlotSlice1[metabolitePlotSlice1$Intensity >= quantile(metabolitePlotSlice1$Intensity, .75, na.rm = T), ]
    #print(paste(maxReference, "is MaxReference"))
    for (i in 1:length(unique(metabolitePlot$File)))
    {
      print(i)
      if (i != maxReference)
      {
        metabolitePlotSlice2 = metabolitePlot[metabolitePlot$File == unique(metabolitePlot$File)[i], ]
        metabolitePlotSlice2BeforeFilter = metabolitePlotSlice2
        
        #take the top of the peak to compare to the reference sample as well.  Thus, we are taking the top of each 'zone' in each File
        #metabolitePlotSlice2 = metabolitePlotSlice2[metabolitePlotSlice2$Intensity > quantile(metabolitePlotSlice2$Intensity, .75, na.rm = T), ]
        if (length(metabolitePlotSlice1$scan) > 1 &
            length(metabolitePlotSlice2$scan) > 1)
        {
          #print("doing ttest at least")
          #implement the t.test() fxn to determine if there if a shift in RT's
          #between the reference sample and the query
          tTestShift = t.test(metabolitePlotSlice1$scan,
                              metabolitePlotSlice2$scan)
          
          #subtract the adjusted shift from the RT data for that sample
          metabolitePlotSlice2BeforeFilter$newRT = metabolitePlotSlice2BeforeFilter$scan + (tTestShift$estimate[1] - tTestShift$estimate[2])

          #create a table recording the adjusted RT shifts for each File
          #in that given zone.
          both = rbind(both, metabolitePlotSlice2BeforeFilter)
        }
        
        #if there's only one data point - not enough to do a t-test, but no worries.
        #we'll just add the point back to the table
        if (length(metabolitePlotSlice2$scan) == 1)
        {
          
          metabolitePlotSlice2BeforeFilter$newRT = metabolitePlotSlice2BeforeFilter$scan
          both = rbind(both,  metabolitePlotSlice2BeforeFilter)
          
        }
      }
    }
    if (nrow(both) > 0)
    {
      #print(nrow(both ))
      #print("is nrow of both")
      bothFull = rbind(bothFull, both)
      #print("binding both")
    }
  }
  if (nrow(bothFull) < 1)
  {
    bothFull =  metabolitePlotFull
  }
  bothFull$scan = bothFull$newRT
  bothFull$newRT = NULL

  
  #perhaps we've extended beyond original bounds ... we may trim this later
  #make sure that the reajuste scans don't exceed the boundaries of the original data

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
  #print("adding cols")
  
  myCompound = myMass
  
  
  #print(megaTableII)
  
  
  #note that this may cause a crash if we don't have a situation where
  #every File has the Mass ...
  #get the Files present in the File
  myFilesInSubset = unique(megaTableII$File)
  
  #print(i)
  #print("is i")
  myMass = as.numeric(myMass)
  
  ppmAllowed = 5
  rangeAllowed = ppmAllowed * .000001 * myMass
  #get the information just for that one Mass
  mySubsetFull = megaTableII
  
  realignedFilesII = megaTableII
  
  myWindowStartVec = c(min(megaTableII$scan))
  myWindowStopVec = c(max(megaTableII$scan))
  
  
  FilesPresent = c(1:length(myFilesInSubset))
  
  
  #Now we finish the RT window pipeline by summing up the signals
  #across all scans for a given Mass from the aligned chromatograms
  
  #print(myWindowStartVec)
  for(coordinates in 1:length(myWindowStartVec))
  {
    myWindowStart =  myWindowStartVec[coordinates]
    myWindowStop = myWindowStopVec[coordinates]
    
    
    #have a table to store all of the information
    allInfoTableFile = matrix(nrow= 0, ncol = length(c("compound","rt", numberOfSampleFiles)))
    allInfoTableFile = as.data.frame(allInfoTableFile)
    
    colnames(allInfoTableFile) =c("compound","rt", as.character(numberOfSampleFiles))
    
    
    #include the Mass info
    allInfoTableFile[1,colnames(allInfoTableFile) == "compound"] = myCompound
    
    #make sure to inclue the RT info as well
    
    allInfoTableFile[1,colnames(allInfoTableFile) == "rt"] = myRT
    
    
    myCompound = myMass
    for(k in 1:length(myFilesInSubset))
    {
      Intensity = 0
      myAdjustedTable = realignedFilesII[realignedFilesII$File == myFilesInSubset[k],]
      #print(nrow(myAdjustedTable ))
      #print("is the nrow of my table subsetted by the File")
      
      if(nrow(myAdjustedTable) > 0)
      {
        
        toSmoothAndSum = myAdjustedTable
        toSmoothAndSum$File = NULL
        colnames(toSmoothAndSum) = c("compound", "Intensity", "scan")
        
        subsetToAdd = toSmoothAndSum[toSmoothAndSum$scan > myWindowStart & toSmoothAndSum$scan < myWindowStop, ]$Intensity
        subsetToAdd = subsetToAdd[subsetToAdd != 0]
        smoothedAddUp = subsetToAdd
        
        if(length(subsetToAdd) > 0)
        {
          
          #print(quantileTest)
          #print(" is quantileTest")
          #print("greater than 10 in length")
          
          #Old Quantile Filter = [smoothedAddUp > quantileTest[1] & smoothedAddUp < quantileTest[2]]
          
          Intensity = sum(smoothedAddUp)
        }
        
        #if we have almost all zeroes, taking the quantile is not going to be useful
        if(length(subsetToAdd) == 0)
        {
          
          Intensity = 0
        }
      }
      
      #print(Intensity)
      #print("is the Intensity")
      allInfoTableFile[1,colnames(allInfoTableFile) == myFilesInSubset[k]] = Intensity
      
    }
    
  }
  
  
  #make sure to turn all NA's to 0
  allInfoTableFile[is.na(allInfoTableFile)] = 0
  #print(head(allInfoTableFile))
  #print("is what we are outputting as a table")
  
  
  return(allInfoTableFile)
}

}


#pipeline for bounds goes here: outputs table where colnames are Files and first row is Intensity (first two columns are compound(Mass), rt
#Before running, need to initialize output table:
FileList = read.table(args[2], header = T)$x




FileList = str_remove(FileList, ".rawoutput.txt")
nameOfallFiles = FileList
numberOfSampleFiles = nameOfallFiles
numberOfSampleFiles = unique(numberOfSampleFiles)
allInfoTableallFiles = matrix(nrow= 0, ncol = length(c("compound","rt", numberOfSampleFiles)))
allInfoTableallFiles = as.data.frame(allInfoTableallFiles)

colnames(allInfoTableallFiles) =c("compound","rt", as.character(numberOfSampleFiles))




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

allFiles = allFiles[order(allFiles$File),]




allFiles$scan = as.numeric(allFiles$scan)
#print(head(allFiles$scan))
#print(" is allFiles$scan")

#print(max(allFiles$scan))
#remove the region of the retention time data referred to as "the void" - typically only low quality signals exist in this region
allFiles = allFiles %>% filter(scan < 100)
#print(paste(length(allFiles$scan), 'is length of scans'))

#Convert Scans From Minutes to milliseconds and round them:
allFiles$scan = round(allFiles$scan * 60 * 1000, digits = 0)
#Bin the "scans" into 3000 buckets. This way each "scan bucket" contains as much information as one "scan" from before.




#Bin the "scans" into 3000 buckets. This way each "scan bucket" contains as much information as one "scan" from before.
scanMax = max(allFiles$scan)
#print(scanMax)
cutValue = scanMax / 3000
toCut = ((scanMax)/600)



allFiles$scan = cut(allFiles$scan, breaks = (max(allFiles$scan)/cutValue))
levels(allFiles$scan) = seq(1,scanMax) 
allFiles$scan = as.numeric(allFiles$scan)


allFiles = allFiles[allFiles$scan > 200,]



allFilesList = unique(allFiles$File) 



dataPerFile = allFiles %>% count(File)


middleOfSignal = round(length(unique(allFiles$File) )/ 2)


#top 3rd file in terms of intensity will be reference file
referenceFile = allFiles[order(allFiles$Intensity, decreasing = TRUE),][3,]$File
testSignal = allFiles[allFiles$File %in% referenceFile,  ]
#plot(testSignal$scan, testSignal$Intensity)



#get the reference signals from the top 5 of the reference
referenceSignals = mean(testSignal$scan[order(testSignal$Intensity, decreasing = TRUE)[1:5]])



toAdjust = allFiles

for(i in 1:length(allFilesList))
{
  #print(i)
  
  fileName = allFilesList[i]
  
  
  #implement the pre-correction
  if(fileName != referenceFile )
  {
    
    #get the top few signals from the main file
    toQuery = allFiles[allFiles$File %in% fileName,  ]
    querySignals = mean(toQuery$scan[order(toQuery$Intensity, decreasing = TRUE)[1:5]])
    
    
    differenceFromReference = referenceSignals -  querySignals
    #print(differenceFromReference)
    
    #adjust the scans now
    toAdjust[toAdjust$File == fileName,]$scan =  toAdjust[toAdjust$File == fileName,]$scan - differenceFromReference
    
  }
  
}


toAdjust = na.omit(toAdjust)

correlatedRegions = boundaryFinder(toAdjust)

bounds = list()


bounds[[1]] = correlatedRegions 
bounds[[2]] = correlatedRegions + 10

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


#print(head(table1Original))
table1Original$polarity = NULL
table1Original$Massspectype = NULL
table1Original$Mstype = NULL
colnames(table1Original) = c("File", "Mass", "Intensity", "scan")


if(!(is.numeric(vecOfStartsFinal))){
  vecOfStartsFinal = 0
}
if(!(is.numeric(vecOfStopsFinal))){
  vecOfStopsFinal = 3000
}

#now that we've got bounds, loop through them all
#in order to call the signal for the mass feature and add it to the table
for(i in 1:length(vecOfStartsFinal))
{
  
  start1 = vecOfStartsFinal[i]
  stop1 =  vecOfStopsFinal[i]

  #for each set of bounds, subset the potential EIC by the bounds 
  segment1 = table1Original[table1Original$scan > start1 & table1Original$scan < stop1,]
  
  #make sure that we have more signal in our region than noise
  #we'll do this by looking at the region outside the bounds of our potential metabolite
  outsideSegment1 = table1Original[table1Original$scan < start1,]
  outsideSegment2 = table1Original[ table1Original$scan > stop1,]
  
  #bind together the regions OUTSIDE our potential mass features
  outsideSegment = rbind(outsideSegment1 , outsideSegment2)
  

  #calculate the signal to noise as signal / (number of scans) for each region (inside the scans AND outside)
  #note that the outside region will be intentionally diluted
  sigToNoise = (sum(segment1$Intensity) / (max(segment1$scan) - min(segment1$scan))) / (sum(outsideSegment$Intensity) / (max(outsideSegment$scan) - min(outsideSegment$scan)))
  
  if(sigToNoise > 2)
  {
    
    #print(head(segment1))
    
    #only use anovalign if actually a large signal that could show high drift
    #if there isn't really any signal there, anyways, we'll just add up the signal
    if(nrow(segment1) > 100)
    {
      segment1 = generalizeWarping(segment1)
    }
    #run anovAlign on the segment now! We want to re-run boundary finder here to get tighter bounds
    
    segment1$scan = round(segment1$scan)
    
    segment1 = na.omit(segment1)
    bounds = boundaryFinder(segment1)
   
    #if don't improve bounds (if they end up as zero), keep the old ones
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



write.csv(allInfoTableallFiles, file =paste(args[3], "/", vecOfMasses[1], ".txt", sep = ""), row.names=F)




