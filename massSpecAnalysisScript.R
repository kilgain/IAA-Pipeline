library(dplyr)
library(data.table)

sink("MassSpecOutput.txt")

#bounds- left and righthand bounds to evaluate
#myData- matrix of data for all samples with no points removed
highQualityPeaks = function(leftBounds, rightBounds, myData) {
  vecOfBoundsToKeep = vector()
  leftBounds = leftBounds
  rightBounds = rightBounds
  for (i in 1:length(leftBounds)) {
    keepPlot = TRUE
    leftHand = leftBounds[i]
    rightHand = rightBounds[i]
    Data = myData
    Data = Data[Data$scan > leftHand & Data$scan < rightHand, ]
    
    # myDataUpper = myData[myData$intensity > (2/3)*max(myData$intensity),3]
    
    DataMiddle = Data[Data$intensity >  (.4) * max(Data$intensity) &
                            Data$intensity < (3 / 4) * max(Data$intensity), ]
    # myDataBottom = myData[myData$intensity < (1/3)*max(myData$intensity),3]
    
    
    if (length(unique(DataMiddle$file)) > (1 / 2) * length(unique(Data$file))) {
      vecOfBoundsToKeep <- c(vecOfBoundsToKeep, keepPlot)
      next
    }
    else{
      keepPlot = FALSE
      vecOfBoundsToKeep <- c(vecOfBoundsToKeep, keepPlot)
      next
    }
    
  }
  return(vecOfBoundsToKeep)
}

#Will take in data matrix (if two peaks the points within the other peak will already be removed) and bounds
#will return the signal/noise value and maximum intensity value within the bounds provided
signalToNoiseCalculator = function(dataMatrix, bounds) {
  myData <- dataMatrix
  leftHand <- bounds[1]
  rightHand <- bounds[2]
  signalRegion = myData[myData$scan > leftHand &
                          myData$scan < rightHand, ]
  signal = mean(signalRegion$intensity[signalRegion$intensity > quantile(signalRegion$intensity, probs =
                                                                           .33) &
                                         signalRegion$intensity < quantile(signalRegion$intensity, probs = .67)])
  signalMax = max(signalRegion$intensity)
  if(signalMax < 0) {
    signalMax = 0
  }
  #exclude anything within this interval from background calculation
  intervalExclude = c(myData$scan > leftHand &
                        myData$scan < rightHand)
  
  #all points outside bounds of any boundaries
  mySubset1Exclude = myData[!intervalExclude, ]
  
  mySubset1Exclude$intensity[is.na(mySubset1Exclude$intensity)] = 0
  # signalBackground = mean(mySubset1Exclude$intensity[mySubset1Exclude$intensity >
  #                                                      1000], na.rm = TRUE)
  signalBackground <-
    mean(mySubset1Exclude$intensity[mySubset1Exclude$intensity < quantile(mySubset1Exclude$intensity, probs = .10)])
  #prevents errors when all points are within bounds, and you cannot calculate a background
  if (is.nan(signalBackground)) {
    signalBackground = 1
  }
  if (is.na(signalBackground)) {
    signalBackground = 1
  }
  if (is.nan(signal)) {
    signal = 1
  }
  if (is.na(signal)) {
    signal = 1
  }
  signalToNoiseRatio = signal / signalBackground
  return(c(signalToNoiseRatio, signalMax))
}




#boundCollapser - if scans are w/in 200 of each other, collapse
#vecOfStarts = left hand bounds
#vecOfStops = right hand bounds
#returns a list of two vectors with any starts/stops within 200 scans merged.
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
  #maxReference = dplyr::count(metabolitePlotFull, file)[order(dplyr::count(metabolitePlotFull, file)$n, decreasing = TRUE), ][3, ]$file
  maxReference = metabolitePlotFull[order(metabolitePlotFull$intensity, decreasing = TRUE), ][3, ]$file
  #this table will contain the alignment reference AND the aligned files for all
  #files in each tertile
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
      metabolitePlot = metabolitePlotFull[metabolitePlotFull$intensity < quantile(metabolitePlotFull$intensity, .25, na.rm = T),]
      # print(i)
    }
    if (i == 2)
    {
      #metabolitePlot = metabolitePlotFull[metabolitePlotFull$Intensity < quantile(metabolitePlotFull$Intensity, .75),]
      metabolitePlot = metabolitePlotFull[metabolitePlotFull$intensity < quantile(metabolitePlotFull$intensity, .75, na.rm =
                                                                                    T) &
                                            metabolitePlotFull$intensity > quantile(metabolitePlotFull$intensity, .25, na.rm = T)  ,]
      # print(i)
    }
    if (i == 3)
    {
      metabolitePlot = metabolitePlotFull[metabolitePlotFull$intensity > quantile(metabolitePlotFull$intensity, .75, na.rm = T),]
      # print(i)
    }
    #for the tertile, pull out the signal belonging to the reference sample to be the reference sample
    metabolitePlotSlice1 =  metabolitePlot[metabolitePlot$file == maxReference,]
    #initialize the column that will contain the adjusted scan info
    metabolitePlotSlice1$newRT = metabolitePlotSlice1$scan
    #initialize the scan adjusted table, which will contain the original reference data
    #as well as the scan adjust info for each file
    both = metabolitePlotSlice1
    #take out the top of the reference peak
    metabolitePlotSlice1 = metabolitePlotSlice1[metabolitePlotSlice1$intensity > quantile(metabolitePlotSlice1$intensity, .75, na.rm = T),]
    for (i in 1:length(unique(metabolitePlot$file)))
    {
      #print(i)
      if (i != maxReference)
      {
        metabolitePlotSlice2 = metabolitePlot[metabolitePlot$file == unique(metabolitePlot$file)[i],]
        metabolitePlotSlice2BeforeFilter = metabolitePlotSlice2
        #take the top of the peak to compare to the reference sample as well.  Thus, we are taking the top of each 'zone' in each file
        metabolitePlotSlice2 = metabolitePlotSlice2[metabolitePlotSlice2$intensity > quantile(metabolitePlotSlice2$intensity, .75, na.rm = T),]
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
          #create a table recording the adjusted RT shifts for each file
          #in that given zone.
          both = rbind(both,  metabolitePlotSlice2BeforeFilter)
        }
      }
    }
    if (nrow(both) > 0)
    {
      # print(nrow(both ))
      # print("is nrow of both")
      bothFull = rbind(bothFull, both)
      # print("binding both")
      # plot(both$scan, both$intensity, main = "after adjusting")
      # plot(both$scan, both$intensity, main = "before adjusting")
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
  bothFull = bothFull[bothFull$scan > minOfSubmittedZone &
                        bothFull$scan < maxOfSubmittedZone, ]
  # print(minOfSubmittedZone)
  # print("is minOfSubmittedZone")
  # print(maxOfSubmittedZone)
  # print("is maxOfSubmittedZone")
  #return the adjusted table
  return(bothFull)
}

#expands boundaries
#MySubset -> original dataset used to extract values of intensity/scan for signal/noise calculation
#SignalToNoiseRatio -> used as starting value for signalToNoise... could calculate this within function
#middleOfPeak -> starting point of bound expansion
#peakSD ->  used for amount of expansion per iteration
#vecOfStarts -> used to check if expanding bounds rub up against established bounds
#vecOfStops -> used to check if expanding bounds rub up against established bounds
boundaryExpander =
  function(mySubset,
           signalToNoiseRatio,
           middleOfPeak,
           peakSD,
           vecOfStarts = vector(),
           vecOfStops = vector()) {
    #initialization of variables
    vecOfStarts <- vecOfStarts
    vecOfStops <- vecOfStops
    leftHandOld = 0
    rightHandOld = 0
    
    #used for identifying if a current bound falls within an established bound
    range <- cbind(vecOfStarts, vecOfStops)
    
    #Want to have a condition to exit if we have hit a left and right bound but not met any other condition to exit
    boundhit1 = FALSE
    boundhit2 = FALSE
    
    #Want to have an iteration limit to prevent over-expansion of bounds
    counter = 0
    
    #initialization
    mySubset <- mySubset
    signalToNoiseRatio <- signalToNoiseRatio
    StartingSignalToNoiseRatio <- signalToNoiseRatio
    middleOfPeak <- middleOfPeak
    peakSD <- peakSD
    
    #This statement should only execute if there are no points left to use in calculation of background
    #causing signal/noise to become x/0 -> NaN
    if (is.na(signalToNoiseRatio) | is.nan(signalToNoiseRatio)) {
      return(c(0, 3000))
    }
    

    #This is where the expansion of bounds occurs, the outer loop condition is that if the signal/noise drops below
    #90% of the original value, exit. This value was chosen for performance.
    while (signalToNoiseRatio > .9 * StartingSignalToNoiseRatio)
    {
      #if we have established that the right bound is within another boundary, we dont want to keep expanding it
      #therefore, only expand if boundhit1 (right) is FALSE
      if (boundhit1 != TRUE) {
        rightHand = middleOfPeak + (2 + .15 * counter) * peakSD
      }
      
      #If we happen to expand the right-hand boundary over the maximum allowed value, set it to the highest value.
      #I did have this at 3000, but the highest scan in the data is 3051.
      
      if (rightHand > 3100) {
        rightHand = 3100
        boundhit1 = TRUE
      }
      
      # This chunk of code checks that the RIGHTHAND expanded bound is not within an already established boundary.
      # If it is, boundhit1 (right) is set to TRUE and the right bound is unified to the lefthand bound of the established boundary
      # differs from boundaryChecker in that this must also return boundhit=TRUE
      if (length(range[, 1]) != 0) {
        for (j in 1:length(range[, 1])) {
          if (rightHand %between% range[j,]) {
            # print(
            #   paste(
            #     "rightHand, ",
            #     rightHand,
            #     "is inside range of previous boundary: ",
            #     range[j, ]
            #   )
            # )
            boundhit1 = TRUE
            #   rightHand <- closest value in vecofStarts OR vecOfstops (NOTE: I added this BEFORE adding the same functionality
            # in multiPeakEvaluator(), where the right bound is always unified to an existing left hand bound to prevent overlap.)

            # if (abs(rightHand - range[j, 1]) < abs(rightHand - range[j, 2])) {
              # print(paste("rightHand was", rightHand))
              rightHand = range[j, 1] - .00001
              # print(paste("rightHand is", rightHand))
            # }
            # else{
              # print(paste("rightHand was", rightHand))
              # rightHand = range[j, 2] + .00001
              # print(paste("rightHand is", rightHand))
            # }
          }
        }
      }
      #if we have established that left bound is within another boudnary, we dont want to keep expanding it
      if (boundhit2 != TRUE) {
        leftHand = middleOfPeak - (2 + .15 * counter) * peakSD
      }
      #if we extend to the left of 0, set bound to 0 and set boundhit2(left) to TRUE
      if (leftHand < 0) {
        leftHand = 0
        boundhit2 = TRUE
      }
      
      # takes a left hand bound that is within an old bound, and sets the new left bound to the old right bound
      if (length(range[, 1]) != 0) {
        for (j in 1:length(range[, 1])) {
          if (leftHand %between% range[j,]) {
            # print(
            #   paste(
            #     "leftHand, ",
            #     leftHand,
            #     "is inside range of previous boundary: ",
            #     as.numeric(range[j, ])
            #   )
            # )
            boundhit2 = TRUE
            #   leftHand <- closest value in vecofStarts OR vecOfstops
            # if (abs(leftHand - range[j, 1]) < abs(leftHand - range[j, 2])) {
              # print(paste("leftHand was", leftHand))
              # leftHand = range[j, 1] - .00001
              # print(paste("leftHand is", leftHand))
            # }
            # else{
              # print(paste("leftHand was", leftHand))
              leftHand = range[j, 2] + .00001
              # print(paste("leftHand is", leftHand))
            # }
          }
        }
      }
      
      #if we have hit both boundaries, return bounds.
      if (boundhit1 & boundhit2) {
        return(c(leftHand, rightHand))
      }
      if (signalToNoiseRatio > 1.1 * StartingSignalToNoiseRatio) {
        # print("signal/noise ratio significantly increased - break")
        # lastPeakGood <<- FALSE
        return(c(leftHand, rightHand))
      }
      if (signalToNoiseRatio < 3) {
        return(c(leftHand, rightHand))
      }
      
      
      #check if expanded region contains a sufficient number of points
      if (counter > 0) {
        # print(paste(leftHandOld, " is old left hand bound before recalculation"))
        leftHandOld <-
          middleOfPeak - (2 + .15 * (counter - 1)) * peakSD
        # print(paste(leftHandOld, " is old left hand bound"))
        rightHandOld <-
          middleOfPeak + (2 + .15 * (counter - 1)) * peakSD
        previousRegion <-
          mySubset$scan[mySubset$scan > leftHandOld &
                          mySubset$scan < rightHandOld]
        numPointsOld <- length(previousRegion)
      }
      else{
        numPointsOld = 0
      }
      signalRegion = mySubset[mySubset$scan > leftHand &
                                mySubset$scan < rightHand, ]
      numPointsNew <- length(signalRegion$scan)
      if (numPointsNew - numPointsOld < 10) {
        return(c(leftHandOld, rightHandOld))
      }
      signal = mean(signalRegion$intensity[signalRegion$intensity > quantile(signalRegion$intensity, probs =
                                                                               .33) &
                                             signalRegion$intensity < quantile(signalRegion$intensity, probs = .67)])
      #exclude anything within this interval
      intervalExclude = c(mySubset$scan > leftHand &
                            mySubset$scan < rightHand)
      ##print(intervalExclude)
      ##print("is intervalExclude ")
      #remove the interval
      mySubset1Exclude = mySubset[!intervalExclude, ]
      #resample after each iteration
      mySubset1Exclude$intensity[is.na(mySubset1Exclude$intensity)] = 0
      signalBackground <-
        mean(mySubset1Exclude$intensity[mySubset1Exclude$intensity < quantile(mySubset1Exclude$intensity, probs = .15)])
      if (is.nan(signalBackground)) {
        signalBackground = 1
      }
      if (is.na(signalBackground)) {
        signalBackground = 1
      }
      if (is.nan(signal)) {
        signal = 1
      }
      if (is.na(signal)) {
        signal = 1
      }
      #        print(paste(signalBackground, "is the background"))
      #        print(paste(signal, " is the signal"))
      #signal to noise ratio
      signalToNoiseRatio = signal / signalBackground
      if (is.na(signalToNoiseRatio)) {
        return(c(0, 3000))
      }
      #             print(paste(signalToNoiseRatio, " is the signal to noise ratio"))
      #iteration limit
      if (counter > 20) {
        return(c(leftHand, rightHand))
      }
      counter = counter + 1
    }
    return(c(leftHand, rightHand))
  }







#Datamatrix -> data file that contains all info- scans/intensity/file/etc
#mode -> multipeak or singlepeak
scanMaxesPerFile = function(DataMatrix, Mode) {
  mySubset = DataMatrix
  vecOFMaxes <- vector()
  for (j in 1:length(unique(mySubset$file)))
  {
    fileLookup =  unique(mySubset$file)[j]
    subsetPlot = mySubset[mySubset$file == fileLookup ,]
    myMax = subsetPlot[order(subsetPlot$intensity, decreasing = TRUE),][1:5,]$scan
    #get the scan with the max signal and store it in the vector
    #in order to determine the parameter of the drift
    vecOFMaxes = c(vecOFMaxes, myMax)
  }
  if (Mode == "singlePeak") {
    vecOFMaxes <- sort(vecOFMaxes)
    return(vecOFMaxes)
  }
  if (Mode == "multiPeak") {
    vecOFMaxes <- sort(vecOFMaxes)
    vecOFMaxes1 <- vecOFMaxes[vecOFMaxes < mean(vecOFMaxes)]
    vecOFMaxes2 <- vecOFMaxes[vecOFMaxes > mean(vecOFMaxes)]
    listOfVecMax <- list(vecOFMaxes1, vecOFMaxes2)
    return(listOfVecMax)
  }
}





#datamatrix -> contains all points that are not within another bound
#vecOfMaxes is from scanMaxesPerFile, it is used to get standard deviation and center of bounds
#vecOfStarts + Stops are used to make sure a new boundary is not overlapping an old one
boundaryScan = function(DataMatrix,
                        vecOfMaxes,
                        vecOfStarts,
                        vecOfStops) {
  vecOFMaxes <- vecOfMaxes
  vecOfStarts <- vecOfStarts
  vecOfStops <- vecOfStops
  
  #used for the plotting, basically saving the data before  truncation
  
  #this is going to be the data that we're going to remove signals from
  mySubset1ExcludedIntervals <- DataMatrix
  middleOfPeak = mean(vecOFMaxes)
  # print(paste(middleOfPeak, "is middle of peak"))
  #determine the range of the drift by the right and left hand extremes
  peakSD = sd(vecOFMaxes)
  if (is.nan(peakSD)) {
    return(c(0, 0))
  }
  if (is.na(peakSD)) {
    return(c(0, 0))
  }
  # print(paste(peakSD, "standard deviation of peak"))
  #helps prevent over expansion of peaks
  if (peakSD > 250) {
    peakSD = 250
  }
  #this helps boundary-expansion to capture all signal for very sharp peaks
  if (peakSD < 50) {
    peakSD = 60
  }
  rightHand = middleOfPeak + 2 * peakSD
  leftHand = middleOfPeak - 2 * peakSD
  
  #isolate region between left and right hand bound of calculated scan to calculate signal
  #another place where signalToNoise() function could be used
  
  # StartingSignalToNoiseRatio = signalToNoiseCalculator(mySubset1ExcludedIntervals, c(leftHand,rightHand))
  signalToNoiseRatio = signalToNoiseCalculator(mySubset1ExcludedIntervals, c(leftHand, rightHand))[1]
  
  #returns left and right hand expanded bounds given background region, signal region, initial signal/noise ratio, peak middle, and SD of peak
  #First variable- data set with identified "good" bounds removed
  #second variable- points outside of current bound to be evaluated (used for noise calculation)
  #(otherwise we would need to include left+righthand bound into function)
  
  bounds <-
    boundaryExpander(
      mySubset1ExcludedIntervals,
      signalToNoiseRatio,
      middleOfPeak,
      peakSD,
      vecOfStarts,
      vecOfStops
    )
  
  

  
  return(bounds)
}


#takes in current bounds to check and old bounds, returns unified bounds if there exists any overlap.
#arg1 = new bounds(left,right)
#arg2+3 = vector of left bounds, vector of right bounds
boundChecker = function(bounds, vecOfStarts, vecOfStops) {
  bounds <- bounds
  range <- cbind(vecOfStarts, vecOfStops)
  if (length(range[, 1]) != 0) {
    for (j in 1:length(range[, 1])) {
      if (bounds[2] %between% range[j,]) {
        # print(paste(
        #   "rightHand, ",
        #   bounds[2],
        #   "is inside range of previous boundary: ",
        #   range[j, ]
        # ))
        #rightHand <- left hand of old bound
        # print(paste("rightHand was", bounds[2]))
        #have to subtract a tiny amount as %inrange% does NOT like ties
        bounds[2] = range[j, 1] - .00001
        # print(paste("rightHand is", bounds[2]))
      }
      if (bounds[1] %between% range[j, ]) {
        # print(paste(
        #   "leftHand, ",
        #   bounds[1],
        #   "is inside range of previous boundary: ",
        #   range[j,]
        # ))
        # leftHand <- righthand of old bound
        # have to add a tiny amount as %inrange% does NOT like ties
        # print(paste("leftHand was", bounds[1]))
        bounds[1] = range[j, 2] + .00001
        # print(paste("leftHand is", bounds[1]))
      }
    }
  }
  return(bounds)
}


#mySubset - dataset to perform peak-finding on
#boundsLeftPeak - the bounds of the left peak to evaluate
#boundsRightPeak - bounds of the right peak to evaluate
#signalOriginal - maximum intensity value within the 3000 scan range for the mass being evaluated
#vecOfStartsGlobal - vector of all previous starts, used for preventing overlaps
#vecOfStops Global - vector of all previous stops, used for preventing overlaps
multiPeakEvaluator = function(mySubset,
                              boundsLeftPeak,
                              boundsRightPeak = c(1, 1),
                              signalOriginal,
                              vecOfStartsGlobal,
                              vecOfStopsGlobal) {
  #Initialize Variables
  vecOfStartsLocal <- vecOfStartsGlobal
  vecOfStopsLocal <- vecOfStopsGlobal
  #If either of these go to false- we will exit the overarching while loop on next iteration
  # peak1Good = TRUE
  # peak2Good = TRUE
  #Only plot bounds that are swapped to "true" that meet specific criteria
  plot1 = FALSE
  plot2 = FALSE
  #Used to determine if bounds should be added to VecOfStarts/Stops
  boundaryAdd = TRUE
  mySubset <- mySubset
  mySubsetCopy <- mySubset
  signalOriginal <- signalOriginal
  
  
  #remove signal from second peak, calculate signal/noise of first peak
  
  #signal of First Peak- could make signalToNoise() function
  mySubset <- mySubset[mySubset$scan > boundsLeftPeak[1] &
                         mySubset$scan < boundsLeftPeak[2], ]
  signal = mean(mySubset$intensity[mySubset$intensity > quantile(mySubset$intensity, probs =
                                                                   .33) &
                                     mySubset$intensity < quantile(mySubset$intensity, probs = .67)])
  signalMax = max(mySubset$intensity)
  if (signalMax < 0) {
    signalMax = 0
  }
  
  
  #Remove all signal within both peaks and calculate background and ratio
  #could make a signalToNoise() function... this is different from other cases due to the %inrange%... removing
  #other multipeak from consideration to get higher signal/noise
  range <- data.frame(matrix(0, ncol = 2, nrow = 3))
  colnames(range) <- c("rangeStart", "rangeEnd")
  range$rangeStart <- c(0.0, boundsLeftPeak[2], boundsRightPeak[2])
  range$rangeEnd <- c(boundsLeftPeak[1], boundsRightPeak[1], 3000)
  mySubsetBackground <-
    mySubsetCopy[as.double(mySubsetCopy$scan) %inrange% range, ]
  signalBackground <-
    mean(mySubsetBackground$intensity[mySubsetBackground$intensity < quantile(mySubsetBackground$intensity, probs = .15)])
  
  if (is.nan(signalBackground)) {
    signalBackground = 1
  }
  if (is.na(signalBackground)) {
    signalBackground = 1
  }
  if (is.nan(signal)) {
    signal = 1
  }
  if (is.na(signal)) {
    signal = 1
  }
  signalToNoiseRatio <- signal / signalBackground
  # print(
  #   paste(
  #     "For left peak, signalMax/signalOriginal",
  #     signalMax / signalOriginal,
  #     "and signalToNoise is ",
  #     signalToNoiseRatio
  #   )
  # )
  #check that no bounds overlap before adding the bound(this needs made into a function)
  
  #overlap checker --> looks to the left hand bound, if that bound is within a previously established bound
  #it is pushed to the right
  boundsLeftPeak <-
    boundChecker(boundsLeftPeak, vecOfStartsLocal, vecOfStopsLocal)
  
  
  
  
  #could be collapsed into addBoundIfGood() function {
  
  #only add the bounds if it meets our criteria
  
  #if it has a good signal to noise ratio OE the signal in the boun is > .90 of the original signal
  #lets us capturepeaks that spread out, because signal to noise will be low, but the signal will be high
  #i.e. lots of intensity, poor chromatography
  if (signalToNoiseRatio > 3 | signalMax > .9 * signalOriginal) {
    plot1 <- TRUE
    if (signalMax < .15 * signalOriginal) {
      plot1 <- FALSE
    }
    
    #if there is barely any signal, ignore (density filter)
    else{
      if (length(mySubset$intensity) < 20) {
        boundaryAdd = FALSE
      }
      
      
      #make sure that we don't have overlaps
      #if the left peak is to the left of one peak
      #but its right peak is to the right,keep the smaller window
      #per the discussion
      if (length(vecOfStarts != 0)) {
        for (k in 1:length(vecOfStarts)) {
          if (boundsLeftPeak[1] < vecOfStarts[k] &
              boundsLeftPeak[2] > vecOfStops[k]) {
            boundaryAdd = FALSE
            plot1 <- FALSE
          }
        }
        
        #if it passes the signal filter, density filter and doesn't encompass another window,
        #we're  good and we add it
        if (boundaryAdd) {
          vecOfStartsLocal <- c(vecOfStartsLocal, boundsLeftPeak[1])
          vecOfStopsLocal <- c(vecOfStopsLocal, boundsLeftPeak[2])
        }
      }
      #always add first boundary found
      else{
        vecOfStartsLocal <- c(vecOfStartsLocal, boundsLeftPeak[1])
        vecOfStopsLocal <- c(vecOfStopsLocal, boundsLeftPeak[2])
      }
    }
  }
  # }end addBoundIfGood() function
  # print(paste(signalToNoiseRatio, " = signal/noise for left peak"))
  #if ratio is small, the next peak is likely to be bad, same if the range of the first peak is very wide
  
  #if the signal to noise ratio is low or the signal is too spread out
  if (signalToNoiseRatio < 2 |
      (boundsLeftPeak[2] - boundsLeftPeak[1] > 800)) {
    plot1 = FALSE
  }
  
  
  
  #remove signal from first peak, calculate signal/noise of second
  #could make this a maxSignal() function
  mySubset <- mySubsetCopy
  mySubset <- mySubset[mySubset$scan > boundsRightPeak[1] &
                         mySubset$scan < boundsRightPeak[2], ]
  signal = mean(mySubset$intensity[mySubset$intensity > quantile(mySubset$intensity, probs =
                                                                   .33) &
                                     mySubset$intensity < quantile(mySubset$intensity, probs = .67)])
  signalMax = max(mySubset$intensity)
  if (signalMax < 0) {
    signalMax = 0
  }
  
  #Remove all signal within both peaks and calculate background and ratio
  #could make a signalToNoise() function to do this
  range <- data.frame(matrix(0, ncol = 2, nrow = 3))
  colnames(range) <- c("rangeStart", "rangeEnd")
  range$rangeStart <- c(0, boundsLeftPeak[2], boundsRightPeak[2])
  range$rangeEnd <- c(boundsLeftPeak[1], boundsRightPeak[1], 3000)
  mySubsetBackground <-
    mySubsetCopy[as.double(mySubsetCopy$scan) %inrange% range, ]
  signalBackground <-
    mean(mySubsetBackground$intensity[mySubsetBackground$intensity < quantile(mySubsetBackground$intensity, probs = .15)])
  
  
  
  boundaryAdd = TRUE
  if (is.nan(signalBackground)) {
    signalBackground = 1
  }
  if (is.na(signalBackground)) {
    signalBackground = 1
  }
  if (is.nan(signal)) {
    signal = 1
  }
  if (is.na(signal)) {
    signal = 1
  }
  signalToNoiseRatio <- signal / signalBackground
  # print(
  #   paste(
  #     "For right peak, signalMax/signalOriginal = ",
  #     signalMax / signalOriginal,
  #     " and signalToNoise is =",
  #     signalToNoiseRatio
  #   )
  # )
  #
  #our previous overlap test made sure the one peak didn't totally encompass another
  #this step, however, makes sure the one region of a peak doesn't overlap with part of
  #another peak at all.  We needed to make sure to stop overextending windows earlier
  
  #re-evaluate this with a fresh pair of eyes to potentially collapse the logic into
  #something more efficient
  
  boundsRightPeak <-
    boundChecker(boundsRightPeak, vecOfStartsLocal, vecOfStopsLocal)
  
  
  
  #could be collapsed into addBoundIfGood() function {
  
  #we're going to evaluate our potential peak on two criteria: look for some quick signal/noise checks that include max signal
  #and also encompass check where we make sure we aren't inside the bounds of another peak
  if (signalToNoiseRatio > 3 | signalMax > .9 * signalOriginal) {
    plot2 = TRUE
    if (signalMax < .15 * signalOriginal) {
      plot2 <- FALSE
      # peak2Good <- FALSE
    }
    else{
      # if (length(mySubset$intensity) < 20) {
      #   boundaryAdd = FALSE
      # }
      
      #this is the encompass check
      if (length(vecOfStarts != 0)) {
        for (k in 1:length(vecOfStarts)) {
          if (boundsRightPeak[1] < vecOfStarts[k] &
              boundsRightPeak[2] > vecOfStops[k]) {
            boundaryAdd = FALSE
            plot2 <- FALSE
          }
        }
        # if (boundaryAdd) {
        #   vecOfStartsLocal <-
        #     c(vecOfStartsLocal, boundsRightPeak[1])
        #   vecOfStopsLocal <- c(vecOfStopsLocal, boundsRightPeak[2])
        # }
      }
      #always add first bound that meets either 90% of max signal or signal/noise > 3
      else{
        vecOfStartsLocal <- c(vecOfStartsLocal, boundsRightPeak[1])
        vecOfStopsLocal <- c(vecOfStopsLocal, boundsRightPeak[2])
      }
    }
  }
  #} end addBoundIfGood() function
  # print(paste(signalToNoiseRatio, " = signal/noise for right peak"))
  #if signal/noise is low, next peak is likely to be bad
  if (signalToNoiseRatio < 2 |
      (boundsRightPeak[2] - boundsRightPeak[1] > 800)) {
    plot2 = FALSE
  }
  return(c(plot1,
           plot2))
  
  
  
}

#for generating figure pdfs
pdf(file = "C:/Users/Louis/OneDrive/Desktop/100Plots/TestDplyrFilter.pdf",
    width = 8,
    height = 5)

#make vector to store masses that have passed quality filters (>50% of samples in top chunk of intensity peak)
#top chunk = top 40-75% of max intensity
massesToKeep = vector()

#This is where the overarching function will be located- thinking about calling it Scanalyze
for (i in 1:100) {
  gc()
  print(i)
  #isolate mass for boundary evaluation.
  myMass <- oneHundredValidatedMasses[i]
  print(myMass)
  ppmAllowed = .000010
  rangeAllowed = ppmAllowed * myMass
  mySubset = BiggestTable %>% filter(mass < (myMass+rangeAllowed) & mass > (myMass-rangeAllowed))
  # [spikeTable$mass <  (myMass + rangeAllowed) &
  #                         spikeTable$mass > (myMass - rangeAllowed),]
  mySubset <- mySubset[order(mySubset$file), ]
  #make sure that the intensity column is numeric
  mySubset$intensity = as.numeric(mySubset$intensity)
  #replace any NA's with 0's
  mySubset$intensity[is.na(mySubset$intensity)] = 0
  #for now, we won't need the Mstype or the polarity info
  #if it is in the table, remove it!
  mySubset$polarity = NULL
  mySubset$massspectype = NULL
  #set colnames to make them consistent with other code
  colnames(mySubset) = c("mass", "intensity", "scan", "file")
  #currently unused, will store final results after running anovalign
  tableOfAnovalignAdjusted = colnames(mySubset)
  
  #remove the region of the retention time data referred to as "the void" - typically only low quality signals exist in this region
  mySubset = mySubset[mySubset$scan > 200, ]
  #copy of data to be used for graphing at the end of the boundary-finding. signals are iteratively removed from MySubset
  #In order to avoid recapturing the same signal again and again.
  mySubsetOG <- mySubset
  
  #exit condition for the while loop.
  #turns to false if signal/noise of new peak is low, or if the highest signal in the peak is 15% of the original signal
  lastPeakGood = TRUE
  
  #storage vectors for bounds identified
  vecOfStarts = vector()
  vecOfStops = vector()
  
  #highest signal present in original sample, used for calculating exit condition
  signalOriginal = max(mySubset$intensity)
  # print(paste(signalOriginal, " is the maximum intensity signal in sample ",i))
  
  #While loops that iteratively isolates peaks until an exit condition is hit (lastPeakGood)
  while (lastPeakGood) {
    #calculate signal background using the bottom 15 percent
    # signalBackground <-
      # mean(mySubset$intensity[mySubset$intensity < quantile(mySubset$intensity, probs = .15)])
    #write out how scanMaxesPerFile works - get the top 5 scans with maximum signal in each samples.
    #Then break this into upper and lower halves based on mean.
    
    
    #then, we;re going to examine the upper and lower halves and figure out if there are potentially two peaks
    
    #if there is just one peak, we'll take that good peak, get its bounds, remove the signal within the the bound
    #and continure searching for other peaks
    
    #if there are two good peaks, we'll go ahead and do the same for both.  and if only one peak is good, we won't remove
    #the signal in whichever half had the bad signal.  And hopefully, in the ext iteration, this "partial" signal
    #will get selected
    
    #Always start by assuming multipeak. Returns a list of 2 scan values that are predictive of peaks.
    # print(length(mySubset$intensity))
    vecOfMaxesList <- scanMaxesPerFile(mySubset, "multiPeak")
    
    #is the "top half" and "bottom half" sufficiently different?
    if (abs(mean(vecOfMaxesList[[1]]) - mean(vecOfMaxesList[[2]])) > 400) {
      # print(paste("sample ", i, " is a multipeak"))
      #if they are, call boundaryScan to find the bounds for peak 1 and peak 2
      
      #peak 1
      boundsLeftPeak <-
        boundaryScan(mySubset,
                     vecOfMaxesList[[1]],
                     vecOfStarts,
                     vecOfStops)
      # print(paste(
      #   "bounds of left peak are ",
      #   boundsLeftPeak[1],
      #   " and ",
      #   boundsLeftPeak[2]
      # ))
      
      #peak 2
      boundsRightPeak <-
        boundaryScan(mySubset,
                     # mySubsetOG,
                     vecOfMaxesList[[2]],
                     vecOfStarts,
                     vecOfStops)
      # print(paste(
      #   "bounds of right peak are ",
      #   boundsRightPeak[1],
      #   " and ",
      #   boundsRightPeak[2]
      # ))
      # 
      
      #returning info about whether or not we've got good peaks
      #if os, we'll remove their signal from the next iteration
      plots <-
        multiPeakEvaluator(
          mySubset,
          boundsLeftPeak,
          boundsRightPeak,
          signalOriginal,
          vecOfStarts,
          vecOfStops
        )
      
      
      #if the left peak is good, keep it and remove the points
      if (plots[1]) {
        vecOfStarts = c(vecOfStarts, boundsLeftPeak[1])
        vecOfStops = c(vecOfStops, boundsLeftPeak[2])
        intervalExclude = c(mySubset$scan > boundsLeftPeak[1] &
                              mySubset$scan < boundsLeftPeak[2])
        mySubset1Exclude = mySubset[!intervalExclude,]
        mySubset <- mySubset1Exclude
      }
      
      #if the right peak is good, keep it and remove the points
      if (plots[2]) {
        vecOfStarts = c(vecOfStarts, boundsRightPeak[1])
        vecOfStops = c(vecOfStops, boundsRightPeak[2])
        intervalExclude = c(mySubset$scan > boundsRightPeak[1] &
                              mySubset$scan < boundsRightPeak[2])
        mySubset1Exclude = mySubset[!intervalExclude,]
        mySubset <- mySubset1Exclude
      }
    
      #looks for two peaks and evaluates them, if either peak is low quality
      #the odds of subsequent peaks being high quality is very low
      #therefore, exit the while loop
      if (plots[1] == FALSE | plots[2] == FALSE) {
        lastPeakGood = FALSE
      }
      if (length(mySubset$intensity < 1000)) {
        lastPeakGood = FALSE
      }
      
      
      
      
    }
    #If top/bottom half of vecOfScans is not sufficiently different, assume one peak and recalculate most common scan value
    else {
      #make singlePeakEvaluator
      #There exist two conditions for keeping a peak for plotting
      #1. Is max intensity >15% of original max intensity.
      #2. Is signal/background of the peak > 3? Is the mean of the middle 2/3 of intensities in the peak 3x the value
      # of the mean of the bottom 15% of points outside the peak?
      plotCheck1 = TRUE
      plotCheck2 = TRUE
      # print(paste("sample ", i, " is a single peak"))
      #this is basically the same, except we assume that there's only one peak
      vecOFMaxes <- scanMaxesPerFile(mySubset, "singlePeak")
      bounds <-
        boundaryScan(mySubset,
                     vecOFMaxes,
                     vecOfStarts,
                     vecOfStops)
      # print(paste("bounds of single are ", bounds[1], " and ", bounds[2]))
      
      #Subset out region inside of identified bounds
      #Maximum signal in the bounds
    
      ratioAndMax <- signalToNoiseCalculator(mySubset, bounds)
      signalToNoiseRatio <- ratioAndMax[1]
      signalMax <- ratioAndMax[2]
      # print(signalMax)
      
      
      #Max function returns negative infinity if there are no points left to evaluate, set signal max to 0 to cause exit
      
      #note: background signal is currently calculated as the mean of the bottom 15% of remaining points once prior peaks
      #have been removed
  
      #should we plot this single peak? Check original signal,
      #signal/noise ratio, and if it overlaps an already existing bound
      #If any of these checks fail, dont plot it
      if (signalMax < .15 * signalOriginal) {
        lastPeakGood = FALSE
        plotCheck1 = FALSE
      }
      #Is mean signal in middle 2/3 of a peak 3 times higher than the mean signal of the bottom 15% of other points?
      if (signalToNoiseRatio < 3) {
        lastPeakGood = FALSE
        plotCheck2 = FALSE
      }
      if (plotCheck1 & plotCheck2) {
        vecOfStarts <- c(vecOfStarts, bounds[1])
        vecOfStops <- c(vecOfStops, bounds[2])
        mySubset <- mySubset[mySubset$scan > bounds[2] | mySubset$scan < bounds[1], ]
      }
    }
  }
  #Ordering the starts and stops to make the subsequent merging work
  vecOfStarts <- vecOfStarts[order(vecOfStarts, decreasing = F)]
  vecOfStops <- vecOfStops[order(vecOfStops, decreasing = F)]
  
  #initialize vectors for final bounds to be plotted
  vecOfStartsFinal <- vector()
  vecOfStopsFinal <- vector()
  
  
  #collapses those bounds that are within 200 scans of one another
  
  
  listOfFinalBounds <- list()
  #boundCollapser - if scans are w/in 200 of each other, collapse
  #vecOfStarts = left hand bounds
  #vecOfStops = right hand bounds
  #returns a list of two vectors with any starts/stops within 200 scans merged.
  listOfFinalBounds <- boundCollapser(vecOfStarts, vecOfStops)
  vecOfStartsFinal <- listOfFinalBounds[[1]]
  vecOfStopsFinal <- listOfFinalBounds[[2]]
  
  if(length(vecOfStartsFinal > 0)) {
  boundsToKeep <- highQualityPeaks(vecOfStartsFinal, vecOfStopsFinal, mySubsetOG)
  vecOfStartsFinal <- vecOfStartsFinal[boundsToKeep]
  vecOfStopsFinal <- vecOfStopsFinal[boundsToKeep]
  }
  # listOfVecOfStarts[[i]] <- vecOfStarts
  # listOfVecOfStops[[i]] <- vecOfStops
  
  #now we have our potential bounds, we're going to iterate through them 
  #and clean things up by sending the data between each set of bounds
  #to the anovalign fxn
  if (length(vecOfStartsFinal > 0)) {
    for (a in 1:length(vecOfStartsFinal)) {
      FinalZone <-
        mySubsetOG[mySubsetOG$scan < vecOfStopsFinal[a] &
                     mySubsetOG$scan > vecOfStartsFinal[a], ]
      adjustedZone <- generalizeWarping(FinalZone)
      tableOfAnovalignAdjusted <-
        rbind(tableOfAnovalignAdjusted, adjustedZone)
    }
    myMass <- oneHundredValidatedMasses[i]
    ppmAllowed = .000015
    rangeAllowed = ppmAllowed * myMass
    mySubsetAnovaAlign = tableOfAnovalignAdjusted[tableOfAnovalignAdjusted$mass <  (myMass + rangeAllowed) &
                                                    tableOfAnovalignAdjusted$mass > (myMass - rangeAllowed),]
  }
  
  

  
  colors <- c("red", "green", "blue", "black", "yellow")
  
  if(length(vecOfStartsFinal > 0)) {
  massesToKeep <- c(massesToKeep, myMass)
  plot(
    mySubsetOG$scan,
    mySubsetOG$intensity,
    xlab = "Scan",
    xlim = c(-10, 3100),
    ylab = "Intensity",
    col = as.factor(mySubsetOG$file),
    main = paste(
      "Left and Right Boundaries of Scans to Be Aligned for Mass",
      myMass
    )
  )
  for (z in 1:length(vecOfStarts)) {
    abline(v = vecOfStartsFinal[z], col = colors[z])
    abline(v = vecOfStopsFinal[z], col = colors[z])
  }
  }
  else{
    # print(paste("sample", i, " was thrown out due to lack of samples in middle 3rd of intensities"))
    plot(
      mySubsetOG$scan,
      mySubsetOG$intensity,
      xlab = "Scan",
      xlim = c(-10, 3100),
      ylab = "Intensity",
      col = as.factor(mySubsetOG$file),
      main = paste(
        "Left and Right Boundaries of Scans to Be Aligned for Mass",
        oneHundredValidatedMasses[i]
      )
    )
  }
  
  
  if (length(vecOfStartsFinal > 0)) {
    plot(
      mySubsetAnovaAlign$scan,
      mySubsetAnovaAlign$intensity,
      xlab = "Scan",
      xlim = c(-10, 3100),
      ylab = "Intensity",
      col = as.factor(mySubsetAnovaAlign$file),
      main = paste(
        "Left and Right Aligned Scans for Mass",
        myMass
      )
    )
    # legend("topleft", legend=as.numeric(unique(as.factor(mySubsetAnovaAlign$file))), pch = 1, col = as.numeric(unique(as.factor(mySubsetAnovaAlign$file))), cex = .4)
    for (y in 1:length(vecOfStarts)) {
      abline(v = vecOfStartsFinal[y], col = colors[y])
      abline(v = vecOfStopsFinal[y], col = colors[y])
    }
  }
  
  
}
sink()
dev.off()