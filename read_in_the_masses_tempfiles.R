library(data.table)
library(dplyr)
library(yaml)
library(stringr)
library(gtools)
#library(VPdtw)
library(minpack.lm)
#library(MassSpecWavelet)
library(splines)



#output File

#collapse this bracket to hide all functions
{
  
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

              boundhit1 = TRUE
              #   rightHand <- closest value in vecofStarts OR vecOfstops (NOTE: I added this BEFORE adding the same functionality
              # in multiPeakEvaluator(), where the right bound is always unified to an existing left hand bound to prevent overlap.)
              
              # if (abs(rightHand - range[j, 1]) < abs(rightHand - range[j, 2])) {
              # print(paste("rightHand was", rightHand))
              rightHand = range[j, 1] - .00001
          
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

              boundhit2 = TRUE

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
  
}



###allen code to sum up the regions

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
  
  print(head(metabolitePlotFull[order(metabolitePlotFull$intensity, decreasing = TRUE),]))
  
  maxReference = metabolitePlotFull[order(metabolitePlotFull$intensity, decreasing = TRUE),][3,]$file
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
      metabolitePlot = metabolitePlotFull[metabolitePlotFull$intensity < quantile(metabolitePlotFull$intensity, .25, na.rm = T), ]
      print(i)
    }
    if (i == 2)
    {
      #metabolitePlot = metabolitePlotFull[metabolitePlotFull$Intensity < quantile(metabolitePlotFull$Intensity, .75),]
      metabolitePlot = metabolitePlotFull[metabolitePlotFull$intensity < quantile(metabolitePlotFull$intensity, .75, na.rm =
                                                                                    T) &
                                            metabolitePlotFull$intensity > quantile(metabolitePlotFull$intensity, .25, na.rm = T)  , ]
      print(i)
    }
    if (i == 3)
    {
      metabolitePlot = metabolitePlotFull[metabolitePlotFull$intensity > quantile(metabolitePlotFull$intensity, .75, na.rm = T), ]
      print(i)
    }
    #for the tertile, pull out the signal belonging to the reference sample to be the reference sample
    metabolitePlotSlice1 =  metabolitePlot[metabolitePlot$file == maxReference, ]
    #initialize the column that will contain the adjusted scan info
    metabolitePlotSlice1$newRT = metabolitePlotSlice1$scan
    #initialize the scan adjusted table, which will contain the original reference data
    #as well as the scan adjust info for each file
    both = metabolitePlotSlice1
    #take out the top of the reference peak
    metabolitePlotSlice1 = metabolitePlotSlice1[metabolitePlotSlice1$intensity > quantile(metabolitePlotSlice1$intensity, .75, na.rm = T), ]
    print(paste(maxReference, "is MaxReference"))
    for (i in 1:length(unique(metabolitePlot$file)))
    {
      #print(i)
      if (i != maxReference)
      {
        metabolitePlotSlice2 = metabolitePlot[metabolitePlot$file == unique(metabolitePlot$file)[i], ]
        metabolitePlotSlice2BeforeFilter = metabolitePlotSlice2
        #take the top of the peak to compare to the reference sample as well.  Thus, we are taking the top of each 'zone' in each file
        metabolitePlotSlice2 = metabolitePlotSlice2[metabolitePlotSlice2$intensity > quantile(metabolitePlotSlice2$intensity, .75, na.rm = T), ]
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
      print(nrow(both ))
      print("is nrow of both")
      bothFull = rbind(bothFull, both)
      print("binding both")
      plot(both$scan, both$intensity, main = "after adjusting")
      plot(both$scan, both$intensity, main = "before adjusting")
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
  #read in the full temp file into a single table
  megaTableII = as.data.table(table)
  megaTableII$file <- str_remove(megaTableII$file, ".rawoutput.txt")
  print("made megatable")
  print(head(megaTableII)) 


  colnames(megaTableII) =  c("file", "compound","intensity", "scan")
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
    allInfoTableFile = matrix(nrow= 0, ncol = length(c("compound","rt", numberOfSampleFiles)))
    allInfoTableFile = as.data.frame(allInfoTableFile)
    
    #note that this could be challenging if we don't have 
    #the files that we need in the subset
    colnames(allInfoTableFile) =c("compound","rt", as.character(numberOfSampleFiles))
    
    
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
        colnames(toSmoothAndSum) = c("compound", "intensity", "scan")
        #smoothedAddUp = toSmoothAndSum$intensity
        
        
        #myAdjustedTable = as.data.table(myAdjustedTable)
        subsetToAdd = toSmoothAndSum[toSmoothAndSum$scan > myWindowStart & toSmoothAndSum$scan < myWindowStop, ]$intensity
        subsetToAdd = subsetToAdd[subsetToAdd != 0]
        smoothedAddUp = subsetToAdd
        
        #quantileTest = quantile(subsetToAdd, c(.25, .85), na.rm = TRUE)
        
        
        if(length(subsetToAdd) > 9)
        {
          
          #print(quantileTest)
          #print(" is quantileTest")
          #print("greater than 10 in length")
          
          #Old Quantile Filter = [smoothedAddUp > quantileTest[1] & smoothedAddUp < quantileTest[2]]
          
          intensity = sum(smoothedAddUp)
        }
        
        #if we have almost all zeroes, taking the quantile is not going to be useful
        if(length(subsetToAdd) < 10)
        {
          
          intensity = sum(subsetToAdd)
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
    }
    
    #if(sigRatioFirst > 3)
    #{
    
    allInfoTableFile[is.na(allInfoTableFile)] = 0
    #allInfoTableAllFiles = rbind(allInfoTableAllFiles, allInfoTableFile)
    #}	
  }
  
  
  #make sure to turn all NA's to 0
  allInfoTableFile[is.na(allInfoTableFile)] = 0
  print( allInfoTableFile)
  
  return(allInfoTableFile)
}


#we're going to have the parent directory to keep track of everyhting

config = read_yaml("/home/lconnelly/Metabolomics/Config.yaml")




args = commandArgs(trailingOnly=TRUE)
vecOfMasses = basename(args)




#have an output directory to put everything
#currentDir = "/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/setaria_plate9_chromatogram_files_final"
#currentDir = "/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/setaria_plate_9_standards"
#currentDir = "/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/all_broad_autocredential"
#currentDir = "/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/HILIC_12_plates_setaria"
#currentDir = "/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/HILIC_12_plates_setaria_standards"
#currentDir = "/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/temp_files_broad_citrulline"
#currentDir = "/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/june16th_bins_filtered_at_50_counts"
#currentDir = "/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/broad_isotopes_multiple"

#currentDir = "/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/DOE_setaria_and_HILIC"
#currentDir = "/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/standards_DOE_HILIC_louis_matches"
#currentDir = "/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step3_convert_get_the_temp_files/standards_louis_more_comprehensive"

currentDir = args

#allDataDirectory = paste(getwd(),"/",currentDir, sep = "")
#makeDirCommandDatFile =  paste("mkdir", allDataDirectory, sep = " ")

#system(makeDirCommandDatFile)


  #options(digits = 20)
  myMassTest = as.numeric(vecOfMasses)
  
  #old adjustment= - (.00000367)*myMassTest
  myMassTest = myMassTest
  myMassTest = round(myMassTest, digits = 5)
  print(myMassTest)
  myMassTestQuery = myMassTest
 

  #set this to 0 if not working on isotopes  
  massSubDirectory = currentDir
  print(massSubDirectory)
  print(" is where we're going to read from")
  #print(massSubDirectory)
  #setwd(massSubDirectory)

  file_list <- list.files(path=massSubDirectory, full.names = TRUE)
  file_list = file_list[sapply(file_list, file.size) > 0]
#  print(file_list)
#  print(" is file_list")

  fullFileName = paste(currentDir, "/",  vecOfMasses,"merged.txt", sep = "")
  mergedExclude = fullFileName


if(!(fullFileName %in% file_list))
{


  allFiles = read.table(file_list[1], header = TRUE, sep = ",")
  #print(head(allFiles))
  #colnames(allFiles) = c("row","File","Mass","Instensity","Scan", "Mode","MsType")

  print(allFiles)
  print(" is allFiles")


  #print(mergedExclude)
  file_list = file_list[-1]
  if(length(file_list) > 0)
  {
  	for (j in 1:length(file_list))
  	{
   		print(j)
   		print("is j")
   		if(file_list[j] != mergedExclude)
   		{   
        
    		#oneFile = tryCatch(read.table(file_list[j], header = TRUE, sep = '|'), error=function(e) NULL)
    			oneFile = read.table(file_list[j], header = TRUE, sep = ",")
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
  #allFiles = allFiles %>% filter(newMass > myMassTest-(myMassTest*.0000025) & newMass < myMassTest+(myMassTest*.0000025))
  

#  write.table("helloAllFilesCreatedNew", "/home/lconnelly/Metabolomics/test00.txt")

  write.table(allFiles, file=fullFileName, row.names = FALSE, sep = "\t")

}else {
#  write.table("helloAllFilesCreatedPreviouslyBeforeRead-in", "/home/lconnelly/Metabolomics/test01.txt")
 
 allFiles = read.table(fullFileName, header = TRUE, sep = "\t")

#  write.table("helloAllFilesCreatedPreviouslyAfterRead-in", "/home/lconnelly/Metabolomics/test02.txt")
}	


write.table("hello", file = "/home/lconnelly/Metabolomics/testThruReadIn.txt")

# write.table("hello", "/home/lconnelly/Metabolomics/test1.txt")



  #initialize QC info for high intensity signals

  #colnames(isoInfo) = c("Mass", "mainPeakHigh")

  #mainPeakHigh = length(unique(allFiles[allFiles$Intensity > quantile(allFiles$Intensity, p = .75),]$File))


  #isoInfoMass = t(data.frame(c(myMassTest, mainPeakHigh)))
  #print(isoInfoMass)
  #print( colnames(isoInfo))
  #colnames(isoInfoMass) = colnames(isoInfo)

  #isoInfo = rbind(isoInfo, isoInfoMass)

  #output file for the table listing the number of samples in high intensity region
  #QCInfo = paste("/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step4_combine_temp_files/QC_info_each_mass","/", myMassTest,"july21st_DOE.txt", sep = "")
  

  #write.table(allFiles,file = fullFileName, sep = "\t", row.names = FALSE) 




  
  #fileToReadIn = read.table(file = "/home/ahubbard/complete_auto_credential_pipeline_start_to_finish/step4_combine_temp_files/michael_wei_lipids_mass_and_RT.txt", sep = "\t", header= TRUE)
  #fileToReadIn$RT = fileToReadIn$RT *60 * 1.7

 
#  name = paste(config$metabolomicsRoot,"/outputPdfs/", myMassTest, "_plotinfo.pdf", sep = "")
#  myMass = myMassTest 
  #mySubsetMetaInfoSubset = fileToReadIn[fileToReadIn$Mass == myMass, ]

  #print("is the match")
  #print(name)
  #print(" is the name of the pdf")

#  pdf(name)
#  plot(allFiles$newMass, allFiles$Intensity, main = "10PPM Window")
#  abline(v = myMass, col = "green")

  #plot(allFiles$scan, allFiles$Intensity)
  #abline(v = myMass, col = "green")

  #ignore RT
  #plot mass intensity
  #next, does it follow a pattern of intensity that is roughly consistent across samples regardless of injection order? 
  #when I look at the top 10% of the merged peak
  #merge everything
  #take top 10% - proxy for retention time
  #pull out intensity from each file and sum them
  #order by injection
   
 




 #myRTRight = mySubsetMetaInfoSubset$RT + 200
  #myRTLeft = mySubsetMetaInfoSubset$RT - 200

 # plot(allFiles$scan, allFiles$Intensity)
  #abline(v = myRTRight, col = "green")
  #abline(v = myRTLeft, col = "green")

#write.table("beforeConsistentInjections", file = paste(config$metabolomicsRoot, "/test2.txt", sep = ""))



###consistentWithInjections = list.files(currentDir)
###x = nchar(basename(currentDir)) + 1
###consistentWithInjections = substring(consistentWithInjections, first = x)
###consistentWithInjections = str_remove(consistentWithInjections, ".rawoutput.txt.hdf5newColumn_plate.hdf5.txt")



###injectionOrders <- list.files(config$injectionLists, full.names = T)
###injectionOrders = mixedsort(injectionOrders)

#if(myMassTest == 182.1234) {
#	write.table(injectionOrders, file = paste(config$metabolomicsRoot, "/injectionOrdersCitrulline.txt", sep = ""))
#}

###allInjections = vector()
###labelBatches = vector()
###for(f in 1:length(injectionOrders)) {
###  injectionOrder = read.table(injectionOrders[f], sep = ",", header = TRUE, stringsAsFactors = FALSE)
###  myOrder = injectionOrder$File.Name
###  myOrder = sub(" ", "_", myOrder)
###  for(a in 1:length(myOrder)) {
###    if(substr(myOrder[a], start = 1, stop = 4) != "HILI") {
###      myOrder[a] = paste("HILIC_", myOrder[a], sep="")
###    }
###  }
###  myOrder = myOrder[myOrder %in% intersect(myOrder, consistentWithInjections)]
###  labelBatches = c(labelBatches, rep(f, times=length(myOrder)))
###  allInjections = c(allInjections, myOrder)
###}

###consistentWithInjections = consistentWithInjections[match(allInjections,consistentWithInjections)]

#write.table(allInjections, file = paste(config$metabolomicsRoot, "/allInjections.txt", sep = ""))




#pipeline for bounds goes here: outputs table where colnames are files and first row is intensity (first three columns are compound(mass), rt, x
#Before running, need to initialize output table:

##This needs  to be the list of all files (should add into yaml config file)

fileList = read.table(config$fileList, header = T)$x

fileList = str_remove(fileList, ".rawoutput.txt")

nameOfAllFiles = fileList
numberOfSampleFiles = nameOfAllFiles

numberOfSampleFiles = unique(numberOfSampleFiles)

allInfoTableAllFiles = matrix(nrow= 0, ncol = length(c("compound","rt", numberOfSampleFiles)))
allInfoTableAllFiles = as.data.frame(allInfoTableAllFiles)

#note that this could be challenging if we don't have 
#the files that we need in the subset
colnames(allInfoTableAllFiles) =c("compound","rt", as.character(numberOfSampleFiles))



#Pipeline for mass features:
for (i in 1:1)
{
  myMass = myMassTest
  # myMass = 182.1230
  print(myMass)
  ppmAllowed = .0000025
  rangeAllowed = ppmAllowed * myMass
    
  myMassTest = myMass

 
  #allDataDirectory = currentDir
  #massSubDirectory = paste(allDataDirectory, "/", myMassTest, sep = "")
  #fullFileName = paste(allDataDirectory, "/", myMassTest,"/", myMassTest,"merged.txt", sep = "")
    

    
  
  mySubset<- allFiles
  
  mySubset <- mySubset %>% filter(newMass > myMassTest-(.0000025)*myMassTest & newMass < myMassTest+(.0000025)*myMassTest)
  
 
  #make sure that the intensity column is numeric
  mySubset$Intensity = as.numeric(mySubset$Intensity)
  #replace any NA's with 0's
  mySubset$Intensity[is.na(mySubset$Intensity)] = 0
  #for now, we won't need the Mstype or the polarity info
  #if it is in the table, remove it!
  mySubset$polarity = NULL
  mySubset$massspectype = NULL
  mySubset$newMass = NULL
  #set colnames to make them consistent with other code
  #colnames(mySubset) = c("mass", "intensity", "scan", "file")
  mySubset$Mstype = NULL
  colnames(mySubset) = c("file", "mass", "intensity", "scan")
  #currently unused, will store final results after running anovalign
  mySubset = mySubset[order(mySubset$file),]
  tableOfAnovalignAdjusted = colnames(mySubset)
  
  #remove the region of the retention time data referred to as "the void" - typically only low quality signals exist in this region
  mySubset = mySubset[mySubset$scan > 200,]
  #copy of data to be used for graphing at the end of the boundary-finding. signals are iteratively removed from MySubset
  #In order to avoid recapturing the same signal again and again.
  mySubsetOG = mySubset
  
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
    print("iterating!")
    
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

      #peak 2
      boundsRightPeak <-
        boundaryScan(mySubset,
                     # mySubsetOG,
                     vecOfMaxesList[[2]],
                     vecOfStarts,
                     vecOfStops)

      
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
        mySubset1Exclude = mySubset[!intervalExclude, ]
        mySubset <- mySubset1Exclude
      }
      
      #if the right peak is good, keep it and remove the points
      if (plots[2]) {
        vecOfStarts = c(vecOfStarts, boundsRightPeak[1])
        vecOfStops = c(vecOfStops, boundsRightPeak[2])
        intervalExclude = c(mySubset$scan > boundsRightPeak[1] &
                              mySubset$scan < boundsRightPeak[2])
        mySubset1Exclude = mySubset[!intervalExclude, ]
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
      print("going into else statement")
      
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
        mySubset <-
          mySubset[mySubset$scan > bounds[2] | mySubset$scan < bounds[1],]
      }
    }
  }
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
  

  table1Original = mySubsetOG
  
  
  if(!(is.numeric(vecOfStartsFinal))){
    vecOfStartsFinal = 0
  }
  if(!(is.numeric(vecOfStopsFinal))){
    vecOfStopsFinal = 3000
  }
  
  print(vecOfStartsFinal)
  print(vecOfStopsFinal)
  
  
  
  #table = tableOriginal
  #now, go through each of the "zones" and determine if they need to be merged or not
  for(i in 1:length(vecOfStartsFinal))
  {
    
    start1 = vecOfStartsFinal[i]
    stop1 =  vecOfStopsFinal[i]
    
    segment1 = table1Original[table1Original$scan > start1 & table1Original$scan < stop1,]
    #segment1 = segment1[segment1$Intensity > quantile(segment1$Intensity, p = .9),]
    segment1 = generalizeWarping(segment1)
    
    print(head(segment1))
    #run anovAlign on the segment now!
    
    #segment1$Mstype = NULL
    #segment1$polarity= NULL
    tableForRegion = callAllOtherFxns(segment1)
    allInfoTableAllFiles = rbind(allInfoTableAllFiles,tableForRegion)
  } 

###three stars = old iteration

name = paste(config$metabolomicsRoot,"/outputPdfs5C/", vecOfMasses, "_plotinfo.pdf", sep = "")

###myMass = myMassTest

pdf(name)
###print(head(segment1))
plot(segment1$scan, segment1$intensity, main = "ANOVALIGN")  

}


write.csv(allInfoTableAllFiles, file = paste("/home/lconnelly/Metabolomics/vecOfIntensitiesForEachMass5C/", vecOfMasses, "-VecOfIntensities.csv", sep = ""), row.names=F)



#Now we have the allInfoTableAllFiles
#cols = compound, rt, (x?), [files in order of injection]
#row1 = sum of intensity within RT bounds + 5ppm



plot(allFiles$newMass, allFiles$Intensity, main = "10PPM Window")
abline(v = myMass, col = "green")

plot(allFiles$scan, allFiles$Intensity)
abline(v = vecOfStartsFinal[1], col = "green")
abline(v = vecOfStopsFinal[1], col = "green")

dev.off()


###start = 6
###for(p in 1:5)
###{
###	if(colnames(allInfoTableAllFiles)[p] != "compound" & colnames(allInfoTableAllFiles)[p] != "rt" & colnames(allInfoTableAllFiles)[p] != "x"){
###		if(start > p){
###		start = p
###		}
###	} 
###}

###vecOfIntensities <- allInfoTableAllFiles[1, start:length(colnames(allInfoTableAllFiles))]
###vecOfIntensities <- as.numeric(vecOfIntensities[1, ])
###vecOfIntensities <- as.data.frame(vecOfIntensities)
###vecOfIntensities$names <- colnames(allInfoTableAllFiles)[start:length(colnames(allInfoTableAllFiles))]
#vecOfIntensities$standard = vecOfMasses
#write.csv(vecOfIntensities, file = paste(config$metabolomicsRoot, "/vecOfIntensities/",vecOfMasses.txt, sep = ""))

#vecOfIntensities$vecOfIntensities <- log10(vecOfIntensities$vecOfIntensities + 1)




#write.csv(vecOfIntensities, file = "/home/lconnelly/Metabolomics/newVecOfIntensitiesCitrullineEndogenou.txt", row.names=F)

###plot(as.numeric(vecOfIntensities$vecOfIntensities), main = "STANDARD BEING TESTED")




###vecOfIntensitiesCitrulline <- read.table(paste(config$metabolomicsRoot, "/newVecOfIntensitiesCitrulline.txt", sep=""), header = T, sep = ",")


#vecOfIntensitiesCitrulline$vecOfIntensities <- log10(vecOfIntensitiesCitrulline$vecOfIntensities + 1)


###vecOfIntensities <- vecOfIntensities[vecOfIntensitiesCitrulline$names %in% vecOfIntensities$names,]
###vecOfIntensities <- vecOfIntensities[vecOfIntensities$names %in% vecOfIntensitiesCitrulline$names,]



###vecOfIntensitiesCitrulline <- vecOfIntensitiesCitrulline[vecOfIntensitiesCitrulline$names %in% vecOfIntensities$names,]
###vecOfIntensitiesCitrulline <- vecOfIntensitiesCitrulline[vecOfIntensities$names %in% vecOfIntensitiesCitrulline$names,]


###print(length(vecOfIntensities$vecOfIntensities))
###print(vecOfIntensities$names[800])

###print(length(vecOfIntensitiesCitrulline$vecOfIntensities))
###print(vecOfIntensitiesCitrulline$names[800])



#print(head(vecOfIntensitiesCitrulline))

###plot(as.numeric(vecOfIntensitiesCitrulline$vecOfIntensities), main = "CITRULLINE CONTROL")

#allFiles = allFiles[allFiles$Intensity > quantile(allFiles$Intensity, probs = .90),]

#allFiles = allFiles %>% filter(newMass > myMassTest-(myMassTest*.0000025) & newMass < myMassTest+(myMassTest*.0000025))

###print(cor(as.numeric(vecOfIntensities$vecOfIntensities), as.numeric(vecOfIntensitiesCitrulline$vecOfIntensities)))



###if(cor(vecOfIntensities$vecOfIntensities, vecOfIntensitiesCitrulline$vecOfIntensities) > .6) {
###	fileOutput = paste(vecOfMasses,"anovalign_july30_all_standards_comprehensive_PASS.txt", sep = "-" )
###	write.table(allInfoTableAllFiles, file = paste("/home/lconnelly/Metabolomics/outputCorrelations/",fileOutput,sep=""))
###}else{
###	fileOutput = paste(vecOfMasses,"anovalign_july30_all_standards_comprehensive_FAIL.txt", sep = "-" )
###	write.table(allInfoTableAllFiles, file = paste("/home/lconnelly/Metabolomics/outputCorrelations/",fileOutput,sep=""))
###
###}



#vecOfIntensities= vector()
#for(i in 1:length(consistentWithInjections)) {
#        file = allFiles[grepl(consistentWithInjections[i], allFiles$File),]
#	intensitySum = sum(file$Intensity)
#	vecOfIntensities = c(vecOfIntensities, intensitySum)
#}


#if(myMassTest == 182.1234) {
#	write.table(vecOfIntensities, file = paste(config$metabolomicsRoot, "/vecOfIntensitiesCitrulline.txt", sep = ""))
#}


###dev.off()

#if(myMassTest == 182.1234) {
#	write.table(consistentWithInjections, file = paste(config$metabolomicsRoot, "/indexForInjectionPlotCitrulline.txt", sep = ""))
#}


  #write.table(isoInfo,file = QCInfo, sep = "\t", row.names = FALSE)









