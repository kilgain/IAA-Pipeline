#collapsing metabolite matcher final
####purpose - to find matches between the two sets of mass-features (IAA & Progenesis)

broadPublishedMassLarge = read.table(file = 'C:/Users/Louis/Downloads/HILp_NN_MN_to_readin.txt', header=TRUE, sep = "\t")
centroidTestBroad = broadPublishedMassLarge[,colnames(broadPublishedMassLarge) %in% c("m.z", "RT")]


#order data by mass so it centroids properly
centroidTestBroad = centroidTestBroad[order(centroidTestBroad$m.z),]

#These are the high confidence XCMS masses- those that are present in >300 samples. Subset data by these masses.
massesInterest = read.table(file = "oct18th_isopairs_first_pass_greater_300_with_isopairs.txt", sep = "\t", header = TRUE)
massesInterest  = massesInterest$Mass
centroidTestBroad = centroidTestBroad[centroidTestBroad$m.z %in% massesInterest, ]


broadACMassLarge = read.table("C:/Users/Louis/Downloads/broad_IAA_intensity_table.csv", header = T, sep = ',')
centroidTestBroadAC = broadACMassLarge[, colnames(broadACMassLarge) %in% c('compound','rt')]
colnames(centroidTestBroadAC) = c('m.z', 'RT')

#order data by mass so it centroids properly:
centroidTestBroadAC = centroidTestBroadAC[order(centroidTestBroadAC$m.z),]
#Need to standardize the retention time measurements (we had 3000 scans in 16 minutes, so this converts to minutes):
centroidTestBroadAC$RT = centroidTestBroadAC$RT / 3000 * 16


#used to tease apart features that share a mass, but have distinct retention times
centroidRTs = function(centroidByRTNow)
{
  
  #we're just going to take as input
  #the table of the RT's
  centroidByRTNow = centroidByRTNow
  
  allRTsInfo = vector()
  toBeCentroided = vector()
  
  #order, now, by the RT
  centroidByRTNow = centroidByRTNow[order(centroidByRTNow$RT),]
  
  
  #iterate out from one scan at a time
  #if any of the times are within 2 minutes, we're going to collapse them
  #into each other
  
  #set the old scan, first of all, to the first scan
  oldScan = centroidByRTNow$RT[1]
  
  #this will be the final vector of centroided scans
  centroidedScans = vector()
  
  #this will the temporary vector of scans that we wil centroid and add to the final
  #vector
  ScansToCentroid = vector()
  
  for(scan in 1:length(centroidByRTNow$RT))
  {
    
    #get the new scan
    newScan = centroidByRTNow$RT[scan]
    
    
    
    #is the difference between the new scan and old scan greater than our cutoff or not?
    gap = newScan - oldScan
    
    
    #look at the difference between this scan and the last scan.  Is it greater than 2 mins
    if(gap > 2)
    {
      
      #if we've moved on to a new scan but we haven't included the newest scan
      
      #print("we've got a gap to centroid now")
      
      #centroid the scans to be centroided, add them to the final vector 
      #AND reset the scans to be centroided
      
      #print(ScansToCentroid)
      centroidedScans = c(centroidedScans,mean(ScansToCentroid))
      
      #don't clear out yet, if we're on the last element
      ScansToCentroid = vector()
      
      
      #it's possible that the last scan will be outside the range of the ones added already
      #in this case, we'll want to have cleared up the ones that need to be centroided
      #BUT also add on the final scan to the list of RTs that will be returned
      
      
      #NOTE: it IS possible that we find ourselves in a situation where we don't get
      #any scans which achieve the gap to trigger centroided.  this, however, is fine and will be dealt
      #with in the final centroid step OUTSIDE this loop
      if(scan == length(centroidByRTNow$RT))
      {
        #add the element
        centroidedScans = c(centroidedScans, newScan )
      }
    }
    
    
    #look at the difference between this scan and the last scan.  Is it less than 2 mins?
    if(gap < 2)
    {
      
      #centroid the scans to be centroided, add them to the final vector 
      #AND reset the scans to be centroided
      ScansToCentroid = c(ScansToCentroid, newScan)
    }
    
    
    oldScan = newScan
  }
  
  #now, do we have any scans still to be centroided?
  #this is possible, because once we get to the last element
  #maybe we haven't had a gap in ANY of the scans
  #yet, which would have caused us to centroid the interval
  
  if(length(ScansToCentroid) > 0)
  {
    centroidedScans = c(centroidedScans,mean(ScansToCentroid))
  }
  
  #collapse all the RT's that are identical
  allRTsInfo = unique(centroidedScans)
  
  return(allRTsInfo)
}










#collapses mass values that are close to one another and have similar retention times
centroidFxn2D = function(centroidTestBroad)
{
  
  #get the m.z. and RT info
  tableCondensed <-as.data.frame(matrix(nrow=0,ncol=2))
  colnames(tableCondensed) = c("m.z", "RT")
  centroidTestBroad = centroidTestBroad
  
  #order the original table by mass here, now
  centroidTestBroad = centroidTestBroad[order(centroidTestBroad$m.z),]
  
  
  #have vectors to store the smoothed masses and RTs
  smoothedMassesAll = vector()
  smoothedRTsAll = vector()
  
  #read in the original vector of the masses
  myMassesOriginal = centroidTestBroad$m.z
  myMassesOriginal = unique(myMassesOriginal)
  
  
  oldMass = myMassesOriginal[1]
  centroidedMasses = vector()
  
  #have a vector to store the potential masses in a region of interest
  #don't take the mean until after you've passed through the the "gap"
  massesToCentroid = vector()
  vecsForROIRT = vector() 
  
  
  
  vecOfMultipleRows = vector()
  
  collapsedRTsEach = 0
  
  #keep track of the masses used for centroiding
  massesUsedForCentroiding = vector()
  
  
  allMassesProcessed = vector()
  totalProcessed = 0
  
  #iterate through the myVecTest of aasses
  
  for(i in 1:length(myMassesOriginal))
  {
    #print(i)
    #print(" is i ")
    mass = myMassesOriginal[i]
    newMass = mass
    
    #look within 5 ppm of the mass
    ppmAllowed = .000001 * 5* newMass
    
    #move forward and if the difference between the current mass and the old mass
    #is less than the threshold, we're going to centroid it eventually
    
    #is there a difference between this current mass and the old one?
    massGap = abs(newMass - oldMass)
    
    #CASE WHERE: mass is split!
    #& newMass != oldMass
    if(massGap <= ppmAllowed)
    {
      
      massesToCentroid = c(massesToCentroid, newMass)
    }
    
    #print(massGap)
    
    if(massGap > ppmAllowed)
    {
      
      
      if(length(massesToCentroid)  == 0)
      {
        massesToCentroid = oldMass
        
      }
      
      #print("we've encountered a split now !!")
      
      #centroid the vector of masses 
      centroidedMass = mean(massesToCentroid)
      
      
      
      #do we need to do the RT centroiding as well?
      toCentroidRT = centroidTestBroad[centroidTestBroad$m.z %in% massesToCentroid,]
      
      
      #resest the vector of massesToCentroid to be empty
      massesToCentroid = vector()
      
      #if we have multiple mass features from the same mass, we'll need to initialize those
      collapsedRTs = centroidRTs(toCentroidRT)
      
      
      
      
      #create the vector of masses for this table column now
      massColumn = rep(unique(centroidedMass), length(collapsedRTs))
      #print(massColumn)
      #create the table by binding the columns
      collapsedMassFeature = cbind(massColumn,collapsedRTs)
      collapsedMassFeature = as.data.frame(collapsedMassFeature)
      
      colnames(collapsedMassFeature) = colnames(toCentroidRT)
      
      #print(collapsedMassFeature)
      #print(" is the collapsedMassFeature")
      
      
      tableCondensed = rbind(tableCondensed, collapsedMassFeature)
      
      
      
    }
    
    
    #if we have the mass to centroid
    
    #if it's the last entry we won't have chance to reset, so include it now
    oldMass = newMass
  }
  
  
  #once we've run through everything in the main loop, it's possible that we'll have 
  #an array of masses to be cleared out
  if(length(massesToCentroid) > 0)
  {
    
    
    #print("we've encountered a split now !!")
    
    centroidedMass = mean(massesToCentroid)
    
    
    #do we need to do the RT centroiding as well?
    toCentroidRT = centroidTestBroad[centroidTestBroad$m.z %in% massesToCentroid,]
    
    #if we have multiple mass features from the same mass, we'll need to initialize those
    
    
    collapsedRTs = centroidRTs(toCentroidRT )
    
    
    
    
    #create the vector of masses for this table column now
    massColumn = rep(centroidedMass, length(collapsedRTs))
    #print(massColumn)
    #create the table by binding the columns
    collapsedMassFeature = cbind(massColumn,collapsedRTs)
    collapsedMassFeature = as.data.frame(collapsedMassFeature)
    
    colnames(collapsedMassFeature) = colnames(toCentroidRT)
    
    #print(collapsedMassFeature)
    #print("is thecollapsedMassFeature")
    #print("binding to tableCondensed!")
    
    tableCondensed = rbind(tableCondensed, collapsedMassFeature)
    
  }
  
  return(tableCondensed)
}








BroadCentroidedToCompare = centroidFxn2D(centroidTestBroad)

BroadCentroidedToCompareAC = centroidFxn2D(centroidTestBroadAC)









massFeatureMatcher = function(set1, set2) {
  matches = 0
  #3.5 minutes and +-15 ppm
  for(mass1 in 1:length(set1$m.z)) {
    if(mass1 %% 1000 == 0) {
      print(mass1)
    }
    massToMatch = set1$m.z[mass1]
    rtToMatch = set1$RT[mass1]
    massTolerance = .000015 * set1$m.z[mass1]
    potentialHits = set2[set2$m.z > massToMatch-massTolerance & set2$m.z < massToMatch+massTolerance & abs(set2$RT-rtToMatch) < 3.5, ]
    if(length(potentialHits$RT) == 0) {
      next
    }
    else{
      matches = matches + 1
    }
  }
  return(matches)
}






massFeatureMatcherMassOnly = function(set1, set2) {
  matches = 0
  #15 ppm
  for(mass1 in 1:length(set1$m.z)) {
    if(mass1 %% 1000 == 0) {
      print(mass1)
    }
    massToMatch = set1$m.z[mass1]
    rtToMatch = set1$RT[mass1]
    massTolerance = .000015 * set1$m.z[mass1]
    potentialHits = set2[set2$m.z > massToMatch-massTolerance & set2$m.z < massToMatch+massTolerance & abs(set2$RT-rtToMatch) < 100000, ]
    if(length(potentialHits$RT) == 0) {
      next
    }
    else{
      matches = matches + 1
    }
  }
  return(matches)
}



numMatches = massFeatureMatcher(BroadCentroidedToCompare,BroadCentroidedToCompareAC)



numMatchesMass = massFeatureMatcherMassOnly(BroadCentroidedToCompare,BroadCentroidedToCompareAC)






install.packages("eulerr")
library(eulerr)


#to make the venn diagram
fit <- euler(c("Autocredential" = 16308, "Progenesis" = 672, "Autocredential&Progenesis" = 4161))
plot(fit,quantities = list(cex = 1.5),fills = list(fill = c("red", "steelblue4"), alpha = 0.75),legend = list(labels = c("IsoLock + Autocredential", "Progenesis"), alpha = 1), main ="High Priority Mass-Intensity Signals at 15 ppm")
