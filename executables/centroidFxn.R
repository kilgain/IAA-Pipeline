
args = commandArgs(trailingOnly = T)




centroidFxn= function(centroidTestBroad)
{
  centroidTestBroad = centroidTestBroad
  #global remade table
  myVecTest = vector()
  myVecTest = centroidTestBroad[order(centroidTestBroad)]
  #order all mass signals by mass
  myVecTest <- myVecTest[order(myVecTest)]
  oldMass = 0
  centroidedMasses = vector()
  #have a vector to store the potential masses in a region of interest
  #don't take the mean until after you've passed through the the "gap"
  #which is the space between signals which would indicate a new mass's gaussian 
  #profile
  #create a vector to hold all masses in a potential region of interest (ROI)
  vecsForROI = vector()
  for(i in 1:length(myVecTest))
  {
    mass = myVecTest[i]
    newMass = mass
    #determine allowable gap between signals as being 3 parts per million
    ppmAllowed = .000001 * 3 * newMass
    #if we are below the threshold for gaps between signals
    #store the mass 
    if(abs(newMass - oldMass) <= ppmAllowed)
    {
      if(length(vecsForROI) == 0)
      {
        #make sure that we inlcude the previous mass in the ROI vector to be averaged
        vecsForROI = c(oldMass)
        #make sure that we don't inclue the old mass ins the centroidedMasses as well!
        centroidedMasses = centroidedMasses[-which(centroidedMasses == oldMass)]
      }
      vecsForROI = c(vecsForROI, mass)
    }
    #we've got a hit that's past the threshold
    #take the average of all the masses in the ROI vector and append 
    #it to a vector of all the fully centroided masses
    if(abs(newMass - oldMass) > ppmAllowed || (i == 1))
    {
      #if the ROI vector has a bunch of masses ( > 1 ) to be centroided
      #take the average of everything in the interval
      if(length(vecsForROI) > 0)
      {
        #take the average of the ROI vector (vecsForROI) and append to the 
        #vector of centroided masses
        centroidedMasses = c(centroidedMasses,mean(vecsForROI))
        #initialize with the first scan
        vecsForROI = vector()
      }
      #if we have a standalone mass, it won't need to be centroided
      #but it won't be a part of any other masses, so we'll just need to add 
      #it to the vector of centroided masses
      if(length(vecsForROI) == 0)
      {
        centroidedMasses = c(centroidedMasses,newMass)
      }
    }
    #if it's the last entry we won't have chance to reset, so include it now
    oldMass = newMass
  }
  return(centroidedMasses)
}



outputPath = args[2]

masses = scan(args[1])

masses = centroidFxn(masses)

#print(length(masses))

write.table(masses, file=paste(outputPath, '/', basename(args[1]),  "centroidMasses.txt", sep=""), row.names=F)

