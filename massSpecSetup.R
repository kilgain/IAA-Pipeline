 

spikeCounts <- list.files("C:/Users/Louis/Downloads/massSpec/spike/", full.names = T)
tableList = list()
for(i in 1:length(spikeCounts))
{
  #clear memory
  gc()
  file = spikeCounts[i]
  table<- read.delim(file, sep = "\t", header = F)
  #removes duplicate masses that contain N/A's
  table<- table[complete.cases(table), ]
  colnames(table) = c("mass", "intensity", "scan", "polarity", "massspectype")
  table$'NA' = NULL
  tableList[[i]] = table
  tableList[[i]]$mass <- round(tableList[[i]]$mass, digits=5)
}

for(i in 1:(length(tableList)-1)) {
  gc()
  if(i == 1) {
    matrix1 <- tableList[[i]]
  }
  matrix2 <- tableList[[i+1]]
  spikeConfirmedMasses <- intersect(matrix1$mass, matrix2$mass)
  #matrix1 will also be used to store results
  matrix1 <- tableList[[1]][tableList[[1]]$mass %in% spikeConfirmedMasses, ]
}
#transfer result to new dataframe as matrix1 gets reused
spikeConfirmedResults <- matrix1[matrix1$mass %in% spikeConfirmedMasses, ]
spikeConfirmedResultsSample2 <- matrix2[matrix2$mass %in% spikeConfirmedMasses, ]
# head(spikeConfirmedResults)
# length(tableList[[1]]$mass)
# length(tableList[[2]]$mass)
# length(spikeConfirmedResults$mass)
# length(spikeConfirmedResultsSample2$mass)
# length(spikeConfirmedMasses)
#  

 

nospikeCounts <- list.files("C:/Users/Louis/Downloads/massSpec/nospike/", full.names = T)
for(i in 1:length(nospikeCounts))
{
  #clear memory
  gc()
  file = nospikeCounts[i]
  table<- read.delim(file, sep = "\t", header = F)
  #removes duplicate masses that contain N/A's
  table<- table[complete.cases(table), ]
  colnames(table) = c("mass", "intensity", "scan", "polarity", "massspectype")
  table$'NA' = NULL
  tableList[[i]] = table
  tableList[[i]]$mass <- round(tableList[[i]]$mass, digits=5)
}

for(i in 1:(length(tableList)-1)) {
  gc()
  if(i == 1) {
    matrix1 <- tableList[[i]]
  }
  matrix2 <- tableList[[i+1]]
  nospikeConfirmedMasses <- intersect(matrix1$mass, matrix2$mass)
  #matrix1 will also be used to store results
  matrix1 <- tableList[[1]][tableList[[1]]$mass %in% nospikeConfirmedMasses, ]
}
#transfer result to new dataframe as matrix1 gets reused
nospikeConfirmedResults <- matrix1[matrix1$mass %in% nospikeConfirmedMasses, ]
nospikeConfirmedResultsSample2 <- matrix2[matrix2$mass %in% nospikeConfirmedMasses, ]
# head(nospikeConfirmedResults)
# length(tableList[[1]]$mass)
# length(tableList[[2]]$mass)
# length(nospikeConfirmedResults$mass)
# length(nospikeConfirmedResultsSample2$mass)
# length(nospikeConfirmedMasses)
 



 
standardConfirmedMasses <- setdiff(spikeConfirmedResults$mass, nospikeConfirmedResults$mass)
standardConfirmed <- spikeConfirmedResults[spikeConfirmedResults$mass %in% standardConfirmedMasses, ]
sorghumConfirmedMasses <- intersect(spikeConfirmedResults$mass, nospikeConfirmedResults$mass)
sorghumConfirmed <- nospikeConfirmedResults[spikeConfirmedResults$mass %in% sorghumConfirmedMasses, ]
sorghumConfirmed2 <- nospikeConfirmedResults[nospikeConfirmedResults$mass %in% sorghumConfirmedMasses, ]
 

 
# #number of elements in noise reduced spike
# length(spikeConfirmedResults$mass)
# 
# #number of elements in noise reduced nonspike
# length(nospikeConfirmedResults$mass)
# 
# #number of unique masses different between spike and nonspike... my guess is that there exists a good chunk of noise in this data still. Hence the further data collection of standards alone we talked about.
# length(standardConfirmedMasses)
# 
# #number of unique masses shared by both spike and nonspike
# length(sorghumConfirmedMasses)
# 
# #number of rows in standards
# length(standardConfirmed$mass)
# 
# #number of rows in sorghum
# length(sorghumConfirmed$mass)
# length(sorghumConfirmed2$mass)
 


# Next Major Task: Collapsing Masses within 1 ppm of one another into a single value that represents one metabolite.
# start with standards
 
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
    #determine allowable gap between signals as being 15 parts per million
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
testcentroidFxn = c(176.1034, 176.1033,176.1032, 182.5)
centroidFxn(testcentroidFxn)

 




 
spikeValidatedMasses <- centroidFxn(standardConfirmedMasses)
nospikeValidatedMasses <- centroidFxn(sorghumConfirmedMasses)
 


 

#  
# Next step: figure out how many of the validated masses are actually real
# choose 100 of the validated masses at random
# loop through the 100 and generate plots of signal vs mass and save as pdf

 
set.seed(777624)
oneHundredValidatedMasses <- sample(spikeValidatedMasses, 100)
oneHundredValidatedMasses
set.seed(NULL)
 

 
spikeCounts <- list.files("C:/Users/Louis/Downloads/massSpec/spike", full.names = T, recursive = T)
spikeTable =  data.frame(matrix(nrow = 0, ncol = 6))
colnames(spikeTable) = c("mass", "intensity", "scan", "polarity", "massspectype", "file")
tableList = list()
for(i in 1:length(spikeCounts))
{
  #clear the memory
  gc()
  file = spikeCounts[i]
  table<- read.delim(file, sep = "\t", header = F)
  table<- table[complete.cases(table), ]
  print(i)
  print(head(file))
  TempTable =  table
  colnames(TempTable) = c("mass", "intensity", "scan", "polarity", "massspectype")
  TempTable$'NA' = NULL
  TempTable$file = spikeCounts[i]
  spikeTable = rbind(spikeTable, TempTable)
}
spikeTable = as.data.frame(spikeTable)




 