# spikeCounts <- list.files("C:/Users/Louis/Downloads/massSpec/spike", full.names = T, recursive = T)
# spikeTable =  data.frame(matrix(nrow = 0, ncol = 6))
# colnames(spikeTable) = c("mass", "intensity", "scan", "polarity", "massspectype", "file")
# tableList = list()
# for(i in 1:length(spikeCounts))
# {
#   #clear the memory
#   gc()
#   file = spikeCounts[i]
#   table<- read.delim(file, sep = "\t", header = F)
#   table<- table[complete.cases(table), ]
#   print(i)
#   print(head(file))
#   TempTable =  table
#   colnames(TempTable) = c("mass", "intensity", "scan", "polarity", "massspectype")
#   TempTable$'NA' = NULL
#   TempTable$file = spikeCounts[i]
#   spikeTable = rbind(spikeTable, TempTable)
# }
# spikeTable = as.data.frame(spikeTable)

# 
# run above code chunk in massSpecSetup.R

#will take all 2,454 masses and write them to individual tables so as to make future subsetting much faster.
spikeTable = readRDS("C:/Users/Louis/Downloads/massSpec/spikeTable.rds")
nospikeTable = readRDS("C:/Users/Louis/Downloads/massSpec/nospikeTable.rds")
BiggestTable = readRDS("C:/Users/Louis/Downloads/massSpec/All36FilesInOneTable.rds")
spikeFinalTable = data.frame(matrix(ncol = 7, nrow=0))
colnames(spikeFinalTable) = c("mass", "intensity", "scan", "polarity", "massspectype", "file", "carbonShift")
spikeDataTableSpecificMass = data.frame(matrix(ncol = 7, nrow=0))
colnames(spikeDataTableSpecificMass) = c("mass", "intensity", "scan", "polarity", "massspectype", "file", "carbonShift")
for(i in 1:1) {
myMass <- 182.1230
ppmAllowed = .000015
rangeAllowed = ppmAllowed * myMass
upperLimit = myMass + rangeAllowed
lowerLimit = myMass - rangeAllowed
spikeDataTableSpecificMass = spikeTable %>% filter(mass <  upperLimit & mass > lowerLimit)
spikeDataTableSpecificMass$carbonShift = 0
spikeFinalTable = rbind(spikeFinalTable, spikeDataTableSpecificMass)
# write.table(spikeDataTableSpecificMass, paste("pointsWithin5PPMof",myMass[i],".txt", sep = ""), row.names = F)
}

spikeTable <- as.data.table(spikeTable)
nospikeTable <- as.data.table(nospikeTable)
list = list(spikeTable, nospikeTable)
BiggestTable <- as.data.frame(rbindlist(list))
saveRDS(BiggestTable, "All36FilesInOneTable.rds")
#Question: we will want to run this on ALL masses, before boundary-finding, right?
#The purpose of isolating the individual masses is so we can more easily separate up work on server for running
#the boundary-finding script I wrote
