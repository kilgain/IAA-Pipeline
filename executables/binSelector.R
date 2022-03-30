library(stringr)


roundPlus = function(x, accuracy, f=ceiling){f(x/ accuracy) * accuracy}


args = commandArgs(trailingOnly=TRUE)

filePath = args[1]

numSamples = as.numeric(args[2])

outputPath = args[3]


rawData = read.delim(file = filePath, sep = ' ')

if(str_detect(outputPath, 'Positive')) {
	rawData = rawData[rawData$polarity == 'Positive',]
} else {
	rawData = rawData[rawData$polarity == 'Negative',]
}



rawData$Mass = round(rawData$Mass, digits = 5)

#rawData = rawData[rawData$Intensity < quantile(rawData$Intensity, probs = (.1)), ]

#this gets us all the masses WITH points. Also need to record masses with 0 points.
buckets = aggregate(rawData$Intensity, by = list(rawData$Mass), FUN = length)
colnames(buckets) = c('mass', 'count')
massesWithPoints = buckets$mass

#all masses with 0 points
allMasses = seq(min(rawData$Mass), max(rawData$Mass)+.00001, .00001)
allMasses = round(allMasses, 5)
massesWithoutPoints = setdiff(allMasses, massesWithPoints)

massesWithoutPoints = data.frame(mass = massesWithoutPoints, count = 0)


#add these masses to the table
buckets = rbind(buckets, massesWithoutPoints)

#buckets = buckets[buckets$count < quantile(buckets$count, probs = (.97)), ]
buckets = buckets[buckets$count < 2, ]

#add a loop for resampling
#get sd
#cutoff can be 2 or 3 SD above the final value

allValues = c()
for(i in 1:100) {

averageNoise = mean(sample(buckets$count, size = .1*length(buckets$mass)))
#averageNoise = mean(sample(buckets$x, size = 2000))
allValues = c(allValues, averageNoise)


}

averageNoise = mean(allValues)
sdNoise = sd(allValues)

averageNoise = averageNoise + (3 * sdNoise)

print(paste("The average noise in the reference file chosen is:", averageNoise))

possibleValues = seq(10,210,10)

totalNoise = averageNoise * numSamples

print(paste("The estimated noise level is:", totalNoise))

totalNoise = roundPlus(totalNoise, 10)

closest = possibleValues[which(abs(possibleValues - totalNoise) == min(abs(possibleValues - totalNoise)))]

write.table(closest, file = outputPath, row.names = F, col.names = F, eol = '')

