#get list of valid masses from massSpecAnalysisScript
#massesToKeep

#for each mass, start by subtracting the mass of one carbon 13 neutron
#1.003355

#then, record the number of points within 3-5 ppm of original mass at that point in the mass spectra data
#of sorghumConfirmedMasses

#iterate the mass of another carbon and record the number of points within 3 ppm of original mass
#compare with previous recorded number
#if higher, keep new value
#if lower, iterate one more time to make sure- but keep old value

#iterate until the value is lower than the previous, then plot all points within 5 ppm of the identified mass
#using the sorghumConfirmed dataset

#iteration limit using mass divided by number of carbon atoms

metaboliteMatcher = function(massesToAnalyze, sorghumMetabolites, listOfVecOfStarts, listOfVecOfStops) {
  # pdf(file = "C:/Users/Louis/OneDrive/Desktop/100Plots/BarPlotsTestGluco.pdf",
  #     width = 8,
  #     height = 5)
listOfVecOfStarts = listOfVecOfStarts
listOfVecOfStops = listOfVecOfStops
massesToAnalyze = massesToAnalyze
sorghumMetabolites = sorghumMetabolites
listOfVecOfNumPoints = list()
for(m in 1:1) {
  vecOfStarts = listOfVecOfStarts[[m]]
  vecOfStops = listOfVecOfStops[[m]]
  # listOfNumPointsAtDifMassShifts = list()
  for(l in 1:length(vecOfStarts)) {
  vecOfNumPointsAtDifMassShifts = vector()
  counter = 1
  massOfInterest <- massesToAnalyze[m]
  maxCounter = massOfInterest / 12
  ppmAllowed = .000005
  rangeAllowed = ppmAllowed * massOfInterest
  while(counter < maxCounter) {
    massOfSorghumMetabolite <- massOfInterest - (counter*1.00336)
    sorghumRegionOfInterest <- sorghumMetabolites %>% filter(mass < (massOfSorghumMetabolite+rangeAllowed)& mass > (massOfSorghumMetabolite-rangeAllowed))
    sorghumRegionOfInterest <- sorghumRegionOfInterest %>% filter( scan > vecOfStarts[l] & scan < vecOfStops[l])
    vecOfNumPointsAtDifMassShifts = c(vecOfNumPointsAtDifMassShifts, length(sorghumRegionOfInterest$mass))
    # listOfNumPointsAtDifMassShifts = vecOfNumPointsAtDifMassShifts
    counter = counter + 1
  }
  }
  listOfVecOfNumPoints[[m]] = vecOfNumPointsAtDifMassShifts
  x <- 1:(counter-1)
  barplot(height = vecOfNumPointsAtDifMassShifts, names=x, xlab = "number of carbon 13's substituted out", ylab = "number of points within 5ppm of yeast standard's mass", main = paste("Metabolite Matcher for Mass", massOfInterest))
  # mySubset = spikeTable %>% filter(mass < (myMass+rangeAllowed) & mass > (myMass-rangeAllowed))
}
# dev.off()
return(listOfVecOfNumPoints)
}

if(length(vecOfStarts = 0)) {vecOfStarts = 0}
if(length(vecOfStops = 0)) {vecOfStops=3000}
numPointsAtEachShift <- metaboliteMatcher(myMass+.0012, BiggestTable, vecOfStarts, vecOfStops)

# metaboliteChecker = function(massesToKeep, numPointsAtEachShift, sorghumMetabolites, standards) {
#   pdf(file = "C:/Users/Louis/OneDrive/Desktop/100Plots/CheckRetentionTimeTest.pdf",
#       width = 8,
#       height = 5) 
#   numPointsAtEachShift = numPointsAtEachShift
#   massesToKeep = massesToKeep
#   standards = standards
#   sorghumMetabolites = sorghumMetabolites
#   for(m in 1:1) {
#     allMassesOfInterest = data.frame(matrix(ncol = 7))
#     colnames(allMassesOfInterest) = c("mass", "intensity", "scan", "polarity", "massspectype", "file", "carbonShift")
#     numPoints = numPointsAtEachShift[[m]]
#     massOfInterest <- massesToKeep[m]
#     ppmAllowed = .000005
#     rangeAllowed = ppmAllowed * massOfInterest
#     for(p in 1:length(numPoints)) {
#       if(numPoints[p] > 200) {
#         massOfSorghumMetabolite <- massOfInterest - (p*1.00336)
#         sorghumRegionOfInterest <- sorghumMetabolites %>% filter(mass < (massOfSorghumMetabolite+rangeAllowed)& mass > (massOfSorghumMetabolite-rangeAllowed))
#         sorghumRegionOfInterest <- sorghumRegionOfInterest %>% filter(scan < (vecOfStops[1])& scan > (vecOfStarts[1]))
#         sorghumRegionOfInterest$carbonShift = p
#         head(sorghumRegionOfInterest)
#         allMassesOfInterest <- rbind(allMassesOfInterest, sorghumRegionOfInterest)
#       }
#     }
#     standardsRegionOfInterest <- standards %>% filter(mass < (massOfInterest+rangeAllowed) & mass > (massOfInterest-rangeAllowed))
#     standardsRegionOfInterest <- standards %>% filter(scan < (vecOfStops[1])& scan > (vecOfStarts[1]))
#     standardsRegionOfInterest$carbonShift = 0
#     head(standardsRegionOfInterest)
#     allMassesOfInterest <- rbind(allMassesOfInterest, standardsRegionOfInterest)
#     # pointsOfInterest <- allMassesOfInterest %>% filter(carbonShift == 0)
#     if(length(allMassesOfInterest$scan) > 300) {
#     plot(
#       allMassesOfInterest$scan,
#       allMassesOfInterest$intensity,
#       xlab = "Scan",
#       xlim = c(-10, 3100),
#       ylab = "Intensity",
#       col = as.factor(allMassesOfInterest$carbonShift),
#       main = paste(
#         "Retention Time for Mass",
#         massesToKeep[m]
#       )
#     )
#       legend("topleft", legend=levels(as.factor(allMassesOfInterest$carbonShift)), text.col=seq_along(levels(factor(allMassesOfInterest$carbonShift))))
#     }
#   }
#   dev.off()
# }
# metaboliteChecker(182.1230+.0012, numPointsAtEachShift, nospikeTable, spikeFinalTable)


p=5
myMass = 121.0874
rangeAllowed = .000005 * myMass
massOfSorghumMetabolite = myMass - (p)*(1.00336)
sorghumRegionOfInterest <- BiggestTable %>% filter(mass < (massOfSorghumMetabolite+rangeAllowed)& mass > (massOfSorghumMetabolite-rangeAllowed))
sorghumRegionOfInterest <- sorghumRegionOfInterest %>% filter(scan < (vecOfStops[1])& scan > (vecOfStarts[1]))
massOfSorghumMetabolite <- mean(sorghumRegionOfInterest$mass)

M0 <- BiggestTable[BiggestTable$mass>((massOfSorghumMetabolite)-((massOfSorghumMetabolite)*(.000005))) & BiggestTable$mass < ((massOfSorghumMetabolite)+((massOfSorghumMetabolite)*(.000005))),]
M0 <- M0 %>% filter(scan < (vecOfStops[1])& scan > (vecOfStarts[1]))
# write.csv(M0, "M0ForAcetamidopropanal.csv", row.names = F)
{
  plot(
    M0$scan,
    M0$intensity,
    xlab = "Scan",
    xlim = c(-10, 3100),
    ylab = "Intensity",
    main = paste(
      "Retention Time for Mass", massOfSorghumMetabolite
    )
  )
  plot(
    M0$mass,
    M0$intensity,
    xlab = "Mass",
    # xlim = c(-10, 3100),
    ylab = "Intensity",
    main = paste(
      "Mass", massOfSorghumMetabolite
    )
  )
}
maxSignalM0 <- max(M0$intensity)


M1 <- BiggestTable[BiggestTable$mass>((massOfSorghumMetabolite+1.00336)-((massOfSorghumMetabolite+1.00336)*(.000005))) & BiggestTable$mass < ((massOfSorghumMetabolite+1.00336)+((massOfSorghumMetabolite+1.00336)*(.000005))),]
M1 <- M1 %>% filter(scan < (vecOfStops[1])& scan > (vecOfStarts[1]))
# write.csv(M1, "M1ForAcetamidopropanal.csv", row.names = F)
# write.table(citrulline1CarbonStandard, "OneCarbon13UsingMass176.104.csv")
#plotter for 176.104+1.00336
{
  plot(
    M1$scan,
    M1$intensity,
    xlab = "Scan",
    xlim = c(-10, 3100),
    ylab = "Intensity",
    main = paste(
      "Retention Time for Mass", massOfSorghumMetabolite, "+1.00336"
    )
  )
  plot(
    M1$mass,
    M1$intensity,
    xlab = "Mass",
    # xlim = c(-10, 3100),
    ylab = "Intensity",
    main = paste(
      "Mass", massOfSorghumMetabolite, "+1.00336"
    )
  )
}
maxSignalM1 <- max(M1$intensity)
print(paste("number of carbons in this compound is equal to:", (maxSignalM1/maxSignalM0)*100/1.109))

vecOfRatios = vector()
for(i in 1:10000)
{
  print(i)
  subsetM0 = M0[sample(nrow(M0), (100)),]
  subsetM1 = M1[sample(nrow(M1), (100)),]
  M1vM0ratio = max(subsetM1$intensity) / max(subsetM0$intensity)
  vecOfRatios = c(vecOfRatios,  M1vM0ratio)
}
vecOfRatios <- sort(vecOfRatios, decreasing = T)
bins <- cut(vecOfRatios, 100)
bins <- sort(summary(bins), decreasing = T)[1]
binBounds <- as.numeric(unlist(regmatches(
  names(bins),
  gregexpr("[[:digit:]]+\\.*[[:digit:]]*", names(bins))
)))
valueToPlot <- mean(binBounds)

hist(vecOfRatios, breaks = 250, main = paste("VectorOfRatios For Mass=", massOfSorghumMetabolite, " carbon shift of ", p, " from ", myMass))
abline(v = p*.0109, col = "green")
abline(v = valueToPlot, col = "red")



p=11
myMass = 216.134
massOfSorghumMetabolite = myMass - (p)*(1.00336)
sorghumRegionOfInterest <- BiggestTable %>% filter(mass < (massOfSorghumMetabolite+rangeAllowed)& mass > (massOfSorghumMetabolite-rangeAllowed))
sorghumRegionOfInterest <- sorghumRegionOfInterest %>% filter(scan < (vecOfStops[1])& scan > (vecOfStarts[1]))
massOfSorghumMetabolite <- mean(sorghumRegionOfInterest$mass)

ElevenCarbon <- BiggestTable[BiggestTable$mass>((massOfSorghumMetabolite)-((massOfSorghumMetabolite)*(.000005))) & BiggestTable$mass < ((massOfSorghumMetabolite)+((massOfSorghumMetabolite)*(.000005))),]
# ElevenCarbon <- ElevenCarbon %>% filter(scan < (vecOfStops[1])& scan > (vecOfStarts[1]))
{
  plot(
    ElevenCarbon$scan,
    ElevenCarbon$intensity,
    xlab = "Scan",
    # xlim = c(-10, 3100),
    ylab = "Intensity",
    main = paste(
      "Retention Time for Mass 205.0981"
    )
  )
  plot(
    ElevenCarbon$mass,
    ElevenCarbon$intensity,
    xlab = "Mass",
    # xlim = c(-10, 3100),
    ylab = "Intensity",
    main = paste(
      "Mass 205.0981"
    )
  )
}
maxSignalElevenCarbonNoC13 <- max(ElevenCarbon$intensity)
ElevenCarbon1C13 <- BiggestTable[BiggestTable$mass>((massOfSorghumMetabolite+1.00336)-((massOfSorghumMetabolite+1.00336)*(.000005))) & BiggestTable$mass < ((massOfSorghumMetabolite+1.00336)+((massOfSorghumMetabolite+1.00336)*(.000005))),]
# ElevenCarbon1C13 <- ElevenCarbon1C13 %>% filter(scan < (vecOfStops[1])& scan > (vecOfStarts[1]))
{
  plot(
    ElevenCarbon1C13$scan,
    ElevenCarbon1C13$intensity,
    xlab = "Scan",
    # xlim = c(-10, 3100),
    ylab = "Intensity",
    main = paste(
      "Retention Time for Mass 205.0981+1.00336"
    )
  )
  plot(
    ElevenCarbon1C13$mass,
    ElevenCarbon1C13$intensity,
    xlab = "Mass",
    # xlim = c(-10, 3100),
    ylab = "Intensity",
    main = paste(
      "Mass 205.0981+1.00336"
    )
  )
}
maxSignalElevenCarbonOneC13 <- max(ElevenCarbon1C13$intensity)

print(paste("number of carbons in this compound is equal to:", (maxSignalElevenCarbonOneC13/maxSignalElevenCarbonNoC13)*100/1.109))


# citrulline1CarbonStandard <- spikeTable[spikeTable$mass>((176.1034+1.00336+.0012)-((176.1034+1.00336+.0006)*(.000005))) & spikeTable$mass < ((176.1034+1.00336+.0006)+((176.1034+1.00336+.0006)*(.000005))),]
# # write.table(citrulline1CarbonStandard, "ZeroCarbon13UsingMass176.1034.csv")
# maxSignalZeroCarbon13 <- max(citrulline1CarbonStandard$intensity)
# #plotter for 176.1034 + .0012 + 1.00336
# {
# plot(
#   citrulline1CarbonStandard$scan,
#   citrulline1CarbonStandard$intensity,
#   xlab = "Scan",
#   # xlim = c(-10, 3100),
#   ylab = "Intensity",
#   main = paste(
#     "Retention Time for Mass 176.1034+1.00336+.0012"
#   )
# )
# plot(
#   citrulline1CarbonStandard$mass,
#   citrulline1CarbonStandard$intensity,
#   xlab = "Mass",
#   # xlim = c(-10, 3100),
#   ylab = "Intensity",
#   main = paste(
#     "Mass 176.1034+1.00336+.0012"
#   )
# )
# }

#MassesThatWork = 121.0874 (acetamidopropanal), 182.1230

pdf(file = "C:/Users/Louis/OneDrive/Desktop/100Plots/CheckMassTestAll15.pdf",
    width = 8,
    height = 5) 
for(p in 1:15) {
  massOfInterest = 182.1230+.0012
  rangeAllowed = .000005 * massOfInterest
  massOfSorghumMetabolite <- massOfInterest - (p*1.00336)
  sorghumRegionOfInterest <- sorghumMetabolites %>% filter(mass < (massOfSorghumMetabolite+rangeAllowed)& mass > (massOfSorghumMetabolite-rangeAllowed))
  sorghumRegionOfInterest <- sorghumRegionOfInterest %>% filter(scan < (vecOfStops[1])& scan > (vecOfStarts[1]))
  # if(length(sorghumRegionOfInterest$intensity) > 0) {
  # signal = max(sorghumRegionOfInterest$intensity)
  # }
  # signalBackground <-
  #   mean(sorghumRegionOfInterest$intensity[sorghumRegionOfInterest$intensity < quantile(sorghumRegionOfInterest$intensity, probs = .1)])
  # if(is.nan(signal/signalBackground)) {
  #   signal = 0
  #   signalBackground = 1
  # }
  if(length(sorghumRegionOfInterest$intensity) > 50) {
    plot(sorghumRegionOfInterest$mass, sorghumRegionOfInterest$intensity, main = paste("number of carbon atoms removed:", p))
  }
}
dev.off()

