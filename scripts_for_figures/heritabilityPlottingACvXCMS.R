#setaria/xcms density plot of heritability
allReadIn <- read.table(file = "C:/Users/Louis/Downloads/prelim_DOE_intensity_table.csv", sep = ",", header = TRUE)
massVec = allReadIn$compound



rownames(allReadIn) = make.unique(paste(allReadIn$compound,allReadIn$rt, sep = "_"))
allReadIn$compound = NULL
allReadIn$rt = NULL



lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
vecOfHeritabilityHits = vector()
unExplainedVariance = vector()
heritability = vector()

fileList = colnames(allReadIn)

fileListSpecies = vector()
fileListSplit = lapply(fileList, `[`,1)
fileListSplit = unlist(fileListSplit)


#fill in the drought info
droughtInfo = fileListSplit
droughtInfo[grepl("40",droughtInfo, fixed = TRUE)] = 1
droughtInfo[grepl("100",droughtInfo, fixed = TRUE)] = 0



#fill in the genotype info
genoInfo = fileListSplit
genoInfo[grepl("TB12.48",genoInfo, fixed = TRUE)] = 1
genoInfo[grepl("A10",genoInfo, fixed = TRUE)] = 2
genoInfo[grepl("TB12.201",genoInfo, fixed = TRUE)] = 3



#fill in the time info
timeInfo = fileListSplit
timeInfo[grepl("T4",timeInfo, fixed = TRUE)] = 1
timeInfo[grepl("T6",timeInfo, fixed = TRUE)] = 2
timeInfo[grepl("T8",timeInfo, fixed = TRUE)] = 3


#combine everything now
forModel = cbind(droughtInfo,genoInfo,timeInfo)

forModel = as.data.frame(forModel)

dataFrameAllFragments = forModel 
dataFrameAllFragments$originalName = fileList

#Match expression values to the metadata
#dataFrameAllFragments <- dataFrameAllFragments[match(colnames(expression),dataFrameAllFragments$originalName),]


#step3- run model
output = as.data.frame(matrix(nrow = 1, ncol = 7))
colnames(output) = c('mass', 'pTime', 'pTreat', 'pGeno', 'pTixTr', 'pTxG', 'heritability')


#length(rownames(expressionSet))
for(m in 1:length(rownames(allReadIn))){
  expressionSub = unlist(allReadIn[m,])
  heritability = vector()
  if(sum(expressionSub) > 0)
  {
    #keep track of all of the dummy (categorical) variables that will be incorporated into the model
    dummyVariablesGenotype = as.factor(dataFrameAllFragments$genoInfo)
    dummyVariablesTreatment = as.factor(dataFrameAllFragments$droughtInfo)
    #dummyVariablesSpecies = as.factor(dataFrameAllFragments$species)
    dummyVariablesTime = as.factor(dataFrameAllFragments$timeInfo)
    #make this a bit more complex
    #include the time (as E)
    #and also the GxE equivalent, or Accession * Time interaction
    #include the effect of genotype (across all times + treatments, does genotype have an effect)
    #effect of time (across all genotypes and treatmens, does the time create an effect)
    #effect of treatment (across all genotypes and times, does treatment exert an effect)
    #as well as the more subtler interaction terms:
    #on average across genotypes, do different genotypes respond differently to drought and control conditions (GxE)
    #or, in our situation, dummyVariablesGenotype*dummyVariablesGenotype
    #on average, across all timepoints, does this GxE change across time?  Model this by dummyVariablesGenotype*dummyVariablesTreatment*dummy$
    
    
    
    model = lm(expressionSub ~  dummyVariablesTime + dummyVariablesTreatment + dummyVariablesGenotype  + dummyVariablesGenotype*dummyVariablesTreatment* dummyVariablesTime)
    
    
    
    pValues = lmp(model)
    anovaTest = anova(model)
    
    
    if(pValues < .05)
    {
      vecOfHeritabilityHits= c(vecOfHeritabilityHits, m)
      unExplainedVariance = c(unExplainedVariance, mean(abs(model$residuals)))
    }
    heritability =  c(heritability,sum(anovaTest$`Sum Sq`[1:5]) / sum(anovaTest$`Sum Sq`))
  }
  if(sum(expressionSub) > 0) {
    
    pvalues = anovaTest$`Pr(>F)`
    informationForOutputRow = c(massVec[m], pvalues[1], pvalues[2], pvalues[3], pvalues[4], pvalues[5], heritability[1])
    
    output = rbind(output, informationForOutputRow)
  }
  
  #if you want to write out the outputs:
  # if(m == length(rownames(allReadIn))) {
  #   output = output[-1,]
  #   write.csv(output, file = "fullHeritabilityOutputTableAllPrelim_AC_7April.csv", row.names=F)
  # }
  
}


ACHeritabilities = output$heritability


outputAC = output



set.seed(1234)

#Bind together the different heritabilities

outputXCMS = read.csv("C:/Users/Louis/Downloads/fullHeritabilityOutputTableAllPrelim_XCMS.csv")
xcmsHeritabilities = outputXCMS$heritability

library(ggplot2)
library(ggpubr)





#ACHeritabilitiesHigh = ACHeritabilities[ACHeritabilities > .2]
ACInfo  = cbind(ACHeritabilities, rep("AC", length(ACHeritabilities )))
xcmsInfo = cbind(xcmsHeritabilities, rep("XCMS", length(xcmsHeritabilities )))

heritabilityTable = data.frame(
  rbind(ACInfo,xcmsInfo )
)
head(heritabilityTable, 4)

colnames(heritabilityTable) = c("Heritability", "Method")

#make sure the heritability table has the heritability column as numeric
heritabilityTable$Heritability = as.numeric(heritabilityTable$Heritability)

library("dplyr")
mu <- heritabilityTable %>% 
  group_by(Method) %>%
  summarise(grp.mean = mean(Heritability))
mu


a <- ggplot(heritabilityTable, aes(x = Heritability)) + ggtitle("Explained Variance in Detected Metabolites Depending on Detection Method") + ylab("Number of Mass Features") + xlab("Proportion of Variance Explained by Model")





# Frequency polygon: 
# Change line colors and types by groups
a + geom_freqpoly( aes(color = Method, linetype = Method),
                   bins = 30, size = 1.5) +
  scale_color_manual(values = c("#00AFBB", "#E7B800"))
# Area plots: change fill colors by Method
# Create a stacked area plots
a + geom_area(aes(fill = Method), color = "white", 
              stat ="bin", bins = 30) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800"))