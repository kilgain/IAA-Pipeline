#Get metadata for broad dataset
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)

metadata = read.delim("C:/Users/Louis/Downloads/hmp2_metadata_to_read_in_just_metabolomics.txt", header = T)
metadata = metadata[metadata$Tube.A..Metabolomics != 'No', ]
metadata = metadata[metadata$data_type == 'metabolomics', ]
metadata = metadata[, colnames(metadata) %in% c('diagnosis', 'Tube.A..Metabolomics')]
colnames(metadata) = c('diagnosis', 'fileList')

rawData = read.csv("C:/Users/Louis/Downloads/HILp_NN_metab_named_for_louis_csv.csv", header = T)
rawData = rawData[ , 1:length(colnames(rawData))-1]






#The following code creates a metadata table used for filtering rows/columns from the raw data
{
  
  #do not want first 2 columns, as they contain mass/rt information
  sampleMap = metadata
  for(i in 1:length(sampleMap$diagnosis)) {
    if(sampleMap$diagnosis[i] == 'nonIBD') {
      sampleMap$diagnosis[i] = 0
    }
    if(sampleMap$diagnosis[i] == 'CD') {
      sampleMap$diagnosis[i] = 1
    }
    if(sampleMap$diagnosis[i] == 'UC') {
      sampleMap$diagnosis[i] = 2
    }
  }
  colnames(sampleMap) = c('diagnosis', 'originalName')
  sampleMap$originalName = str_replace(sampleMap$originalName, '-', '.')
}
#output is in sampleMap
head(sampleMap)


#Need to remove all rows that do not have a metabolite name- uninterested
rawData = rawData[rawData$Metabolite != '', ]
#Need to removal all N/A's and replace with 0's
rawData[is.na(rawData)] = 0
#resets rownames to 0-187
rownames(rawData) = NULL


#isolate the nonIBD patients, UC patients, and CD patients
controlSampleNames = sampleMap$originalName[sampleMap$diagnosis == 0]
CDSampleNames = sampleMap$originalName[sampleMap$diagnosis == 1]
UCSampleNames = sampleMap$originalName[sampleMap$diagnosis == 2]

#selects mass, rt, and all drought samples and puts  them in a new transposed dataframe (for format reasons)
rawDataNonIBD = data.frame((rawData[colnames(rawData) %in% colnames(rawData)[3:4] | colnames(rawData) %in% controlSampleNames | colnames(rawData) %in% colnames(rawData[6])]))
rawDataCD = data.frame((rawData[colnames(rawData) %in% colnames(rawData)[3:4] | colnames(rawData) %in% CDSampleNames | colnames(rawData) %in% colnames(rawData[6])]))
rawDataUC = data.frame((rawData[colnames(rawData) %in% colnames(rawData)[3:4] | colnames(rawData) %in% UCSampleNames | colnames(rawData) %in% colnames(rawData[6])]))


#stores plots
listOfPlotsPro = list()
listOfDensityPlots = list()

#iteratively loops through metabolitesOfInterest to generate a line graph of metabolite intensities by genotype
#as scaled by controlMedian
for(i in 1:length(rownames(rawDataUC))) {
  #name of metabolite
  metabolite = rawData$Metabolite[i]
  #gets median of control group metabolite
  zeroes = rawDataNonIBD[i, 4:length(colnames(rawDataNonIBD))] > 0
  controlMedian = as.numeric(rawDataNonIBD[i, 4:length(colnames(rawDataNonIBD))])
  controlMedian = median(controlMedian[zeroes])
  
    
  controlMetabolite = data.frame(t(rawDataNonIBD[i, 4:length(colnames(rawDataNonIBD))]))
  controlMetabolite$diagnosis = 'NonIBD'
  controlMetabolite$sample = rownames(controlMetabolite)
  
  
  CDMetabolite = data.frame(t(rawDataCD[i, 4:length(colnames(rawDataCD))]))
  CDMetabolite$diagnosis = 'CD'
  CDMetabolite$sample = rownames(CDMetabolite)
  
  
  UCMetabolite = data.frame(t(rawDataUC[i, 4:length(colnames(rawDataUC))]))
  UCMetabolite$diagnosis = 'UC'
  UCMetabolite$sample = rownames(UCMetabolite)
  
  
  #I want the final data to be in format (SampleName / Value / Diagnosis)
  finalTable = rbind(controlMetabolite, CDMetabolite, UCMetabolite)
  colnames(finalTable) = c('Intensity', 'Diagnosis', 'SampleName')
  finalTable$AdjustedIntensity = log10(finalTable$Intensity / controlMedian)

  
  a <- ggplot(finalTable, aes(x = AdjustedIntensity)) + ggtitle(paste("metabolite mass/name =", rawData$m.z[i], '/', metabolite, "rt =", rawData$RT..min.[i])) + ylab("Number of Samples") + xlab("log10(intensity / medianIntensity(Control))") + xlim(-3, 3) 
  
  b = a + geom_density(aes(color = Diagnosis, linetype = Diagnosis, fill = Diagnosis, alpha = .5))
  b = b + scale_color_manual(values = c("#CB4C4E", "#00AFBB", "#E7B800"))
  b = b + scale_fill_manual(values = c("#CB4C4E", "#00AFBB", "#E7B800"))
  
  
  
  # Frequency polygon:
  # Change line colors and types by groups
  a = a + geom_freqpoly( aes(color = Diagnosis, linetype = Diagnosis),
                         bins = 20, size = 1.5) +
    scale_color_manual(values = c("#CB4C4E", "#00AFBB", "#E7B800"))
  
  
  listOfPlotsPro[[i]] = a
  listOfDensityPlots[[i]] = b
}

figure = ggarrange(listOfPlotsPro[[143]],listOfPlotsPro[[119]],listOfPlotsPro[[57]],listOfPlotsPro[[28]], ncol = 1, nrow = 4, common.legend = TRUE, legend = 'bottom')
figureDensity = ggarrange(listOfDensityPlots[[143]],listOfDensityPlots[[119]],listOfDensityPlots[[57]],listOfDensityPlots[[28]], ncol = 1, nrow = 4, common.legend = TRUE, legend = 'bottom')

annotate_figure(figure, top = text_grob("PROG log-fold change under drought for metabolites across samples", color = "red", face = 'bold', size = 14))




