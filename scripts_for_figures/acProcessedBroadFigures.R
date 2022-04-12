#Get metadata for broad dataset
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)


#Read in data:
metadata = read.delim("C:/Users/Louis/Downloads/hmp2_metadata_to_read_in_just_metabolomics.txt", header = T)

metadata = metadata[metadata$Tube.A..Metabolomics != 'No', ]

metadata = metadata[metadata$data_type == 'metabolomics', ]

metadata = metadata[, colnames(metadata) %in% c('diagnosis', 'Tube.A..Metabolomics')]

colnames(metadata) = c('diagnosis', 'fileList')

files = list.files('C:/Users/Louis/OneDrive/Desktop/broadDataNamedMetabolites')

rawData = read.csv(paste('C:/Users/Louis/OneDrive/Desktop/broadDataNamedMetabolites/', files[1], sep = ''))

for(f in 2:4) {
  rawData = rbind(rawData, read.csv(paste('C:/Users/Louis/OneDrive/Desktop/broadDataNamedMetabolites/', files[f], sep = '')))
}


#formatting
columns = str_which(colnames(rawData[, 1:length(colnames(rawData))]), "SM.")
columns = c(1,2,columns)
rawData = rawData[ ,columns]
names = str_extract(colnames(rawData[, 3:length(colnames(rawData))]), 'SM......')
colnames(rawData) = c('compound', 'rt', names)




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


#resets rownames
rownames(rawData) = NULL


#isolate the nonIBD patients, UC patients, and CD patients
controlSampleNames = sampleMap$originalName[sampleMap$diagnosis == 0]
CDSampleNames = sampleMap$originalName[sampleMap$diagnosis == 1]
UCSampleNames = sampleMap$originalName[sampleMap$diagnosis == 2]

#selects mass, rt, and all drought samples and puts  them in a new transposed dataframe (for format reasons)
rawDataNonIBD = data.frame(rawData[colnames(rawData) %in% colnames(rawData)[1:2] | colnames(rawData) %in% controlSampleNames])
rawDataCD = data.frame((rawData[colnames(rawData) %in% colnames(rawData)[1:2] | colnames(rawData) %in% CDSampleNames]))
rawDataUC = data.frame((rawData[colnames(rawData) %in% colnames(rawData)[1:2] | colnames(rawData) %in% UCSampleNames]))


#stores plots
listOfPlotsAC = list()
listOfDensityPlots = list()
#iteratively loops through metabolitesOfInterest to generate a line graph of metabolite intensities by genotype
#as scaled by controlMedian
for(i in 1:length(rownames(rawDataUC))) {

  #gets median of control group metabolite
  zeroes = rawDataNonIBD[i, 3:length(colnames(rawDataNonIBD))] > 0
  controlMedian = as.numeric(rawDataNonIBD[i, 3:length(colnames(rawDataNonIBD))])
  controlMedian = median(controlMedian[zeroes])
  
  #Control data
  controlMetabolite = data.frame(t(rawDataNonIBD[i, 3:length(colnames(rawDataNonIBD))]))
  controlMetabolite$diagnosis = 'NonIBD'
  controlMetabolite$sample = rownames(controlMetabolite)
  
  #CD data
  CDMetabolite = data.frame(t(rawDataCD[i, 3:length(colnames(rawDataCD))]))
  CDMetabolite$diagnosis = 'CD'
  CDMetabolite$sample = rownames(CDMetabolite)
  
  #UC data
  UCMetabolite = data.frame(t(rawDataUC[i, 3:length(colnames(rawDataUC))]))
  UCMetabolite$diagnosis = 'UC'
  UCMetabolite$sample = rownames(UCMetabolite)
  
  
  #final data table in format (Intensity / Diagnosis / SampleName)
  finalTable = rbind(controlMetabolite, CDMetabolite, UCMetabolite)
  colnames(finalTable) = c('Intensity', 'Diagnosis', 'SampleName')
  finalTable$AdjustedIntensity = log10(finalTable$Intensity / controlMedian)
  
  
  a <- ggplot(finalTable, aes(x = AdjustedIntensity)) + ggtitle(paste("metabolite mass/name =", rawData$compound[i], '/', metabolite, "rt =", rawData$rt[i])) + ylab("Number of Samples") + xlab("log10(intensity / medianIntensity(Control))") + xlim(-3, 3) 
  
  #first representation of data
  b = a + geom_density(aes(color = Diagnosis, linetype = Diagnosis, fill = Diagnosis, alpha = .5))
  b = b + scale_color_manual(values = c("#CB4C4E", "#00AFBB", "#E7B800"))
  b = b + scale_fill_manual(values = c("#CB4C4E", "#00AFBB", "#E7B800"))
  
  
  #second representation of data
  # Frequency polygon:
  # Change line colors and types by groups
  a = a + geom_freqpoly( aes(color = Diagnosis, linetype = Diagnosis),
                         bins = 20, size = 1.5) +
    scale_color_manual(values = c("#CB4C4E", "#00AFBB", "#E7B800"))
  
  
  listOfPlotsAC[[i]] = a
  listOfDensityPlots[[i]] = b
}

figureAC = ggarrange(listOfPlotsAC[[2]],listOfPlotsAC[[3]], listOfPlotsAC[[6]], listOfPlotsAC[[1]], ncol = 1, nrow = 4, common.legend = TRUE, legend = 'bottom')
figureDensity = ggarrange(listOfDensityPlots[[2]],listOfDensityPlots[[3]],listOfDensityPlots[[6]],listOfDensityPlots[[1]], ncol = 1, nrow = 4, common.legend = TRUE, legend = 'bottom')

annotate_figure(figureAC, top = text_grob("log-fold change under drought for metabolites across samples", color = "red", face = 'bold', size = 14))




