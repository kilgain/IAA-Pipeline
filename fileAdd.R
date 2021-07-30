
library(yaml)

file = commandArgs(trailingOnly = T)
print(file)
print("is file")

config <- read_yaml("/home/lconnelly/Metabolomics/Config.yaml")

file = paste(config$samplesDirectory, "/", file, sep="") 

table =  read.delim(file, sep = "\t", header = F)
table = table[complete.cases(table), ]

if(ncol(table) < 6) {
 table = cbind(basename(file), table)
}

colnames(table) = c("File", "Mass", "Intensity", "scan", "polarity", "Mstype")





write.table(table, file=file, row.names = F)
