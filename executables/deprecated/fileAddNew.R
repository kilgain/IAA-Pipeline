library(stringr)

file = commandArgs(trailingOnly = T)
print(file)
print("is file")

table =  read.delim(file, sep = "\t", header = F)
table = table[complete.cases(table), ]

if(ncol(table) < 6) {
 table = cbind(basename(file), table)
}

if(str_detect(table[1, 5], 'PolarityType')) {
 table[,5] = str_replace(table[,5], 'PolarityType.', '')
}


if(str_detect(table[1,6], 'MsOrderType')){
 table[,6] = str_replace(table[,6], 'MsOrderType.', '')
}

colnames(table) = c("File", "Mass", "Intensity", "scan", "polarity", "Mstype")





write.table(table, file=file, row.names = F)
