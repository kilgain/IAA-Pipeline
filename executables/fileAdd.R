

file = commandArgs(trailingOnly = T)
#print(file)
#print("is file")

#Read in file
table =  read.delim(file, sep = "\t", header = F)
table = table[complete.cases(table), ]

#check if the 'file' column is present, if not add it
if(ncol(table) < 6) {
 table = cbind(basename(file), table)
}

#rename columns
colnames(table) = c("File", "Mass", "Intensity", "scan", "polarity", "Mstype")




#writeout
write.table(table, file=file, row.names = F)
