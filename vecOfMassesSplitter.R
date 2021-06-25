


spikeValidatedMasses <- read.csv("C:/Users/Louis/OneDrive/Desktop/vecOfMasses.csv")
chunksOfMasses <- split(spikeValidatedMasses$mass, ceiling(seq_along(spikeValidatedMasses$mass)/50))
i = 1
for(chunk in chunksOfMasses) {
  chunkToString = vector()
  chunkToString = chunk
  chunkToRunString = paste(chunkToString, collapse = ",")
  print(i)
  #so my R script should also take "i" as an argument to be able to print out progress, presumably
  theString = paste("bash /home/lconnelly/matcher_parallelize.sh", i, chunkToRunString)
  print(theString)
  #so this is what actually executes the R script with the arguments provided in the above lines
  system(theString)
}

#I think I am missing a chunk of code between this ^ and the make_job() script that takes theString and 
#runs it using the job scheduler?

#How does this script relate to the make_job() script? How do they interact?
#what type of script is make_job? what is {vecName}?
#what would be my equivalent to step5_call_all_the_mass_features/parallelize_the_masses_dynamic_time_warping_files.sh?
#it would be a bash script that runs my R script with the specified arguments

