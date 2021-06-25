#!/bin/bash
#have an R script basically run this iteratively
theNum=$1
vecName=${theNum}
stringRun=$2 
make_job ()
{
  echo "
  Log        = ${vecName}.log
  Output     =  ${vecName}_tempFiles.out
  Error	   = ${vecName}.error
  request_memory = 10GB
  executable                =/usr/bin/Rscript
  arguments               =/path/to/louisRscript.R ${stringRun}
  queue" > ${vecName}_condor.sh
  #condor_submit ${vecName}_condor.sh
  #get the basename here
}
make_job