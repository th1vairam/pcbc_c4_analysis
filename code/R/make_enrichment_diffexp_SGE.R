#!usr/bin/env Rscript

# Submission Script in R
# Clear R console screen output
cat("\014")

# Clear R workspace
setwd('/mnt/Github/pcbc_c4_analysis/code/R')

# Load libraries
library(synapseClient)
library(CovariateAnalysis)

# Login to synapse
# apikey = read.table('/mnt/synapseAPIKey', stringsAsFactors = F)
# synapseLogin(username = 'th_vairam', apiKey = apikey$V1)
synapseLogin()

# Get differential methylation (ALL)
diffexp.id = 'syn5706668'
diffexp = downloadFile(diffexp.id)
ind = which(diffexp$setname %in% setdiff(diffexp$setname, tmp$file.name))

# Make directory and write shell scripts for running these files
system('mkdir sgeEnrichDiffexp')
fp_all = file(paste('sgeEnrichDiffexp/allSubmissions.sh'),'w+')    
cat('#!/bin/bash',file=fp_all,sep='\n')
close(fp_all)

for (id in ind) {
  fp = file (paste('/mnt/Github/pcbc_c4_analysis/code/R/sgeEnrichDiffexp/SUB',id,'sh',sep='.'), "w+")
  cat('#!/bin/bash', 
      'sleep 30', 
      paste('Rscript /mnt/Github/pcbc_c4_analysis/code/R/enrichment_diffexp_SGE.R',diffexp.id,id), 
      file = fp,
      sep = '\n')
  close(fp)
  
  fp_all = file(paste('sgeEnrichDiffexp/allSubmissions.sh'),'a+')    
  cat(paste('sh', paste('/mnt/Github/pcbc_c4_analysis/code/R/sgeEnrichDiffexp/SUB',id,'sh',sep='.')),
      file=fp_all,
      sep='\n')
  close(fp_all)
}