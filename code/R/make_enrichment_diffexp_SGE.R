#!usr/bin/env Rscript

# Submission Script in R
# Clear R console screen output
cat("\014")

# Clear R workspace
setwd('/shared/Github/pcbc_c4_analysis/code/R')

# Load libraries
library(synapseClient)

# Login to synapse
apikey = read.table('/shared/synapseAPIKey')
synapseLogin(apiKey = apikey)

# Get differential methylation (ALL)
diffexp.id = 'syn5706668'
diffexp = downloadFile(diffexp.id)

# Make directory and write shell scripts for running these files
system('mkdir sgeEnrichDiffexp')
fp_all = file(paste('sgeEnrichDiffexp/allSubmissions.sh'),'w+')    
cat('#!/bin/bash',file=fp_all,sep='\n')
close(fp_all)

for (id in 1:10) {
  fp = file (paste('/shared/Github/pcbc_c4_analysis/code/R/sgeEnrichDiffexp/SUB',id,sep='.'), "w+")
  cat('#!/bin/bash', 
      'sleep 30', 
      paste('Rscript /shared/Github/pcbc_c4_analysis/code/R/enrichment_diffexp_SGE.R',id), 
      file = fp,
      sep = '\n')
  close(fp)
  
  fp_all = file(paste('sgeEnrichDiffexp/allSubmissions.sh'),'a+')    
  cat(paste('qsub','-cwd','-V',paste('/shared/Github/pcbc_c4_analysis/code/R/sgeEnrichDiffexpDiffexp/SUB',id,sep='.'),
            '-o',paste('/shared/Github/pcbc_c4_analysis/code/R/sgeEnrichDiffexp/SUB',id,'o',sep='.'),
            '-e',paste('/shared/Github/pcbc_c4_analysis/code/R/sgeEnrichDiffexp/SUB',id,'e',sep='.')),
      file=fp_all,
      sep='\n')
  close(fp_all)
}