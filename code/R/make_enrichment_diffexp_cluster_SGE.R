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
synapseLogin()

# Get differential methylation (ALL)
diffexp.id = 'syn5231514'
diffexp.genesets = downloadFile(diffexp.id)
clusters = diffexp.genesets$cluster

# Make directory and write shell scripts for running these files
system('mkdir sgeEnrichDiffexp')
fp_all = file(paste('sgeEnrichDiffexp/allSubmissions.sh'),'w+')    
cat('#!/bin/bash',file=fp_all,sep='\n')
close(fp_all)

for (clust in clusters) {
  fp = file (paste('/mnt/Github/pcbc_c4_analysis/code/R/sgeEnrichDiffexp/SUB',clust,'sh',sep='.'), "w+")
  cat('#!/bin/bash', 
      'sleep 30', 
      paste('Rscript /mnt/Github/pcbc_c4_analysis/code/R/enrichment_diffexp_clusters_SGE.R','syn5231514',clust,'syn5834610'), 
      file = fp,
      sep = '\n')
  close(fp)
  
  fp_all = file(paste('sgeEnrichDiffexp/allSubmissions.sh'),'a+')    
  cat(paste('sh', paste('/mnt/Github/pcbc_c4_analysis/code/R/sgeEnrichDiffexp/SUB',clust,'sh',sep='.')),
      file=fp_all,
      sep='\n')
  close(fp_all)
}