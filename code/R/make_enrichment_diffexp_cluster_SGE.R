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

# Get github links for provenance
thisFileName <- 'make_diffexp_clusters_SGE.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/pcbc_c4_analysis", 
                    ref="branch", 
                    refName='enrich')

thisFile <- getPermlink(repository = thisRepo,
                        repositoryPath=paste0('code/R/', thisFileName))

# Make submissions for mRNA
mrna.id = 'syn5231514'
mrna.genesets = downloadFile(mrna.id)

allInput = ddply(mrna.genesets, .(cluster),.fun = function(x){
  data.frame(genesets = x$feature,
             cluster = x$cluster,
             parentId = 'syn5834610',
             executed = thisFile,
             source.file = mrna.id,
             fileName = paste('Cluster', x$cluster))
})

# Make submissions for miRNA
mirna.id = 'syn5834565'
mirna.genesets = downloadFile(mirna.id)

allInput = rbind( allInput, 
  ddply(mirna.genesets, .(cluster),.fun = function(x){
    data.frame(genesets = x$feature,
               cluster = x$cluster,
               parentId = 'syn5836396',
               executed = thisFile,
               source.file = mirna.id,
               fileName = paste('Cluster', x$cluster))
    }))

# Make submissions for methyl
methyl.id = 'syn5834569'
methyl.genesets = downloadFile(methyl.id)

allInput = rbind( allInput, 
                  ddply(methyl.genesets, .(cluster),.fun = function(x){
                    data.frame(genesets = x$feature.concordant,
                               cluster = x$cluster,
                               parentId = 'syn5838771',
                               executed = thisFile,
                               source.file = methyl.id,
                               fileName = paste('Cluster', x$cluster, 'Concordant'))
                  }))

allInput = rbind( allInput, 
                  ddply(methyl.genesets, .(cluster),.fun = function(x){
                    data.frame(genesets = x$feature.disconcordant,
                               cluster = x$cluster,
                               parentId = 'syn5838771',
                               executed = thisFile,
                               source.file = methyl.id,
                               fileName = paste('Cluster', x$cluster, 'Discordant'))
                  }))

# Make submissions for splicing
splicing.id = 'syn5834573'
splicing.genesets = downloadFile(splicing.id)

allInput = rbind( allInput, 
                  ddply(splicing.genesets, .(cluster),.fun = function(x){
                    data.frame(genesets = x$feature,
                               cluster = x$cluster,
                               parentId = 'syn5839082',
                               executed = thisFile,
                               source.file = splicing.id,
                               fileName = paste('Cluster', x$cluster))
                  }))

# Make directory and write shell scripts for running these files
system('mkdir sgeEnrichDiffexp')
fp_all = file(paste('sgeEnrichDiffexp/allSubmissions.sh'),'w+')    
cat('#!/bin/bash',file=fp_all,sep='\n')
close(fp_all)

for (i in 1:dim(allInput)[1]) {
  fp = file (paste('/mnt/Github/pcbc_c4_analysis/code/R/sgeEnrichDiffexp/SUB',i,'sh',sep='.'), "w+")
  cat('#!/bin/bash', 
      'sleep 30', 
      paste('Rscript /mnt/Github/pcbc_c4_analysis/code/R/enrichment_diffexp_clusters_SGE.R', 
            paste0('"',paste(t(allInput[i,]), collapse = '" "'),'"')), 
      file = fp,
      sep = '\n')
  close(fp)
  
  fp_all = file(paste('sgeEnrichDiffexp/allSubmissions.sh'),'a+')    
  cat(paste('sh', paste('/mnt/Github/pcbc_c4_analysis/code/R/sgeEnrichDiffexp/SUB',i,'sh',sep='.')),
      file=fp_all,
      sep='\n')
  close(fp_all)
}