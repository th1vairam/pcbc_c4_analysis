#!/usr/bin/env Rscript

# Function to consolidate all enrichment results 
# Clear R console screen output
cat("\014")

# Clear R workspace
setwd('/mnt/Github/pcbc_c4_analysis/code/R')
############################################################################################################

############################################################################################################
#### Libraries ####
library(synapseClient)
library(knitr)
library(githubr)
library(CovariateAnalysis)
library(data.table)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)

synapseLogin()
############################################################################################################

############################################################################################################
# Consolidate all enrichment results based on category
# Get differential expression (ALL)
diffexp.id = 'syn5706668'
diffexp = downloadFile(diffexp.id) %>%
  tbl_df

# Get enrichment analysis results
fileNames = synQuery(paste0('select * from file where parentId == "syn5762641"')) %>%
  dplyr::rename(setname = file.name) %>%
  left_join(diffexp)

# Get enrichment of diffstate specific up and down regulated genes (Only molecular pathways)
orderedDiffState = c('SC', 'DE', 'MESO5', 'ECTO', 'MESO15', 'MESO30', 'EB')
for (reg in c("none")){
  for (diffState in orderedDiffState){
    for (dir in c('UP', 'DOWN')){
      tmp = filter(fileNames, fromstate == diffState, assay == reg, direction == dir)
      enrich = mapply(function(x,y){
        results = downloadFile(x) %>%
          filter(CategoryName %in% c("BioCarta_2015", "GO_Biological_Process", "KEGG_2015", "Panther", 
                                     "Reactome_2015", "WikiPathways_2015")) %>%
          dplyr::rename(GeneSetName = SetName)
        
        # Remove mouse related gene sets
        ind1 = setdiff(1:dim(results)[1], grep('Mus musculus', results$GeneSetName))
        ind2 = setdiff(1:dim(results)[1], grep('mm9', results$GeneSetName))
        results = results[intersect(ind1,ind2),]  
        
        results = results %>%
          dplyr::mutate(FDR = p.adjust(pval)) %>%
          filter(FDR <= 0.01) %>%
          dplyr::select(CategoryName, GeneSetName, Odds.Ratio) %>%
          arrange(desc(Odds.Ratio)) %>%
          dplyr::slice(1:100) %>% 
          plyr::rename(c('Odds.Ratio' = paste(y, x, sep = '.')))
        
        }, tmp$file.id, tmp$tostate, SIMPLIFY = F) %>%
        join_all(type = 'full', match = 'all') 
      enrich[is.na(enrich)] = 0
      
      write.xlsx(enrich, file = paste0('EnrichmentResults_regulatedBy_',reg,'.xlsx'), 
                 sheetName = paste(diffState, dir, 'FDR_0.01', 'TOP_100_PATHWAYS'),
                 row.names =F, append = T)
    }
  }
  obj = File(paste0('EnrichmentResults_regulatedBy_',reg,'.xlsx'), name = paste('Enrichment Results Regulated By',reg), parentId = 'syn5752514')
  obj = synStore(obj, used = c('syn5706668', 'syn5762641'), activityName = 'Consolidate Enrichment Results')
}