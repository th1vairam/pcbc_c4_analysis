#!/usr/bin/env Rscript
# Calculate correlation between mrna, mirna, and methylation expressions at any 2 differentiation states

## It is assumed your working directory is where this file is
setwd('/mnt/Github/pcbc_c4_analysis/code/R/')

# Get arguments from command line
args = commandArgs(TRUE)

# Clear R console screen output
cat("\014")

# Load libraries
library(plyr)
library(dplyr)
library(data.table)
library(stringr)
library(tidyr)
library(reshape2)

library(matrixStats)
library(biomaRt)

library(knitr)
library(knit2synapse)
library(synapseClient)
library(rGithubClient) ## Needs the dev branch

synapseLogin()

# Source needed files from lib folder
source('./lib/get450KProbeMapping.R')
source('./lib/rownameToFirstColumn.R')

# Synapse store parameters
parentId = "syn5194922"
ALL_USED_IDs = c()

# Create folder to store results in synapse
CODE <- Folder(name = 'Coexpression Between Different Assays At Any Two Differentiation Stage', parentId = parentId)
CODE <- synStore(CODE)

ActivityName <- 'Calculate co-expression between assays'

ThisFileName <- 'CalculateCoExp_Any2DiffState.R'

# Github link
ThisRepo <- getRepo(repository = "th1vairam/pcbc_c4_analysis", 
                    ref="branch", 
                    refName='discordant_anal')

ThisFile <- getPermlink(repository = ThisRepo,
                        repositoryPath=paste0('code/R/', ThisFileName))    

#### Block for all related functions ####
# Utility function to download files from synapse
downloadFile <- function(id){
  tmp = fread(synGet(id)@filePath, data.table=F, header=T)
}

# Function to calculate co-expression
calculateCor <- function(metaData, omics1Mat, omics2Mat, Interactions){
  compName = paste(unique(metaData$Diffname_short), collapse = '_vs_')
  
  ## For group1
  omics1Mat = omics1Mat[,metaData$biologicalSampleName]
  omics2Mat = omics2Mat[,metaData$biologicalSampleName]
  
  # Calcualte correlation between features at group1
  correlation = apply(Interactions, 1, function(x, omics1Mat, omics2Mat){
    WGCNA::bicor(t(omics1Mat[as.character(x[1]),]), t(omics2Mat[as.character(x[2]),]), use = "pairwise.complete.obs")  
  }, omics1Mat, omics2Mat)
  
  correlation = Interactions %>%
    dplyr::mutate(correlation = correlation)
  
  writeLines(paste('Completed',compName))
  return(correlation)
}

#### Get expression matrices and metadata ####
# Get counts matrices
counts_id = c(mrna = 'syn5011095', mirna = 'syn5014454', methyl = 'syn4487642')
ALL_USED_IDs = c(ALL_USED_IDs, counts_id)

counts = lapply(counts_id, function(id){
  Counts = downloadFile(id)
  
  row.names(Counts) = Counts[,1]
  # colnames(Counts) = gsub('-','.',colnames(Counts))
  Counts = as.matrix(Counts[,-(1)])
})

# Get splicing matrices
spliceJunction_id = 'syn5048714'
ALL_USED_IDs = c(ALL_USED_IDs, spliceJunction = spliceJunction_id)

PSI = downloadFile(spliceJunction_id)
rownames(PSI) = PSI$`Minor-Isoform`
PSI = PSI[,-c(1,227:236)]

counts = c(counts, list(spliceJunction = data.matrix(PSI)))

#### Get metadata matrices ####
# Get metadata matrices
metadata_id = c(mrna = 'syn3156503', mirna = 'syn3219876', 
                methyl = 'syn3156828', spliceJunction = 'syn3156503')
ALL_USED_IDs = c(ALL_USED_IDs, metadata_id)

# Get metadata from synapse
metadata = lapply(metadata_id, function(id){
  # Get metadata from synapse tables
  MetaData <- synTableQuery(sprintf('SELECT * FROM %s', id))@values
  MetaData[MetaData == 'N/A'] = NA
  
  # Replace all special characters with blank
  myFix <- function(x) str_replace_all(x, '[^[:alnum:]]', '')
  MetaData <- MetaData %>%
    dplyr::mutate_each(funs(myFix), -UID, -C4_Cell_Line_ID, -biologicalSampleName,
                       -public, -pass_qc, -exclude) # fix them but don't touch some columns
  
  # MetaData$UID = gsub('-','.',MetaData$UID)
  rownames(MetaData) = MetaData$UID
  return(MetaData)
})

# Arrange metadata rows wrt counts columns
metadata = mapply(function(Counts, MetaData){
  ind = intersect(colnames(Counts),rownames(MetaData))
  MetaData = MetaData[ind,]
  MetaData = MetaData[colnames(Counts),]
}, counts, metadata, SIMPLIFY = F)

#### Download all interactions from synapse ####
ALL_USED_IDs = c(ALL_USED_IDs, all.interactions = 'syn5643683')
allInteractions = downloadFile('syn5643683') %>%
  dplyr::mutate(Assay = paste(feature.assay,target.assay,sep = '_'))
allInteractions = split(allInteractions, allInteractions$Assay)

#### Calculate co-expression in all given category of interactions at any given 2 diffstates ####
# combine samples to unique biological sample name (take median values for replicates)
interactions = allInteractions[[args[1]]]
Counts = counts
MetaData = metadata

feature.assay = unique(interactions$feature.assay)
target.assay = unique(interactions$target.assay)
  
# Match up biological sample between assays
biosampleInBoth <- intersect(MetaData[[feature.assay]]$biologicalSampleName, 
                             MetaData[[target.assay]]$biologicalSampleName)

# Filter metadata
feature.metadata <- MetaData[[feature.assay]] %>%
  filter(biologicalSampleName %in% biosampleInBoth)

target.metadata <- MetaData[[target.assay]] %>%
  filter(biologicalSampleName %in% biosampleInBoth)
  
# Take the median across multiple biological samples per feature
feature.median <- Counts[[feature.assay]] %>%
  rownameToFirstColumn('feature') %>%
  dplyr::select(feature, one_of(feature.metadata$UID)) %>%
  dplyr::filter(feature %in% unique(interactions$feature)) %>%
  melt %>%
  dplyr::rename(UID = variable, expression = value) %>%
  left_join(feature.metadata %>% dplyr::select(UID, biologicalSampleName)) %>%
  group_by(feature, biologicalSampleName) %>% 
  summarize(median_expression=median(expression)) %>% 
  dcast(feature ~ biologicalSampleName)

rownames(feature.median) = feature.median$feature
feature.median$feature = NULL
  
target.median <- Counts[[target.assay]] %>%
  rownameToFirstColumn('feature') %>%
  dplyr::select(feature, one_of(target.metadata$UID)) %>%
  dplyr::filter(feature %in% unique(interactions$target)) %>%
  melt %>%
  dplyr::rename(UID = variable, expression = value) %>%
  left_join(target.metadata %>% dplyr::select(UID, biologicalSampleName)) %>%
  group_by(feature, biologicalSampleName) %>% 
  summarize(median_expression=median(expression)) %>% 
  dcast(feature ~ biologicalSampleName)

rownames(target.median) = target.median$feature
target.median$feature = NULL
  
# Get unique metadata
unique.metadata <- feature.metadata %>%
  dplyr::select(biologicalSampleName, Diffname_short, Originating_Lab, Cell_Type,
                Cell_Line_Type, Cell_Line_of_Origin, Tissue_of_Origin,
                Reprogramming_Gene_Combination, Culture_Conditions,
                Cell_Type_of_Origin_Level2, Reprogramming_Vector_Type_Level2) %>%
  dplyr::filter(biologicalSampleName %in% biosampleInBoth) %>%
  dplyr::mutate(Diffname_short = gsub('-','',Diffname_short)) %>% 
  unique
  
feature.median = feature.median[,unique.metadata$biologicalSampleName]
target.median = target.median[,unique.metadata$biologicalSampleName]
  
# Split metadata in terms of diffstate
unique.metadata <- split(unique.metadata, unique.metadata$Diffname_short)
  
# Combine metadata for all combinations of any 2 diffstate
all.combination = combn(names(unique.metadata),2)
unique.metadata = apply(all.combination, 2, function(x, MetaData){
  y = MetaData[x]
  y = rbindlist(y)
  return(y)
}, unique.metadata)
names(unique.metadata) = apply(all.combination,2,paste, collapse = '_vs_')
  
# Calculate co-expression
correlation = calculateCor(unique.metadata[[args[2]]],
                           feature.median, target.median, 
                           interactions) %>%
  plyr::rename(c('correlation' = args[2]))

# Store coexpression matrix in synapse
write.table(correlation, file = paste('correaltion',args[1],args[2],'tsv',sep='.'), sep = '\t', quote=F, row.names=F)
obj = File(paste('correaltion',args[1],args[2],'tsv',sep='.'), 
           name = paste('correaltion',args[1],args[2]), parentId = CODE$properties$id)
annotations(obj) = list(interactionAssay = args[1], diffStateComparison = args[2])
obj = synStore(obj, used = as.character(ALL_USED_IDs), executed = ThisFile, activityName = ActivityName)