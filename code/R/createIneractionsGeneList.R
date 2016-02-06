## It is assumed your working directory is where this file is

# Clear R console screen output
cat("\014")

# Load libraries
library(plyr)
library(dplyr)
library(data.table)
library(stringr)
library(tidyr)
library(reshape2)
library(igraph)

library(tools)
library(ggplot2)
library(RColorBrewer)
library(matrixStats)
library(biomaRt)

library(knitr)
library(knit2synapse)
library(synapseClient)
library(rGithubClient) ## Needs the dev branch

synapseLogin()

# source utility files from ./lib folder
source('./lib/rownameToFirstColumn.R')
source('./lib/get450KProbeMapping.R')

### Create folder in synapse to store results
# Synapse store parameters
parentId = "syn5618211"
ALL_USED_IDs = c()

ActivityName <- 'Convert all interactions to gct format'

ThisFileName <- 'createIneractionsGeneList.R'

# Github link
ThisRepo <- getRepo(repository = "th1vairam/pcbc_c4_analysis", 
                    ref="branch", 
                    refName="diff_exp")

ThisFile <- getPermlink(repository = ThisRepo,
                        repositoryPath=paste0('code/Rmd/', ThisFileName))    

### Functions 
downloadFile <- function(id){
  tmp = fread(synGet(id)@filePath, data.table=F, header=T)
}

### Get expression matrices, metadata and calculate median values at each diff state
#### Get counts matrices
counts_id = c(mrna_count = 'syn5011095', mirna_count = 'syn5014454', methyl_beta = 'syn4487642')
ALL_USED_IDs = c(ALL_USED_IDs, counts_id)

counts = lapply(counts_id, function(id){
  Counts = downloadFile(id)
  
  row.names(Counts) = Counts[,1]
  Counts = as.matrix(Counts[,-(1)])
})

#### Get splicing matrices
splicing_id = 'syn5048714'
ALL_USED_IDs = c(ALL_USED_IDs, splicing_psi = splicing_id)

PSI = downloadFile(splicing_id)
rownames(PSI) = PSI$`Minor-Isoform`
PSI.mapping = PSI[,c('Minor-Isoform', 'Symbol')]
colnames(PSI.mapping) = c('feature','target')
PSI = PSI[,-c(1,227:236)]

counts = c(counts, list(splicing_psi = data.matrix(PSI)))

#### Get metadata matrices
metadata_id = c(mrna_metadata = 'syn5011149', mirna_metadata = 'syn5014460', 
                methyl_metadata = 'syn4487669', splicing_metadata = 'syn5048719')
ALL_USED_IDs = c(ALL_USED_IDs, metadata_id)

metadata = lapply(metadata_id, function(id){
  Covariates = downloadFile(id)
  
  row.names(Covariates) = Covariates[,1]
  Covariates = Covariates[,-(1)]
})

# Find median values of counts, beta and psi 
counts = mapply(function(Covariates, Counts, Assay){
  Covariates = split(Covariates, Covariates$Diffname_short)
  Counts = sapply(Covariates, function(covariates, counts){
    tmp = rowMedians(counts[, rownames(covariates)], na.rm=T)
    names(tmp) = rownames(counts)
    return(tmp)
  }, Counts)
  
  Counts = rownameToFirstColumn(Counts, 'feature') %>%
    dplyr::mutate(Assay = Assay) 
  return(Counts)
}, metadata, counts, c('mrna','mirna','methyl','splicing'), SIMPLIFY = F)
counts = rbindlist(counts, fill = TRUE) %>%
  dplyr::select(feature, Assay, one_of('SC','DE','MESO5','ECTO','MESO15','MESO30','EB'))
WGCNA::collectGarbage()


### Download mappings
#### Get TF-DNA mapping from Enrichr genesets
# Download TFs from Enrichr genesets
load(synGet('syn4867851')@filePath)
ALL_USED_IDs = c(ALL_USED_IDs,genesets = 'syn4867851')

# Get unique TF - gene mapping from three data bases: ChEA, TRANSFAC&JASPAR and ENCODE
ind = grep('HUMAN', names(GeneSets$ChEA))
TFsMapping1 = mapply(function(x, y){
  TFname = str_split(x, '-')[[1]][1]
  return(list(data.frame(feature = TFname, target = y)))
}, names(GeneSets$ChEA)[ind], GeneSets$ChEA[ind]) %>%
  rbindlist %>%
  unique

ind = grep('human', names(GeneSets$TRANSFAC_and_JASPAR_PWMs))
TFsMapping2 = mapply(function(x, y){
  TFname = str_split(x, ' ')[[1]][1]
  return(list(data.frame(feature = TFname, target = y)))
}, names(GeneSets$TRANSFAC_and_JASPAR_PWMs)[ind], GeneSets$TRANSFAC_and_JASPAR_PWMs[ind]) %>%
  rbindlist %>%
  unique

ind = grep('hg19', names(GeneSets$"ENCODE_TF_ChIP-seq_2015"))
TFsMapping3 = mapply(function(x, y){
  TFname = str_split(x, '_')[[1]][1]
  return(list(data.frame(feature = TFname, target = y)))
}, names(GeneSets$"ENCODE_TF_ChIP-seq_2015")[ind], GeneSets$"ENCODE_TF_ChIP-seq_2015"[ind]) %>%
  rbindlist %>%
  unique

TFsMapping = rbindlist(list(TFsMapping1, TFsMapping2, TFsMapping3)) %>% unique           
rm(list = c('GeneSets','TFsMapping1','TFsMapping2','TFsMapping3','ind'))
WGCNA::collectGarbage()

#### Get mRNA-miRNA mapping from synapse (Lorena version)
# Get miRNA mapping files
miRNA.mRNA.id = 'syn3461627'
ALL_USED_IDs <- c(ALL_USED_IDs, mirna.mrna.map = miRNA.mRNA.id)
miRNA.mRNA = fread(synGet(miRNA.mRNA.id)@filePath, data.table=F, header=F)
setnames(miRNA.mRNA, c('V1','V2'), c('hsa_id','ensembl_gene_id'))

# Get human related mapping
Hs = useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org") # use this one when biomart.org is down
Hs = useDataset("hsapiens_gene_ensembl", Hs)
human_ensg2symbol = getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                          filters = "ensembl_gene_id",                         
                          values = unique(miRNA.mRNA$ensembl_gene_id),
                          mart = Hs)

miRNA.mRNA <- left_join(miRNA.mRNA, human_ensg2symbol, by = 'ensembl_gene_id') %>%
  dplyr::rename(feature=hsa_id, target = hgnc_symbol) %>% 
  dplyr::select(feature, target)

#### Get mRNA-miRNA mapping from synapse (exp. validated version)
# Get miRNA mapping files
miRNA.mRNA2.id = 'syn5049680'
ALL_USED_IDs <- c(ALL_USED_IDs, mirna.mrna.map2 = miRNA.mRNA2.id)
miRNA.mRNA2 = fread(synGet(miRNA.mRNA2.id)@filePath, data.table=F, header=T)

miRNA.mRNA2 <- miRNA.mRNA2 %>%
  dplyr::select(one_of(c("miRNA", "Target Gene"))) %>%
  plyr::rename(c("miRNA" = "feature", "Target Gene" = "target")) %>%
  unique

miRNA.mRNA = unique(bind_rows(miRNA.mRNA, miRNA.mRNA2))

#### Get miRNA-methylation mapping from synapse (Lorena version)
# Get miRNA methyl mapping files
methyl.mirna.id = 'syn4895962'
ALL_USED_IDs <- c(ALL_USED_IDs, methyl.mirna = methyl.mirna.id)
methyl.mirna = fread(synGet(methyl.mirna.id)@filePath, data.table=F, header=F)

methyl.mirna <- methyl.mirna %>%  
  plyr::rename(c("V1" = "feature", "V2" = "target")) %>%
  unique

#### Get methylation-mrna mapping
# Get mRNA methyl mapping
methyl.mrna <- get450KProbeMapping(counts$feature[counts$Assay == 'methyl'])$Annotation %>%  
  dplyr::select(methProbeIDs, nearestTx) %>%
  dplyr::rename(feature = methProbeIDs, target = nearestTx) %>%
  unique

#### Combine all interactions together
allInteractions = rbindlist(list(TFsMapping, miRNA.mRNA, miRNA.mRNA2, methyl.mirna, methyl.mrna, PSI.mapping)) %>%
  unique %>%
  dplyr::filter(feature %in% counts$feature, target %in% counts$feature)
allInteractions1 = data.frame(target = allInteractions$feature, feature = allInteractions$target)
allInteractions = rbind(allInteractions, allInteractions1) %>% unique

#### Reshape interactions ad gct files
library(doParallel)
registerDoParallel(cores=16)
newInteractions = ddply(allInteractions[1:1000,], 'feature',.fun=function(x){
  target = paste(unique(x$target),collapse=',')
  return(target)
},.parallel = T);
colnames(newInteractions) = c("feature", "targets")

# Store median counts in synapse
write.table(counts, file='medianCounts.tsv', sep = '\t', quote=F, row.names=F)
obj = File('medianCounts.tsv', name = 'Median Counts, Beta, PSI', parentId = parentId)
obj = synStore(obj, used = as.character(ALL_USED_IDs[c('mrna_count', 'mirna_count', 'methyl_beta', 'splicing_psi',
                                                       "mrna_metadata", "mirna_metadata", "methyl_metadata", "splicing_metadata")]), 
               executed = ThisFile, activityName = ActivityName)

# Store interaction gct files in synapse
write.table(newInteractions, file='interactionsGCT.tsv', sep = '\t', quote=F, row.names=F)
obj = File('interactionsGCT.tsv', name = 'All Interactions (gct format)', parentId = parentId)
obj = synStore(obj, used = as.character(ALL_USED_IDs[c('genesets', 'mirna.mrna.map' ,'mirna.mrna.map2', 'methyl.mirna', 
                                                       'methyl_beta', 'splicing_psi')]), 
               executed = ThisFile, activityName = ActivityName)