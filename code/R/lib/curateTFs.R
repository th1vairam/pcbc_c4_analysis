## It is assumed your working directory is where this file is

# Clear R console screen output
cat("\014")

# Load libraries
library(CovariateAnalysis)
library(data.table)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)

library(synapseClient)
library(knitr)
library(githubr)

synapseLogin()

ALL_USED_IDs = c(ALL_USED_IDs, TFdb = 'syn5672453', map = 'syn5986999', Enrichr = 'syn4867851')

# Get TF checkpoints db from synapse
TFs = read.table(synGet('syn5672453')@filePath, header=T, sep= '\t', fill = TRUE) %>%
  dplyr::filter(entrez_human != 0) %>%
  dplyr::select(entrez_human) %>%
  setnames('entrez_human','entrezgene') %>%
  dplyr::mutate(Assay = 'TF')

# Get HGNC mapping from synapse
human_entrz2symbol = read.table(synGet('syn5986999')@filePath, fill=NA, quote='', header=T, sep = '\t') %>%
  dplyr::select(Approved.Symbol, Entrez.Gene.ID) %>%
  dplyr::rename(hgnc_symbol = Approved.Symbol, entrezgene = Entrez.Gene.ID)

# Convert entrex gene ids to gene symbols
TFs = left_join(TFs, human_entrz2symbol) %>%
  dplyr::select(hgnc_symbol, Assay) %>%
  setnames('hgnc_symbol','feature') %>%
  filter(!is.na(feature)) %>%
  unique

# Get TFs from enrichr
load(synGet('syn4867851')@filePath) # Will load robj called GeneSets
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

TFsMapping = rbindlist(list(TFsMapping1, TFsMapping2, TFsMapping3)) %>% 
  unique %>%   
  dplyr::mutate(feature.assay = 'mrna', target.assay = 'mrna', 
                dataSource = 'ENCODE_TF_ChIP-seq_2015, TRANSFAC_and_JASPAR_PWMs, ChEA')
rm(list= 'GeneSets')
gc()

TFs = rbind(TFs,
            data.frame(feature = setdiff(TFsMapping$feature %>% unique, TFs$feature),
                       Assay = 'TF'))

#### Store results in synapse
# Create folder to store results in synapse
ActivityName <- 'Consolidate TFs from enrichr and checkpointsTFDB'

ThisFileName <- 'curateTFs.R'

# Github link
ThisRepo <- getRepo(repository = "th1vairam/pcbc_c4_analysis", 
                    ref="branch", 
                    refName='discordant_anal')

ThisFile <- getPermlink(repository = ThisRepo,
                       repositoryPath=paste0('code/R/lib/', ThisFileName))

# Write rank list to synapse
write.table(TFs, file = 'TF.tsv', sep = '\t', row.names=F, quote=F)
obj = File('TF.tsv', name = 'Transcription Factors', parentId = 'syn4598159')
obj = synStore(obj, used = as.character(ALL_USED_IDs), executed = ThisFile, activityName = ActivityName)