#!/usr/bin/env Rscript

# Function to perform enrichment analysis of differential expression modules
# Get arguments from command line
args = commandArgs(TRUE)

# Clear R console screen output
cat("\014")

# Clear R workspace
setwd('/mnt/Github/pcbc_c4_analysis/code/R')
.libPaths('/mnt/mylibs/')
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

# Login to synapse
synapseLogin()
############################################################################################################

############################################################################################################
#### Github commit ####
# Get github links for provenance
thisFileName <- 'enrichment_diffexp_clusters_SGE.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/pcbc_c4_analysis", 
                    ref="branch", 
                    refName='enrich')

thisFile <- getPermlink(repository = thisRepo,
                        repositoryPath=paste0('code/R/', thisFileName))

# Synapse specific parameters
activityName = 'Enrichment Analysis'
activityDescription = 'Enrichment analysis based on enrichr genesets'
############################################################################################################

############################################################################################################
#### Get gene sets ####
# Download enrichr gene sets from synapse
GL_OBJ = synGet('syn4867851')
ALL_USED_IDs = GL_OBJ$properties$id
load(GL_OBJ@filePath) # This RData file will load a list names GeneSets

gsets = c("BioCarta_2015", "Cancer_Cell_Line_Encyclopedia", "ChEA", "Chromosome_Location", "Cross_Species_Phenotype",
          "Disease_Signatures_from_GEO_down", "Disease_Signatures_from_GEO_up", "Drug_Perturbations_from_GEO",
          "ENCODE_Histone_Modifications_2015", "ENCODE_TF_ChIP-seq_2015", "ESCAPE", "Epigenomics_Roadmap_HM_ChIP-seq",
          "GO_Biological_Process", "GO_Cellular_Component", "GO_Molecular_Function", "GeneSigDB", "Genome_Browser_PWMs.1",
          "HomoloGene", "HumanCyc", "Human_Gene_Atlas", "Human_Phenotype_Ontology", "KEA", "KEGG_2015", "MSigDB_Oncogenic_Signatures",
          "Mouse_Gene_Atlas", "OMIM_Disease","OMIM_Expanded", "PPI_Hub_Proteins", "Panther", "Reactome_2015", 
          "TF-LOF_Expression_from_GEO", "TRANSFAC_and_JASPAR_PWMs", "TargetScan_microRNA", "Transcription_Factor_PPIs","WikiPathways_2015")
GeneSets = GeneSets[gsets]
############################################################################################################

############################################################################################################
#### Get differentialy expressed gene sets ####
ALL_USED_IDs = c(ALL_USED_IDs, args[1])
diffexp.genesets = downloadFile(args[1])
cluster = args[3]

genesToTest = str_split(diffexp.genesets$feature[diffexp.genesets$cluster == cluster], pattern = ',') %>%
  unlist
############################################################################################################

############################################################################################################
#### Background gene list ####
counts.id = 'syn5011095'
ALL_USED_IDs = c(ALL_USED_IDs, counts.id)
counts = downloadFile(counts.id)
backGroundGenes = unique(counts$GeneName)
############################################################################################################

############################################################################################################
#### Filter gene list ####
GeneSets = filterGeneSets(GeneSets, backGroundGenes, minSize = 10, maxSize = 5000)
############################################################################################################

############################################################################################################
#### Perform enrichment analysis ####
enrichResults = lapply(GeneSets,
                       function(geneSetsList, genesToTest, genesInBackground){
                         tmp = lapply(geneSetsList, 
                                      function(genesInGeneSet, genesInSignificantSet, genesInBackground){
                                        fisherEnrichment(genesInSignificantSet, genesInGeneSet, genesInBackground)
                                        }, genesToTest, genesInBackground) %>%
                           rbindlist(idcol = 'SetName', use.names=T, fill=T)
                         }, genesToTest, backGroundGenes) %>%
  rbindlist(idcol = 'CategoryName', use.names=T, fill=T)

enrichResults$FDR = p.adjust(enrichResults$pval,'fdr')
enrichResults$adj.pval = p.adjust(enrichResults$pval,'bonferroni')
write.table(enrichResults, file = paste('Cluster',cluster,'tsv', sep='.'), sep='\t', row.names=F)
############################################################################################################

############################################################################################################
#### Write to synapse ####
# Write results to synapse
ENR_OBJ = File(paste('Cluster',cluster,'tsv', sep='.'),
               name = paste('Cluster',cluster), 
               parentId = args[2])

annotations(ENR_OBJ) = list(algo = 'Fisher')

ENR_OBJ = synStore(ENR_OBJ, 
                   executed = thisFile,
                   used = ALL_USED_IDs,
                   activityName = activityName,
                   activityDescription = activityDescription)
############################################################################################################