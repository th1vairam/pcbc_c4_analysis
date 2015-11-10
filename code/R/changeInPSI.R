# title: Get mean change in PSI for each differential comparisons
# author: "Thanneer Perumal"

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
library(biomaRt)
library(tools)

library(knitr)
library(knit2synapse)
library(synapseClient)
library(rGithubClient) ## Needs the dev branch

synapseLogin()

### Download differential expression results
# Get all comparison names from synapse
compNames <- synTableQuery("SELECT * FROM syn4483642")@values
ALL_USED_IDs = "syn4483642"

# Add direction to comparison names
compNames = rbind(compNames %>%
                    mutate(direction = 'up',
                           comparisonName = str_c(comparison,direction,sep='__')),
                  compNames %>%
                    mutate(direction = 'down',
                           comparisonName = str_c(comparison,direction,sep='__')))

# Add change in PSI values for each comparison
PSI_ID = 'syn5048714'
METADATA_ID = 'syn3156503'
ALL_USED_IDs = c(ALL_USED_IDs, PSI_ID, METADATA_ID)

# Download PSI values
PSI_OBJ = synGet(PSI_ID)
PSI = fread(getFileLocation(PSI_OBJ), data.table=F, header=T)
rownames(PSI) = PSI$"Minor-Isoform"

# Seperate PSI and PSI annotation
PSI.ANNOT = dplyr::select(PSI,  one_of(c("Symbol", "Description", "Minor-Isoform", "Major Isoform", "AltExons", "PME",
                                         "dPSI", "rho", "Max Inclusion PSI", "Coordinates", "feature")))
PSI = PSI[, setdiff(colnames(PSI), colnames(PSI.ANNOT))]
colnames(PSI) = gsub('.bed', '', colnames(PSI))

# Specify factor and continuous covarites pool (adjusted covariates will only be the subset of these covariates)
FactorCovariates = c('Diffname_short', 'run', 'lane', 'Cell_Line_Type', 'Cell_Line_of_Origin', 'Tissue_of_Origin', 'Reprogramming_Gene_Combination', 'Culture_Conditions', 'Donor_Life_Stage', 'Race', 'Ethnicity' , 'Gender', 'Disease', 'Originating_Lab', 'Donor_ID', 'Cell_Type_of_Origin_Level2', 'Reprogramming_Vector_Type')
ContCovariates = c('PassageAtThaw', 'PassageAtHarvest')

# Get metadata
METADATA_OBJ = synTableQuery(paste('SELECT * FROM',METADATA_ID,sep=' '))
ALL_USED_IDs[length(ALL_USED_IDs)+1] = METADATA_OBJ@schema
METADATA = METADATA_OBJ@values

# Preprocess metadata
METADATA[METADATA == 'N/A'] = NA

# Replace all special characters with blank
myFix <- function(x) str_replace_all(x, '[^[:alnum:]]', '')
METADATA <- METADATA %>%
  dplyr::mutate_each(funs(myFix), -UID, -C4_Cell_Line_ID, -biologicalSampleName,
                     -public, -pass_qc, -exclude) # fix them but don't touch some columns

# Set rownames
rownames(METADATA) = METADATA$UID
METADATA = METADATA[colnames(PSI),]

# Find inter relation between factor covariates
COVARIATES = METADATA[,c(FactorCovariates,ContCovariates)]

# Convert factor covariates to factors
COVARIATES[,FactorCovariates] = lapply(COVARIATES[,FactorCovariates], factor)
COVARIATES[,ContCovariates] = lapply(COVARIATES[,ContCovariates], as.numeric)

# Get differential splicing results (this is used to filter only those comparisons which are of use)
splicingIds = c(splicingId.all = "syn5049321", splicingId.DE = "syn5065271", splicingId.EB = "syn5065380", 
                splicingId.ECTO = "syn5065336", splicingId.MESO5 = "syn5065297", splicingId.SC = "syn5065245")
ALL_USED_IDs = c(ALL_USED_IDs, as.character(splicingIds))
splicing.diffExp = lapply(splicingIds, function(id){
  tmp = fread(synGet(id)@filePath, data.table=F, header=T)
})
diffExp.comparison = rbindlist(splicing.diffExp) %>%
  dplyr::select(Comparison) %>% unlist %>% unique

compNames.new = filter(compNames, comparisonName %in% diffExp.comparison, !(class %in% c("Cell_Type_of_Origin")))

WGCNA::collectGarbage()
changePSI = apply(compNames.new, 1, function(x, covariates, psiFraction){  
  ind1 = which(covariates[,as.character(x['class'])] == as.character(x["variable1Short"]))
  ind2 = which(covariates[,as.character(x['class'])] == as.character(x["variable2Short"]))
  
  if (length(ind1) > 0 && length(ind2) > 0 ){
    var1Mean = rowMeans(psiFraction[,covariates$UID[ind1], drop=F], na.rm=T)
    var2Mean = rowMeans(psiFraction[,covariates$UID[ind2], drop=F], na.rm=T)
    changePSI = var1Mean - var2Mean
    
    changePSI = rownameToFirstColumn(changePSI, 'JunctionIDs')    
  } else {
    changePSI = data.frame(methProbeIDs =NA, DF=NA)    
  }
  colnames(changePSI)[2] = x["comparisonName"]
  return(changePSI)
}, COVARIATES, PSI)
changePSI = plyr::join_all(changePSI)

### Synapse store
activityName = "Calculate change in PSI for all differential comparisons"

parentId = 'syn4991628'

thisFileName <- 'changePSI.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/pcbc_c4_analysis", 
                    ref="branch", 
                    refName='splicing')

thisFile <- getPermlink(repository = thisRepo,
                        repositoryPath=paste0('code/R/', thisFileName))

# Write large matrix to RData file
save(list = "changePSI" , file = "changeInPSI.RData")

# Write file to synapse
file.obj = File("changeInPSI.RData", name = "Change in PSI", parentId = parentId)
file.obj = synStore(file.obj, activityName = activityName, used = as.character(ALL_USED_IDs), executed = thisFile)