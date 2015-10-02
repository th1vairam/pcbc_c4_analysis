#### Code to merge enrichment results ####

# Clear R console screen output
cat("\014")

# Load required libraries
library(data.table)
library(dplyr)
library(tools)
library(synapseClient)

## Needs the dev branch
library(rGithubClient)

synapseLogin()

# Get github links for provenance
thisFileName <- 'mergeEnrichmentResults_mRNA_SGE.R'
thisRepo <- getRepo(repository = "th1vairam/pcbc_c4_analysis", ref="branch", refName='enrich')
thisFile <- getPermlink(repository = thisRepo, repositoryPath=paste0('code/R/', thisFileName))

# Get all the enrichment files
All.Files = synQuery('select * from file where parentId == "syn4905034"')

# Get all the comparison names from synapse table
compNames <- synTableQuery("SELECT * FROM syn4483642")@values
compNames = rbind(compNames %>% mutate(direction = 'up mRNA', comparisonName = str_c(comparison,direction,sep='__')),
                  compNames %>% mutate(direction = 'down mRNA', comparisonName = str_c(comparison,direction,sep='__')))
compNames = split(compNames, compNames$class, drop = F)
compNames = lapply(compNames, function(x){ x = split(x, x$dataRestrictionShort, drop = F)})

for (class in names(compNames)){
  for(dataRestrictionShort in names(compNames[[class]])){
    file.ids = All.Files$file.id[All.Files$file.name %in% compNames[[class]][[dataRestrictionShort]]$comparisonName]
    if (length(file.ids) != 0){
      results = lapply(file.ids, function(id){x = fread(synGet(id)@filePath, data.table=F, header=T)})
      results = ldply(results)
      write.table(results, file = paste(dataRestrictionShort,class,'mRNA',sep="."), sep="\t")
      file = File(path = paste(dataRestrictionShort,class,'mRNA',sep="."), 
                  name = paste(paste(dataRestrictionShort,class,sep="."),'mRNA'),
                  parentId = "syn4896411")
      file = synStore(file, used = file.ids, executed = thisFile, activityName = "Merge enrichmnet results")
    }
  }
}