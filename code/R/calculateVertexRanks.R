## It is assumed your working directory is where this file is
## Function to rank the vertices of networks

# Clear R console screen output
cat("\014")

#### Load libraries
loadLibs <- function(){
  library(CovariateAnalysis)
  library(data.table)
  library(tidyr)
  library(plyr)
  library(dplyr)
  library(stringr)

  library(synapseClient)
  library(knitr)
  library(githubr)

  library(igraph)

  library(doParallel)
  library(foreach)
}
loadLibs()

cl <- makeCluster(2)
registerDoParallel(cl)

synapseLogin()

#### Source required R files
source('./lib/rankVertices.R')

#### Get differential expression results
diffexp = downloadFile('syn6039731')
ALL_USED_IDs = 'syn6039731'

#### Get all active interactions
all.interactions = downloadFile('syn6039630') %>%
  unite(Comparison, from.state, to.state, sep = '_vs_')
ALL_USED_IDs = c(ALL_USED_IDs, 'syn6039630')

#### Get regulatory influence scores for each node
rankList = all.interactions %>%
  ddply(.(Comparison), .fun = function(edges, diffexp, rankVertices, loadLibs){
    loadLibs()
    vertices = filter(diffexp,
                      feature %in% unique(c(edges$feature, edges$target)),
                      Comparison == unique(edges$Comparison)) %>%
      data.frame
    rownames(vertices) = vertices$feature
    
    edges$target[is.na(edges$target)] = 'NULL'
    g = igraph::graph_from_edgelist(edges %>%
                                      dplyr::select(feature, target) %>% 
                                      as.matrix,
                                    directed = T)
    E(g)$weight = abs(edges$coexpression)
    g = g - V(g)[V(g)$name == 'NULL']
    
    V(g)$fc = abs(vertices[V(g)$name, 'logFC'])
    V(g)$fdr = vertices[V(g)$name, 'adj.P.value']
    
    rk = rankVertices(g, 3) %>%
      rownameToFirstColumn('feature') %>%
      plyr::rename(c('DF' = 'rank')) %>%
      dplyr::mutate(Comparison = unique(edges$Comparison))
  }, diffexp, rankVertices, loadLibs, .progress = 'text', .parallel = T)

#### Store results in synapse
# Create folder to store results in synapse
ActivityName <- 'Rank nodes based on differential expresison and network structure'

ThisFileName1 <- 'calculateVertexRanks.R'
ThisFileName2 <- 'lib/rankVertices.R'

# Github link
ThisRepo <- getRepo(repository = "th1vairam/pcbc_c4_analysis", 
                    ref="branch", 
                    refName='discordant_anal')

ThisFile1 <- getPermlink(repository = ThisRepo,
                        repositoryPath=paste0('code/R/', ThisFileName1))
ThisFile2 <- getPermlink(repository = ThisRepo,
                         repositoryPath=paste0('code/R/', ThisFileName2))

# Write rank list to synapse
write.table(rankList, file = 'allNodeRankings.tsv', sep = '\t', row.names=F, quote=F)
obj = File('allNodeRankings.tsv', name = 'All Node Rankings', parentId = 'syn5996097')
obj = synStore(obj, used = as.character(ALL_USED_IDs), executed = list(ThisFile1,ThisFile2), activityName = ActivityName)
