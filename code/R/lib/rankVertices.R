rankVertices <- function(g, order){
  
  Gx = igraph::V(g)$fc * -log10(igraph::V(g)$fdr)
  names(Gx) = igraph::V(g)$name
  
  neighbor.graph = igraph::make_ego_graph(g, order = order, nodes = igraph::V(g), mode = "out")
  names(neighbor.graph) = igraph::V(g)$name
  
  D = igraph::distances(g)
  
  S = igraph::strength(g, mode = "out")
  
  Nx = mapply(function(x, y){
    SP = igraph::shortest_paths(y, from = x, to = V(y), mode = "out", predecessors = T)
    
    Nx = (Gx[igraph::V(y)$name]*(1/D[x,igraph::V(y)$name])*(1/S[igraph::as_ids(SP$predecessors)]))
    names(Nx) = V(y)$name
    return(sum(Nx[names(Nx) != x]))
  }, V(g)$name[S!=0], neighbor.graph[S!=0])
  Nx = c(Nx, S[S == 0])
  
  return(sort(rank(Gx*Nx[names(Gx)])))
}