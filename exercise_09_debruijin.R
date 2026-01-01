library(Biostrings)

# 2.1 Vytvoření uzlů a hran
DeBruijnGraph <- function(kmers){
  nodes <- unique(unlist(lapply(kmers, function(x) c(substr(x, 1, nchar(x)-1),
                                                     substr(x, 2, nchar(x))))))
  
  edges <- lapply(kmers, function(kmer){
    from <- substr(kmer, 1, nchar(kmer)-1)
    to   <- substr(kmer, 2, nchar(kmer))
    c(from, to)
  })
  
  list(nodes = nodes, edges = edges)
}

CanHaveEulerianPath <- function(graph){
  # graph = list(nodes = ..., edges = ...)
  
  # 1. vytvoříme adjacency list
  out_deg <- setNames(rep(0, length(graph$nodes)), graph$nodes)
  in_deg  <- setNames(rep(0, length(graph$nodes)), graph$nodes)
  
  for(edge in graph$edges){
    from <- edge[1]
    to   <- edge[2]
    out_deg[from] <- out_deg[from] + 1
    in_deg[to]   <- in_deg[to] + 1
  }
  
  # 2. spočítáme rozdíly
  start_nodes <- sum((out_deg - in_deg) == 1)
  end_nodes   <- sum((in_deg - out_deg) == 1)
  
  # 3. všechny ostatní uzly musí mít out_deg == in_deg
  balanced_nodes <- all((out_deg - in_deg)[!(out_deg - in_deg %in% c(1,-1))] == 0)
  
  if ((start_nodes == 1 && end_nodes == 1 && balanced_nodes) ||
      (start_nodes == 0 && end_nodes == 0 && balanced_nodes)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

kmers <- c("ATT","TTC","TCG","CGA","GAT")
dbg <- DeBruijnGraph(kmers)
CanHaveEulerianPath(dbg)



# 2.2 Eulerian path
# jednoduchá verze: Hierholzerův algoritmus
FindEulerianPath <- function(graph){
  # vytvoříme adjacency list
  adj <- list()
  for (node in graph$nodes) adj[[node]] <- character(0)
  for (e in graph$edges){
    adj[[e[1]]] <- c(adj[[e[1]]], e[2])
  }
  
  path <- character(0)
  stack <- character(0)
  curr <- graph$edges[[1]][1]  # startovní uzel
  
  while(length(stack) > 0 || length(adj[[curr]]) > 0){
    if(length(adj[[curr]]) == 0){
      path <- c(path, curr)
      curr <- tail(stack,1)
      stack <- head(stack,-1)
    } else {
      stack <- c(stack, curr)
      next_node <- adj[[curr]][1]
      adj[[curr]] <- adj[[curr]][-1]
      curr <- next_node
    }
  }
  path <- c(path, curr)
  rev(path)
}

# 2.3 Assembly z Eulerian path
AssembleFromPath <- function(path){
  seq <- path[1]
  for(i in 2:length(path)){
    seq <- paste0(seq, substr(path[i], nchar(path[i]), nchar(path[i])))
  }
  seq
}

kmers <- c("ATT","TTC","TCG","CGA","GAT")
dbg <- DeBruijnGraph(kmers)
path <- FindEulerianPath(dbg)
assembly <- AssembleFromPath(path)

assembly

