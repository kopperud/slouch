parent <- function(phy, x){
  m <- which(phy$edge[, 2] == x)
  return(phy$edge[m, 1])
}

lineage.nodes <- function(phy, x){
  k <- x
  N <- length(phy$tip.label)
  while(x != N + 1){
    k <- c(k, parent(phy, x))
    x <- tail(k, n = 1)
  }
  return(k)
}

lineage.constructor <- function(phy, e, anc_maps, regimes, ace){
  nodes <- lineage.nodes(phy, e)
  min_age <- min(node.depth.edgelength(phy)[nodes])
  
  if(anc_maps == "regimes"){
    which.regimes <- lapply(levels(regimes), function(x) {res <- match(regimes[nodes], x); res[is.na(res)] <- 0; return(res) })
    times <-  ape::node.depth.edgelength(phy)[nodes]
    timeflip <- times[1] - times ## Time from tip to node(s)
  }else if(anc_maps == "ace"){
    stopifnot(!is.null(ace))
    tip_regimes <- regimes[1:length(phy$tip.label)]
    which.tips <-  sapply(levels(regimes), function(x) {res <- match(tip_regimes, x); res[is.na(res)] <- 0; return(res) })
    which.internal <- ace$lik.anc
    
    which.regimes <- rbind(which.tips, which.internal)[nodes,]
    which.regimes <- lapply(1:ncol(which.regimes), function(i) which.regimes[,i])
    
    times <-  ape::node.depth.edgelength(phy)[nodes]
    timeflip <- times[1] - times ## Time from tip to node(s)
  }else if(anc_maps == "simmap"){
    ## Simmap splits up each edge into sub-edges, depending on the split. So, we use edges instead of nodes, and introduce sub-edges
    edge_is <- which(phy$edge[,2] %in% nodes)
    subedges <- unlist(lapply(edge_is, function(i) phy$maps[[i]]))
    simmap_regimes <- rev(names(subedges))
    
    which.regimes <- lapply(levels(regimes), function(x) {res <- match(simmap_regimes, x); res[is.na(res)] <- 0; return(res)})
    # Problem. simmap does not provide root estimate. Assuming root estimate is equal to the oldest branch estimate
    root <- lapply(which.regimes, function(e) tail(e, n= 1))
    which.regimes <- lapply(1:3, function(x) c(which.regimes[[x]], root[[x]]))
  
    timeflip <- cumsum(c(min_age, unname(subedges)))
    times <- rev(timeflip)
  }
  
  names(which.regimes) <- levels(regimes)
  
  t_end <- tail(timeflip, n = -1) ## Time from tip to end of segment(s)
  t_beginning <- head(timeflip, n = -1) ## Time from tip to beginning of segment(s)
  regime_time <- c(t_end - t_beginning, min_age)
  
  return(list(nodes = nodes, 
              times = times,
              t_end = t_end,
              t_beginning = t_beginning,
              regime_time = regime_time,
              which.regimes = which.regimes))
}

weights_segments <- function(a, lineage){
  res <- c(exp(-a * lineage$t_beginning) - exp(-a * lineage$t_end), 
           exp(-a * lineage$times[1]))
  return(res)
}

weights_regimes <- function(a, lineage) {
  #nt <- lineage$nodes_time
  res <- weights_segments(a, lineage) ## Rcpp wrapper, equivalent to above commented out code
  w <- vapply(lineage$which.regimes, function(e) sum(e*res), FUN.VALUE = 0) ## Sum coefficients for which the respective regime is equal
  return(w)
}

weights_regimes_brown <- function(lineage){
  w <- sapply(lineage$which.regimes, function(e) sum(e*lineage$regime_time))
  return(w)
}

weight.matrix <- function(phy, a, lineages){
  if(a > 300000000000) a <- 300000000000
  res <- t(vapply(lineages, function(x) weights_regimes(a, x), 
                  FUN.VALUE = numeric(length(lineages[[1]]$which.regimes))) ## Specify type of output
           )

  rownames(res) <- phy$tip.label
  return(res)
}

weight.matrix.brown <- function(lineages){
  res <- t(sapply(lineages, weights_regimes_brown))
  return(res)
}

## Thanks to user "snaut" at stackoverflow, http://stackoverflow.com/users/1999873/snaut

concat.factor <- function(...){
  as.factor(do.call(c, lapply(list(...), as.character)))
}