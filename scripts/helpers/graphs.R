makeGraph <- function(rho, theta, beta = 2) {
  adj <- abs(rho)^beta
  adj[upper.tri(adj, diag = TRUE)] <- 0
  adj[adj < theta] <- 0
  gr <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE, weighted = TRUE)
}

makeUnionGraph <- function(x, y) {
  ug <- igraph::union(x, y)
  igraph::E(ug)$weight <- apply(cbind(igraph::E(ug)$weight_1, igraph::E(ug)$weight_2), 1, max, na.rm = TRUE)
  ug <- igraph::simplify(ug)
  ug
}

findCommunities <- function(graph) {
  out <- cluster_infomap(graph)
  out <- communities(out)
  out <- out[sapply(out, length) >= 10]
  out <- out[order(- sapply(out, length))]
  names(out) <- paste0("Community", seq_along(out))
  out
}

printCommunities <- function(x, geneid) {
  for (i in names(x)) {
    cat(i)
    cat("\n")
    cat(sort(geneid[x[[i]]]$GeneSymbol))
    cat("\n\n")
  }
}
