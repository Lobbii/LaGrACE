#' Load reference sample causal network
#'
#' this functions loads reference sample causal network (FGES) from TXT file
#' @param file Path to TXT file containing the reference sample network
#' @return A list containing edge, nodes and adjacency information of causal graph learned by FGES
#' @export

read_causality_graph<-function(file){
  if (!file.exists(file)) {
    stop("Cannot find file!")
  }
  tmp <- readLines(file)
  if (tmp[1] != "Graph Nodes:") {    stop("file does not contain a compatible graph.")}

  if (grepl(";", tmp[2])){
    split <- ";"
  }  else {split <- ","}

  nodes <- unlist(strsplit(tmp[2], split = split))
  nodes <- sort(nodes)
  if (tmp[4] != "Graph Edges:") {stop("file does not contain a compatible graph.")}

  tmp <- tmp[-(1:4)]
  edges <- c()
  line <- tmp[1]
  while (!is.na(line) && line != "") {
    line <- sub("[0-9]+\\. ", "", line)
    edge <- unlist(strsplit(line, split = " "))
    if (is.na(sum(match(c(edge[1], edge[3]), nodes)))) {
      print(line)
      stop("file contains a malformed causality graph.")
    }
    edges <- c(edges, c(edge[1], edge[3], edge[2]))
    tmp <- tmp[-1]
    line <- tmp[1]
  }
  edges<-matrix(edges, ncol=3, byrow=TRUE)

  adjacencies <- lapply(nodes, function(node) {
    neighborhood <- c()
    # loop over the edges
    for (i in 1:nrow(edges)) {
      # if a node is in an edge, find its partner (neighbor)
      if( node %in% edges[i, 1:2]) {
        neighbor <- edges[i, c(node != edges[i, 1:2], F)]
        if (!(neighbor %in% neighborhood))
          neighborhood <- c(neighborhood, neighbor)
      }
    }
    return(neighborhood)
  })
  names(adjacencies) <- nodes

  return(list('edges'=edges,'nodes'=nodes,'adjacencies'=adjacencies))
}



