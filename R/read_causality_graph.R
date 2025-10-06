#' Load reference sample causal network
#'
#' Loads the FGES causal network from a TXT file with the expected format and
#' returns edges, nodes, and adjacencies.
#'
#' @param file Path to TXT file containing the reference sample network
#' @return A list with `edges`, `nodes`, and `adjacencies`
#' @export

read_causality_graph <- function(file) {
  if (!file.exists(file)) {
    stop(sprintf("Cannot find file: %s", file))
  }
  tmp <- readLines(file, warn = FALSE)
  if (length(tmp) < 4L) {
    stop("Network file is too short to be valid.")
  }
  if (tmp[1] != "Graph Nodes:") {
    stop("File does not contain a compatible graph: missing 'Graph Nodes:' header.")
  }

  split <- if (grepl(";", tmp[2], fixed = TRUE)) ";" else ","
  nodes <- unlist(strsplit(tmp[2], split = split, fixed = TRUE))
  nodes <- sort(nodes)
  if (tmp[4] != "Graph Edges:") {
    stop("File does not contain a compatible graph: missing 'Graph Edges:' header.")
  }

  tmp <- tmp[-(1:4)]
  edges <- character(0)
  line <- tmp[1]
  while (!is.na(line) && nzchar(line)) {
    line <- sub("[0-9]+\\. ", "", line)
    edge <- unlist(strsplit(line, split = " "))
    if (length(edge) < 3L || any(is.na(match(c(edge[1], edge[3]), nodes)))) {
      stop(sprintf("File contains malformed graph line: '%s'", line))
    }
    edges <- c(edges, c(edge[1], edge[3], edge[2]))
    tmp <- tmp[-1]
    line <- tmp[1]
  }
  edges <- matrix(edges, ncol = 3, byrow = TRUE)

  adjacencies <- lapply(nodes, function(node) {
    neighborhood <- character(0)
    for (i in seq_len(nrow(edges))) {
      # if a node is in an edge, find its partner (neighbor)
      if (node %in% edges[i, 1:2]) {
        neighbor <- edges[i, c(node != edges[i, 1:2], FALSE)]
        if (!(neighbor %in% neighborhood)) {
          neighborhood <- c(neighborhood, neighbor)
        }
      }
    }
    neighborhood
  })
  names(adjacencies) <- nodes

  return(list(edges = edges, nodes = nodes, adjacencies = adjacencies))
}



