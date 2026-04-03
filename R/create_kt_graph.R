#' Create Network Graph from KT Alignments
#'
#' Converts KT alignment data into an undirected igraph object where nodes represent
#' protein accessions and edges represent alignments. Automatically handles reciprocal
#' alignments by selecting the best alignment per protein pair (highest bitscore).
#'
#' @param alignment_data A data.table from load_kt_alignment() or filter_alignments()
#' @param layout_method Character, layout method for positioning nodes.
#'   Options: "fruchterman.reingold", "kamada.kawai", "circle", "random"
#'   Default: "fruchterman.reingold"
#' @param layout_iterations Integer, number of iterations for force-directed layouts
#'   Default: 1000
#'
#' @return An undirected igraph object with:
#' \describe{
#'   \item{Nodes}{Unique accessions with family, organism attributes}
#'   \item{Edges}{Deduplicated alignments with identity, length, evalue, bitscore}
#'   \item{Layout}{X,Y coordinates for node positioning}
#' }
#'
#' @details
#' The function automatically detects reciprocal alignments (A->B and B->A) and
#' deduplicates them for the undirected graph by selecting the alignment with the
#' highest bitscore. This ensures each protein pair has only one edge while
#' preserving the best alignment statistics.
#'
#' @examples
#' \dontrun{
#' # Create graph from filtered alignments
#' hq_alignments <- filter_alignments(alignments, min_identity = 30)
#' g <- create_kt_graph(hq_alignments)
#'
#' # With custom layout
#' g <- create_kt_graph(hq_alignments, layout_method = "kamada.kawai")
#' }
#'
#' @importFrom igraph graph_from_data_frame V E layout_with_fr layout_with_kk
#' @importFrom igraph layout_in_circle layout_randomly
#' @importFrom data.table data.table
#' @export
create_kt_graph <- function(alignment_data,
                           layout_method = "fruchterman.reingold",
                           layout_iterations = 10000) {

  if (!inherits(alignment_data, "data.table")) {
    stop("Input must be a data.table from load_kt_alignment() or filter_alignments()")
  }

  if (nrow(alignment_data) == 0) {
    stop("No alignments found in input data")
  }

  # Create node attributes data frame
  query_nodes <- alignment_data[, .(
    accession = query_accession,
    family = query_kt_family,
    organism = query_organism
  )]

  subject_nodes <- alignment_data[, .(
    accession = subject_accession,
    family = subject_kt_family,
    organism = subject_organism
  )]

  # Combine and get unique nodes with their attributes
  all_nodes <- unique(rbind(query_nodes, subject_nodes))

  message(paste("Creating network with", nrow(all_nodes), "nodes"))

  # Create edge list and handle reciprocal alignments
  initial_edge_list <- alignment_data[, .(
    from = query_accession,
    to = subject_accession,
    identity = percent_identity,
    length = alignment_length,
    evalue = evalue,
    bitscore = bitscore
  )]

  # Check for reciprocal alignments and deduplicate for undirected graph
  # Create a canonical pair identifier (alphabetically sorted)
  initial_edge_list[, pair_id := paste(
    pmin(from, to),
    pmax(from, to),
    sep = "_"
  )]

  # Count reciprocal pairs
  total_edges <- nrow(initial_edge_list)
  unique_pairs <- uniqueN(initial_edge_list$pair_id)
  reciprocal_ratio <- total_edges / unique_pairs

  message(paste("Initial edges:", total_edges))
  message(paste("Unique pairs:", unique_pairs))
  message(paste("Reciprocal ratio:", round(reciprocal_ratio, 2)))

  # For undirected graph, keep only one edge per pair (the best one)
  if (reciprocal_ratio > 1.0) {
    message("Reciprocal alignments detected - selecting best alignment per pair")

    # Select the best alignment per pair (highest bitscore, then highest identity)
    edge_list <- initial_edge_list[
      order(-bitscore, -identity),
      .SD[1],
      by = pair_id
    ][, pair_id := NULL]

    message(paste("Deduplicated to", nrow(edge_list), "edges for undirected graph"))
  } else {
    message("No significant reciprocal alignments found - using original edge list")
    edge_list <- initial_edge_list[, pair_id := NULL]
  }

  # Create igraph object as undirected graph
  g <- graph_from_data_frame(d = edge_list, vertices = all_nodes, directed = FALSE)

  # Node attributes are automatically assigned from vertices data frame
  # Additional node properties
  V(g)$id <- V(g)$name
  V(g)$label <- V(g)$name

  # Apply layout
  layout_coords <- switch(layout_method,
    "fruchterman.reingold" = layout_with_fr(g, niter = layout_iterations),
    "kamada.kawai" = layout_with_kk(g),
    "circle" = layout_in_circle(g),
    "random" = layout_randomly(g),
    layout_with_fr(g, niter = layout_iterations) # default fallback
  )

  # Add coordinates to graph
  V(g)$x <- layout_coords[, 1]
  V(g)$y <- layout_coords[, 2]

  # Add layout method as graph attribute
  g$layout_method <- layout_method

  message("Network graph created successfully")
  message(paste("Families represented:", length(unique(V(g)$family))))
  message(paste("Organisms represented:", length(unique(V(g)$organism))))

  return(g)
}

#' Get Network Statistics
#'
#' Calculate comprehensive network statistics for a KT graph
#'
#' @param graph An igraph object from create_kt_graph()
#'
#' @return A list with network statistics and family analysis
#'
#' @examples
#' \dontrun{
#' g <- create_kt_graph(hq_alignments)
#' stats <- get_kt_graph_stats(g)
#' print(stats$basic_stats)
#' }
#'
#' @importFrom igraph vcount ecount edge_density transitivity average.path.length
#' @importFrom igraph degree cluster_louvain modularity membership
#' @importFrom igraph as_data_frame
#' @export
get_kt_graph_stats <- function(graph) {

  # Basic network statistics
  basic_stats <- list(
    nodes = vcount(graph),
    edges = ecount(graph),
    density = round(edge_density(graph), 4),
    clustering = round(transitivity(graph, type = "average"), 4),
    avg_path_length = round(average.path.length(graph), 2)
  )

  # Degree distribution
  degree_dist <- degree(graph)
  degree_stats <- list(
    avg_degree = round(mean(degree_dist), 2),
    max_degree = max(degree_dist),
    top_nodes = names(sort(degree_dist, decreasing = TRUE)[1:10])
  )

  # Community detection
  communities <- cluster_louvain(graph)
  community_stats <- list(
    num_communities = length(communities),
    modularity = round(modularity(communities), 3),
    membership = membership(communities)
  )

  # Add community info to graph
  V(graph)$community <- community_stats$membership

  # Family-based analysis
  family_counts <- table(V(graph)$family)

  # Analyze connectivity within vs between families
  edge_df <- as_data_frame(graph, what = "edges")
  node_df <- as_data_frame(graph, what = "vertices")

  # Add family info to edges by matching node names
  edge_df$from_family <- node_df$family[match(edge_df$from, node_df$name)]
  edge_df$to_family <- node_df$family[match(edge_df$to, node_df$name)]

  intra_family <- sum(edge_df$from_family == edge_df$to_family, na.rm = TRUE)
  inter_family <- sum(edge_df$from_family != edge_df$to_family, na.rm = TRUE)

  # Top connected families
  family_node_degrees <- aggregate(degree_dist, by = list(V(graph)$family), FUN = mean)
  colnames(family_node_degrees) <- c("family", "avg_degree")
  family_node_degrees <- family_node_degrees[order(family_node_degrees$avg_degree, decreasing = TRUE), ]

  family_stats <- list(
    family_counts = sort(family_counts, decreasing = TRUE),
    intra_family_edges = intra_family,
    inter_family_edges = inter_family,
    inter_intra_ratio = round(inter_family/intra_family, 2),
    top_families_by_degree = head(family_node_degrees, 10)
  )

  return(list(
    basic_stats = basic_stats,
    degree_stats = degree_stats,
    community_stats = community_stats,
    family_stats = family_stats,
    graph_updated = graph
  ))
}