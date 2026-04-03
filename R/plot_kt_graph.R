#' Apply Visual Properties to KT Graph
#'
#' Applies visual styling to a KT graph including edge colors, widths,
#' and node properties based on network attributes.
#'
#' @param graph An igraph object from create_kt_graph()
#' @param edge_color_by Character, edge attribute to use for coloring.
#'   Default: "identity"
#' @param color_palette Function or vector, color palette for edges.
#'   Default: viridis::viridis
#' @param edge_width_range Numeric vector of length 2, range for edge widths.
#'   Default: c(0.5, 3)
#' @param node_size Numeric, size for all nodes. Default: 8
#' @param node_color_by Character, node attribute for coloring or single color.
#'   Default: "lightblue"
#'
#' @return Modified igraph object with visual properties applied
#'
#' @examples
#' \dontrun{
#' g <- create_kt_graph(hq_alignments)
#' g <- apply_kt_graph_visuals(g)
#'
#' # Custom styling
#' g <- apply_kt_graph_visuals(g,
#'                            edge_color_by = "bitscore",
#'                            color_palette = rainbow,
#'                            node_color_by = "family")
#' }
#'
#' @importFrom viridis viridis
#' @importFrom scales rescale
#' @importFrom igraph V E vertex_attr edge_attr
#' @export
apply_kt_graph_visuals <- function(graph,
                                  edge_color_by = "identity",
                                  color_palette = viridis,
                                  edge_width_range = c(0.5, 3),
                                  node_size = 8,
                                  node_color_by = "lightblue") {

  if (!inherits(graph, "igraph")) {
    stop("Input must be an igraph object")
  }

  # Apply edge colors based on specified attribute
  if (edge_color_by %in% names(igraph::edge_attr(graph))) {
    edge_values <- igraph::edge_attr(graph, edge_color_by)
    value_range <- range(edge_values, na.rm = TRUE)

    message(paste("Edge color range (", edge_color_by, "):",
                  round(value_range[1], 2), "to", round(value_range[2], 2)))

    # Generate color palette
    if (is.function(color_palette)) {
      colors <- color_palette(100)
    } else {
      colors <- color_palette
    }

    # Map values to colors
    color_indices <- as.numeric(cut(edge_values, breaks = 100))
    E(graph)$color <- colors[color_indices]

    # Set edge width based on the same attribute
    E(graph)$width <- rescale(edge_values, to = edge_width_range)
  } else {
    warning(paste("Edge attribute", edge_color_by, "not found. Using default colors."))
    E(graph)$color <- "gray"
    E(graph)$width <- 1
  }

  # Apply node properties
  V(graph)$size <- node_size

  # Node colors
  if (length(node_color_by) == 1 && node_color_by %in% names(igraph::vertex_attr(graph))) {
    # Color by attribute (e.g., family)
    node_values <- igraph::vertex_attr(graph, node_color_by)
    unique_values <- unique(node_values)
    n_colors <- length(unique_values)

    if (is.function(color_palette)) {
      node_colors <- color_palette(n_colors)
    } else {
      node_colors <- rainbow(n_colors) # fallback
    }

    color_map <- setNames(node_colors, unique_values)
    V(graph)$color <- color_map[node_values]

    message(paste("Nodes colored by", node_color_by, "with", n_colors, "categories"))
  } else {
    # Single color for all nodes
    V(graph)$color <- node_color_by
  }

  message("Visual properties applied to graph")
  return(graph)
}

#' Plot KT Graph
#'
#' Create a plot of the KT graph using base R graphics
#'
#' @param graph An igraph object with visual properties applied
#' @param main Character, title for the plot. Default: "KT Protein Network"
#' @param vertex_label_size Numeric, size of vertex labels. Default: 0.8
#' @param show_labels Logical, whether to show vertex labels. Default: FALSE
#' @param legend_position Character, position for legend if coloring by attribute.
#'   Options: "topleft", "topright", "bottomleft", "bottomright", "none"
#'   Default: "topright"
#'
#' @return Plots the graph (invisible return)
#'
#' @examples
#' \dontrun{
#' g <- create_kt_graph(hq_alignments)
#' g <- apply_kt_graph_visuals(g, node_color_by = "family")
#' plot_kt_graph(g, show_labels = TRUE)
#' }
#'
#' @importFrom igraph plot.igraph V
#' @export
plot_kt_graph <- function(graph,
                         main = "KT Protein Network",
                         vertex_label_size = 0.8,
                         show_labels = FALSE,
                         legend_position = "topright") {

  if (!inherits(graph, "igraph")) {
    stop("Input must be an igraph object")
  }

  # Set up plot parameters
  layout_coords <- cbind(V(graph)$x, V(graph)$y)

  # Plot parameters
  plot_params <- list(
    x = graph,
    layout = layout_coords,
    main = main,
    vertex.size = V(graph)$size,
    vertex.color = V(graph)$color,
    vertex.label.cex = vertex_label_size,
    edge.color = E(graph)$color,
    edge.width = E(graph)$width
  )

  # Handle labels
  if (show_labels) {
    plot_params$vertex.label <- V(graph)$label
  } else {
    plot_params$vertex.label <- NA
  }

  # Create plot
  do.call(plot, plot_params)

  # Add legend if nodes are colored by attribute
  if ("family" %in% names(igraph::vertex_attr(graph)) && legend_position != "none") {
    families <- unique(V(graph)$family)
    family_colors <- unique(V(graph)$color[!is.na(V(graph)$family)])

    if (length(families) <= 20 && length(families) == length(family_colors)) {
      legend(legend_position,
             legend = families,
             fill = family_colors,
             title = "KT Family",
             cex = 0.7)
    }
  }

  invisible(NULL)
}

#' Create KT Graph Summary Plot
#'
#' Creates a multi-panel summary plot showing network statistics
#'
#' @param graph An igraph object from create_kt_graph()
#' @param stats_list Output from get_kt_graph_stats()
#'
#' @return Creates a multi-panel plot (invisible return)
#'
#' @examples
#' \dontrun{
#' g <- create_kt_graph(hq_alignments)
#' stats <- get_kt_graph_stats(g)
#' plot_kt_graph_summary(g, stats)
#' }
#'
#' @importFrom graphics par hist barplot
#' @importFrom igraph degree
#' @export
plot_kt_graph_summary <- function(graph, stats_list) {

  # Set up multi-panel plot
  old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  on.exit(par(old_par))

  # 1. Degree distribution
  degree_dist <- degree(graph)
  hist(degree_dist,
       main = "Degree Distribution",
       xlab = "Node Degree",
       ylab = "Frequency",
       col = "lightblue",
       border = "darkblue")

  # 2. Family distribution
  family_counts <- stats_list$family_stats$family_counts
  top_families <- head(family_counts, 15)
  barplot(top_families,
          main = "Top Families by Node Count",
          xlab = "Family",
          ylab = "Number of Nodes",
          col = "lightgreen",
          las = 2,
          cex.names = 0.7)

  # 3. Edge identity distribution
  edge_identity <- E(graph)$identity
  hist(edge_identity,
       main = "Alignment Identity Distribution",
       xlab = "Percent Identity",
       ylab = "Frequency",
       col = "lightcoral",
       border = "darkred")

  # 4. Network overview with community colors
  if ("community" %in% names(igraph::vertex_attr(graph))) {
    communities <- V(graph)$community
    n_communities <- max(communities)
    community_colors <- rainbow(n_communities)[communities]
    V(graph)$color <- community_colors
  }

  plot_kt_graph(graph,
               main = paste("Network Overview\n(",
                           stats_list$basic_stats$nodes, "nodes,",
                           stats_list$basic_stats$edges, "edges)"),
               show_labels = FALSE,
               legend_position = "none")

  invisible(NULL)
}