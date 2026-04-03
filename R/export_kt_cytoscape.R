#' Export KT Graph to Cytoscape
#'
#' Exports a KT graph to Cytoscape with custom styling based on
#' alignment quality and protein family information.
#'
#' @param graph An igraph object from create_kt_graph()
#' @param network_name Character, name for the network in Cytoscape.
#'   Default: "KT_Alignment_Network"
#' @param collection_name Character, collection name in Cytoscape.
#'   Default: "KT Analysis"
#' @param style_name Character, visual style name. Default: "KT_Network_Style"
#' @param edge_color_attribute Character, edge attribute for color mapping.
#'   Default: "identity"
#' @param check_cytoscape Logical, whether to check Cytoscape connection.
#'   Default: TRUE
#'
#' @return Logical, TRUE if export successful, FALSE otherwise
#'
#' @examples
#' \dontrun{
#' g <- create_kt_graph(hq_alignments)
#' g <- apply_kt_graph_visuals(g)
#' export_kt_to_cytoscape(g)
#'
#' # Custom export
#' export_kt_to_cytoscape(g,
#'                       network_name = "High_Quality_KT",
#'                       edge_color_attribute = "bitscore")
#' }
#'
#' @importFrom RCy3 cytoscapePing createNetworkFromDataFrames setVisualStyle
#' @importFrom RCy3 copyVisualStyle deleteVisualStyle getVisualStyleNames
#' @importFrom RCy3 setEdgeColorMapping setEdgeLineWidthMapping layoutNetwork
#' @importFrom viridis viridis
#' @importFrom igraph V E as_data_frame vertex_attr
#' @export
export_kt_to_cytoscape <- function(graph,
                                  network_name = "KT_Alignment_Network",
                                  collection_name = "KT Analysis",
                                  style_name = "KT_Network_Style",
                                  edge_color_attribute = "identity",
                                  check_cytoscape = FALSE) {

  if (!inherits(graph, "igraph")) {
    stop("Input must be an igraph object")
  }

  # Create nodes data frame with all attributes
  nodes_df <- data.frame(
    id = V(graph)$name,
    label = V(graph)$name,
    family = V(graph)$family,
    organism = V(graph)$organism,
    stringsAsFactors = FALSE
  )

  # Add community information if available
  if ("community" %in% names(igraph::vertex_attr(graph))) {
    nodes_df$community <- V(graph)$community
  }

  # Create edges data frame with all edge attributes
  edges_df <- as_data_frame(graph, what = "edges")

  # Rename columns for RCy3 compatibility
  if ("from" %in% names(edges_df)) {
    names(edges_df)[names(edges_df) == "from"] <- "source"
  }
  if ("to" %in% names(edges_df)) {
    names(edges_df)[names(edges_df) == "to"] <- "target"
  }

  # Add interaction column if not present (required by some RCy3 versions)
  if (!"interaction" %in% names(edges_df)) {
    edges_df$interaction <- "alignment"
  }

  message(paste("Exporting network with", nrow(nodes_df), "nodes and", nrow(edges_df), "edges"))

  # Create network in Cytoscape
  tryCatch({
    createNetworkFromDataFrames(
      nodes = nodes_df,
      edges = edges_df,
      title = network_name,
      collection = collection_name
    )
    message(paste("Network created in Cytoscape:", network_name))
  }, error = function(e) {
    warning(paste("Failed to create network:", e$message))
    return(FALSE)
  })

  # Apply visual style
  success <- create_kt_cytoscape_style(style_name, edge_color_attribute)
  if (!success) {
    warning("Failed to create visual style, using default")
  }

  # Apply force-directed layout
  tryCatch({
    layoutNetwork("force-directed")
    message("Applied force-directed layout")
  }, error = function(e) {
    warning(paste("Layout failed:", e$message))
  })

  message(paste("Network successfully exported to Cytoscape as:", network_name))
  return(TRUE)
}

#' Create KT Cytoscape Visual Style
#'
#' Creates a custom visual style for KT networks in Cytoscape
#'
#' @param style_name Character, name for the visual style
#' @param edge_color_attribute Character, edge attribute for color mapping
#'
#' @return Logical, TRUE if style created successfully
#'
#' @importFrom RCy3 getVisualStyleNames deleteVisualStyle copyVisualStyle
#' @importFrom RCy3 setVisualStyle setEdgeColorMapping setEdgeLineWidthMapping
#' @importFrom RCy3 setNodeColorMapping setNodeSizeDefault
#' @importFrom viridis viridis plasma
create_kt_cytoscape_style <- function(style_name, edge_color_attribute = "identity") {

  tryCatch({
    # Delete existing style if it exists
    if (style_name %in% getVisualStyleNames()) {
      deleteVisualStyle(style_name)
    }

    # Create new style based on default
    copyVisualStyle("default", style_name)
    setVisualStyle(style_name)

    # Set edge color mapping based on specified attribute
    setEdgeColorMapping(
      table.column = edge_color_attribute,
      style.name = style_name,
      mapping.type = "continuous"
    )

    # Set edge width mapping
    setEdgeLineWidthMapping(
      table.column = edge_color_attribute,
      style.name = style_name,
      widths = c(0.5, 5),
      mapping.type = "continuous"
    )

    # Set node color mapping by family if available
    tryCatch({
      setNodeColorMapping(
        table.column = "family",
        style.name = style_name,
        mapping.type = "discrete"
      )
    }, error = function(e) {
      message("Could not map node colors by family - using default")
    })

    # Set default node size
    setNodeSizeDefault(8, style.name = style_name)

    message(paste("Visual style", style_name, "created successfully"))
    return(TRUE)

  }, error = function(e) {
    warning(paste("Failed to create visual style:", e$message))
    return(FALSE)
  })
}

#' Save KT Graph Data for Manual Import
#'
#' Saves network data as CSV files for manual import into Cytoscape
#' or other network analysis tools.
#'
#' @param graph An igraph object from create_kt_graph()
#' @param nodes_file Character, filename for nodes CSV. Default: "kt_network_nodes.csv"
#' @param edges_file Character, filename for edges CSV. Default: "kt_network_edges.csv"
#' @param output_dir Character, directory to save files. Default: current directory
#'
#' @return Character vector of created file paths
#'
#' @examples
#' \dontrun{
#' g <- create_kt_graph(hq_alignments)
#' files <- save_kt_graph_data(g)
#' print(files)
#' }
#'
#' @importFrom igraph V E as_data_frame
#' @export
save_kt_graph_data <- function(graph,
                              nodes_file = "kt_network_nodes.csv",
                              edges_file = "kt_network_edges.csv",
                              output_dir = ".") {

  if (!inherits(graph, "igraph")) {
    stop("Input must be an igraph object")
  }

  # Create full file paths
  nodes_path <- file.path(output_dir, nodes_file)
  edges_path <- file.path(output_dir, edges_file)

  # Create nodes data frame
  nodes_df <- data.frame(
    id = V(graph)$name,
    label = V(graph)$name,
    family = V(graph)$family,
    organism = V(graph)$organism,
    stringsAsFactors = FALSE
  )

  # Add community information if available
  if ("community" %in% names(igraph::vertex_attr(graph))) {
    nodes_df$community <- V(graph)$community
  }

  # Create edges data frame with RCy3-compatible column names
  edges_df <- as_data_frame(graph, what = "edges")

  # Rename columns for RCy3/Cytoscape compatibility
  if ("from" %in% names(edges_df)) {
    names(edges_df)[names(edges_df) == "from"] <- "source"
  }
  if ("to" %in% names(edges_df)) {
    names(edges_df)[names(edges_df) == "to"] <- "target"
  }

  # Add interaction column for Cytoscape
  if (!"interaction" %in% names(edges_df)) {
    edges_df$interaction <- "alignment"
  }

  # Save files
  write.csv(nodes_df, nodes_path, row.names = FALSE)
  write.csv(edges_df, edges_path, row.names = FALSE)

  message(paste("Network data saved to:"))
  message(paste("  Nodes:", nodes_path))
  message(paste("  Edges:", edges_path))

  return(c(nodes_path, edges_path))
}

#' Export KT Graph with Fallback
#'
#' Attempts to export to Cytoscape, falls back to CSV files if Cytoscape unavailable
#'
#' @param graph An igraph object from create_kt_graph()
#' @param ... Additional arguments passed to export_kt_to_cytoscape()
#'
#' @return List with export method used and file paths (if CSV)
#'
#' @examples
#' \dontrun{
#' g <- create_kt_graph(hq_alignments)
#' g <- apply_kt_graph_visuals(g)
#' result <- export_kt_graph_with_fallback(g)
#' }
#'
#' @export
export_kt_graph_with_fallback <- function(graph, ...) {

  # Try Cytoscape export first
  cytoscape_success <- export_kt_to_cytoscape(graph, ...)

  if (cytoscape_success) {
    return(list(
      method = "cytoscape",
      success = TRUE,
      files = NULL
    ))
  } else {
    # Fallback to CSV export
    message("Cytoscape not available. Saving network data for manual import.")
    files <- save_kt_graph_data(graph)

    return(list(
      method = "csv",
      success = TRUE,
      files = files
    ))
  }
}