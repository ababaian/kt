#' Plot Cluster Size Distribution
#'
#' Creates a histogram showing the distribution of cluster sizes (number of
#' members per cluster centroid) from UCLUST clustering results.
#'
#' @param clustering_data A data.table from load_uclust()
#' @param title Character, plot title. Default: "Cluster Size Distribution"
#' @param bins Integer, number of histogram bins. Default: 30
#' @param log_scale Logical, whether to use log10 scale for x-axis. Default: FALSE
#' @param color Character, fill color for histogram. Default: "steelblue"
#' @param alpha Numeric, transparency for histogram bars. Default: 0.7
#' @param show_stats Logical, whether to show summary statistics. Default: TRUE
#' @param facet_by Character, variable to facet by ("kt_family", "none"). Default: "none"
#' @param min_cluster_size Integer, minimum cluster size to include. Default: 0
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' clustering_data <- load_uclust("clusters.uc", include_hits = TRUE)
#' plot_cluster_counts(clustering_data)
#'
#' # Log scale for large range
#' plot_cluster_counts(clustering_data, log_scale = TRUE)
#'
#' # By family
#' plot_cluster_counts(clustering_data, facet_by = "kt_family")
#'
#' # Filter small clusters
#' plot_cluster_counts(clustering_data, min_cluster_size = 2)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_histogram labs theme_minimal theme
#' @importFrom ggplot2 element_text geom_vline annotate scale_x_log10 facet_wrap
#' @importFrom dplyr count filter summarise
#' @importFrom data.table data.table
#' @export
plot_cluster_counts <- function(clustering_data,
                               title = "Cluster Size Distribution",
                               bins = 30,
                               log_scale = FALSE,
                               color = "steelblue",
                               alpha = 0.7,
                               show_stats = TRUE,
                               facet_by = "none",
                               min_cluster_size = 0) {

  if (!inherits(clustering_data, "data.table")) {
    stop("Input must be a data.table from load_uclust()")
  }

  # Check required columns
  required_cols <- c("record_type", "cluster_id")
  missing_cols <- setdiff(required_cols, names(clustering_data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  # Count members per cluster (H records = cluster members)
  cluster_sizes <- clustering_data[record_type == "H", .N, by = .(cluster_id, kt_family)]

  # Add clusters with zero members (singletons)
  all_clusters <- clustering_data[record_type == "C", .(cluster_id, kt_family)]
  cluster_sizes <- merge(all_clusters, cluster_sizes, by = c("cluster_id", "kt_family"), all.x = TRUE)
  cluster_sizes[is.na(N), N := 0]  # Replace NA with 0 for singleton clusters

  # Filter by minimum cluster size
  if (min_cluster_size > 0) {
    cluster_sizes <- cluster_sizes[N >= min_cluster_size]
    message(paste("Filtered to", nrow(cluster_sizes), "clusters with >=", min_cluster_size, "members"))
  }

  # Summary statistics
  total_clusters <- nrow(cluster_sizes)
  total_members <- sum(cluster_sizes$N)
  mean_size <- round(mean(cluster_sizes$N), 2)
  median_size <- median(cluster_sizes$N)
  max_size <- max(cluster_sizes$N)
  singletons <- sum(cluster_sizes$N == 0)

  message(paste("Cluster Summary:"))
  message(paste("  Total clusters:", total_clusters))
  message(paste("  Total members:", total_members))
  message(paste("  Singleton clusters:", singletons, paste0("(", round(100*singletons/total_clusters, 1), "%)")))
  message(paste("  Mean cluster size:", mean_size))
  message(paste("  Median cluster size:", median_size))
  message(paste("  Largest cluster:", max_size, "members"))

  # Create base plot
  p <- ggplot(cluster_sizes, aes(x = N))

  # Add histogram
  p <- p + geom_histogram(bins = bins, fill = color, alpha = alpha, color = "white")

  # Apply log scale if requested
  if (log_scale) {
    p <- p + scale_x_log10(labels = scales::comma_format()) +
      labs(x = "Cluster Size (log10 scale)")
  } else {
    p <- p + labs(x = "Cluster Size (number of members)")
  }

  # Add faceting if requested
  if (facet_by == "kt_family") {
    p <- p + facet_wrap(~kt_family, scales = "free_y")
    plot_title <- paste(title, "by Family")
  } else {
    plot_title <- title
  }

  # Add statistics if requested
  if (show_stats && facet_by == "none") {
    p <- p +
      geom_vline(xintercept = mean_size, color = "red", linetype = "dashed", size = 1) +
      geom_vline(xintercept = median_size, color = "blue", linetype = "dashed", size = 1) +
      annotate("text",
               x = Inf, y = Inf,
               label = paste0("Total clusters: ", total_clusters, "\n",
                             "Singletons: ", singletons, " (", round(100*singletons/total_clusters, 1), "%)\n",
                             "Mean: ", mean_size, "\n",
                             "Median: ", median_size, "\n",
                             "Max: ", max_size),
               hjust = 1.1, vjust = 1.1,
               color = "black", size = 3.5,
               fontface = "bold")
  }

  # Apply theme and final labels
  p <- p +
    labs(
      title = plot_title,
      y = "Number of Clusters"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      strip.text = element_text(size = 10, face = "bold")
    )

  return(p)
}

#' Get Detailed Cluster Statistics
#'
#' Calculates comprehensive statistics about cluster size distribution
#' from UCLUST clustering results.
#'
#' @param clustering_data A data.table from load_uclust()
#'
#' @return A list with detailed clustering statistics
#'
#' @examples
#' \dontrun{
#' clustering_data <- load_uclust("clusters.uc", include_hits = TRUE)
#' stats <- get_cluster_statistics(clustering_data)
#' print(stats$summary)
#' }
#'
#' @importFrom dplyr summarise group_by arrange desc
#' @importFrom data.table data.table
#' @export
get_cluster_statistics <- function(clustering_data) {

  if (!inherits(clustering_data, "data.table")) {
    stop("Input must be a data.table from load_uclust()")
  }

  # Count members per cluster
  cluster_sizes <- clustering_data[record_type == "H", .N, by = .(cluster_id, kt_family)]

  # Add singletons (clusters with no members)
  all_clusters <- clustering_data[record_type == "C", .(cluster_id, kt_family)]
  cluster_sizes <- merge(all_clusters, cluster_sizes, by = c("cluster_id", "kt_family"), all.x = TRUE)
  cluster_sizes[is.na(N), N := 0]

  # Overall statistics
  summary_stats <- list(
    total_clusters = nrow(cluster_sizes),
    total_members = sum(cluster_sizes$N),
    singletons = sum(cluster_sizes$N == 0),
    mean_size = round(mean(cluster_sizes$N), 2),
    median_size = as.integer(median(cluster_sizes$N)),
    min_size = min(cluster_sizes$N),
    max_size = max(cluster_sizes$N),
    q25_size = as.integer(quantile(cluster_sizes$N, 0.25)),
    q75_size = as.integer(quantile(cluster_sizes$N, 0.75))
  )

  # Size distribution
  size_distribution <- as.data.frame(table(cluster_sizes$N))
  colnames(size_distribution) <- c("cluster_size", "frequency")
  size_distribution$cluster_size <- as.numeric(as.character(size_distribution$cluster_size))
  size_distribution <- size_distribution[order(size_distribution$cluster_size), ]

  # Family-wise statistics
  family_stats <- cluster_sizes[, .(
    clusters = .N,
    total_members = sum(N),
    singletons = sum(N == 0),
    mean_size = round(mean(N), 2),
    median_size = as.integer(median(N)),
    max_size = max(N)
  ), by = kt_family][order(-clusters)]

  # Top largest clusters
  top_clusters <- cluster_sizes[order(-N)][1:min(10, nrow(cluster_sizes))]

  return(list(
    summary = summary_stats,
    size_distribution = size_distribution,
    family_stats = family_stats,
    top_clusters = top_clusters,
    cluster_sizes = cluster_sizes
  ))
}

#' Plot Cluster Size vs Identity
#'
#' Creates a scatter plot showing the relationship between cluster size
#' and clustering identity threshold effects.
#'
#' @param clustering_data A data.table from load_uclust()
#' @param min_cluster_size Integer, minimum cluster size to include. Default: 1
#' @param alpha Numeric, point transparency. Default: 0.6
#' @param color_by Character, variable to color points by. Default: "kt_family"
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' clustering_data <- load_uclust("clusters.uc", include_hits = TRUE)
#' plot_cluster_size_identity(clustering_data)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_minimal theme
#' @importFrom ggplot2 element_text scale_color_viridis_d scale_y_log10
#' @importFrom dplyr left_join
#' @export
plot_cluster_size_identity <- function(clustering_data,
                                      min_cluster_size = 1,
                                      alpha = 0.6,
                                      color_by = "kt_family") {

  # Get cluster sizes
  cluster_sizes <- clustering_data[record_type == "H", .N, by = .(cluster_id, kt_family)]

  # Add singletons
  all_clusters <- clustering_data[record_type == "C", .(cluster_id, kt_family)]
  cluster_sizes <- merge(all_clusters, cluster_sizes, by = c("cluster_id", "kt_family"), all.x = TRUE)
  cluster_sizes[is.na(N), N := 0]

  # Filter by minimum size
  cluster_sizes <- cluster_sizes[N >= min_cluster_size]

  # Get average identity for each cluster
  cluster_identity <- clustering_data[
    record_type == "H" & !is.na(percent_identity),
    .(avg_identity = mean(percent_identity, na.rm = TRUE),
      min_identity = min(percent_identity, na.rm = TRUE)),
    by = .(cluster_id, kt_family)
  ]

  # Merge cluster sizes with identity data
  plot_data <- merge(cluster_sizes, cluster_identity, by = c("cluster_id", "kt_family"), all.x = TRUE)

  # Create plot
  p <- ggplot(plot_data[!is.na(avg_identity)],
              aes(x = avg_identity, y = N + 1)) +  # +1 to handle log scale with zeros
    geom_point(alpha = alpha, size = 2) +
    scale_y_log10(labels = function(x) x - 1) +  # Adjust labels back
    labs(
      title = "Cluster Size vs Average Identity",
      x = "Average Percent Identity (%)",
      y = "Cluster Size (log scale)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

  # Add coloring if requested
  if (color_by == "kt_family") {
    p <- p + aes(color = kt_family) +
      scale_color_viridis_d(name = "KT Family") +
      theme(legend.position = "right")
  }

  return(p)
}