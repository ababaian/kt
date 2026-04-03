#' Load UCLUST clustering output
#'
#' Reads UCLUST clustering output (.uc files) and alignment data into a structured data.frame
#' UCLUST format has 10 tab-separated columns with clustering and alignment information.
#'
#' @param file_path Path to the UCLUST .uc file
#' @param include_seeds Logical, whether to include seed/representative sequences (default TRUE)
#' @param include_hits Logical, whether to include hit/member sequences (default TRUE)
#'
#' @return A data.table with columns:
#' \describe{
#'   \item{record_type}{Record type: 'S' for seed/representative, 'H' for hit/member, 'C' for centroid}
#'   \item{cluster_id}{Cluster identifier (0-based)}
#'   \item{sequence_length}{Length of the sequence}
#'   \item{percent_identity}{Percent identity to cluster representative (NA for seeds)}
#'   \item{strand}{Strand orientation ('.' for protein)}
#'   \item{query_start}{Query alignment start position (NA for seeds)}
#'   \item{target_start}{Target alignment start position (NA for seeds)}
#'   \item{alignment}{Alignment string (NA for seeds)}
#'   \item{query_label}{Query sequence identifier}
#'   \item{target_label}{Target sequence identifier (same as query for seeds)}
#'   \item{search_group}{Parsed search group from sequence identifier}
#'   \item{kt_family}{Parsed KT family from sequence identifier}
#'   \item{organism}{Parsed organism from sequence identifier}
#'   \item{accession}{Parsed accession from sequence identifier}
#' }
#'
#' @examples
#' \dontrun{
#' # Load UCLUST clustering output
#' clustering_data <- load_uclust("inst/extdata/kt1/kt1.id90.uc")
#'
#' # Get only cluster representatives
#' seeds_only <- load_uclust("inst/extdata/kt1/kt1.id90.uc", include_hits = FALSE)
#'
#' # Basic exploration
#' table(clustering_data$record_type)
#' table(clustering_data$kt_family)
#' hist(clustering_data$percent_identity, na.rm = TRUE)
#' }
#'
#' @importFrom data.table data.table fread uniqueN
#' @export
load_uclust <- function(file_path, include_seeds = TRUE, include_hits = TRUE) {

  # Check if file exists
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  if (!include_seeds && !include_hits) {
    stop("At least one of include_seeds or include_hits must be TRUE")
  }

  message("Reading UCLUST file: ", basename(file_path))

  # Read the UCLUST file
  # UCLUST format: 10 tab-separated columns
  col_names <- c("record_type", "cluster_id", "sequence_length", "percent_identity",
                 "strand", "query_start", "target_start", "alignment",
                 "query_label", "target_label")

  dt <- data.table::fread(file_path,
                         header = FALSE,
                         sep = "\t",
                         col.names = col_names,
                         colClasses = list(character = c(1, 4:10),
                                         integer = c(2, 3)))

  message("Read ", nrow(dt), " records")

  # Filter by record type if requested
  if (!include_seeds) {
    dt <- dt[record_type != "S"]
  }
  if (!include_hits) {
    dt <- dt[record_type != "H"]
  }

  message("Kept ", nrow(dt), " records after filtering")

  # Process percent identity (convert "*" to NA and to numeric)
  dt$percent_identity[dt$percent_identity == "*"] <- NA
  dt$percent_identity <- as.numeric(dt$percent_identity)

  # Process alignment positions (convert "*" to NA and to numeric)
  dt$query_start[dt$query_start == "*"] <- NA
  dt$query_start <- as.numeric(dt$query_start)

  dt$target_start[dt$target_start == "*"] <- NA
  dt$target_start <- as.numeric(dt$target_start)

  # Process alignment string (convert "*" to NA)
  dt$alignment[dt$alignment == "*"] <- NA

  # Process target label (convert "*" to NA)
  dt$target_label[dt$target_label == "*"] <- NA

  # Parse sequence identifiers (expected format: kt;family;organism;accession)
  # Use query_label for parsing since it's always present
  message("Parsing sequence identifiers...")

  header_parts <- strsplit(dt$query_label, ";", fixed = TRUE)

  # Initialize parsed columns
  search_group <- character(nrow(dt))
  kt_family <- character(nrow(dt))
  organism <- character(nrow(dt))
  accession <- character(nrow(dt))

  # Helper operator for NULL default values
  `%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
  }

  for (i in seq_along(header_parts)) {
    parts <- header_parts[[i]]

    if (length(parts) < 4) {
      warning("Query label ", i, " has fewer than 4 parts: ", dt$query_label[i])
      # Fill missing parts with NA
      parts <- c(parts, rep(NA, 4 - length(parts)))
    }

    search_group[i] <- parts[1] %||% NA_character_
    kt_family[i] <- parts[2] %||% NA_character_
    organism[i] <- parts[3] %||% NA_character_
    accession[i] <- parts[4] %||% NA_character_
  }

  # Add parsed columns
  dt$search_group <- search_group
  dt$kt_family <- kt_family
  dt$organism <- organism
  dt$accession <- accession

  message("Summary:")
  message("  Total records: ", nrow(dt))
  message("  Seeds (representatives): ", sum(dt$record_type == "S"))
  message("  Hits (members): ", sum(dt$record_type == "H"))
  message("  Centroids: ", sum(dt$record_type == "C"))
  message("  Unique clusters: ", data.table::uniqueN(dt$cluster_id))
  message("  Unique families: ", length(unique(dt$kt_family[!is.na(dt$kt_family)])))
  message("  Identity range: ",
          ifelse(any(!is.na(dt$percent_identity)),
                 paste(round(min(dt$percent_identity, na.rm = TRUE), 1), "-",
                       round(max(dt$percent_identity, na.rm = TRUE), 1), "%"),
                 "N/A (seeds only)"))

  return(dt)
}


#' Get clustering summary statistics
#'
#' Provides summary statistics for UCLUST clustering results
#'
#' @param uclust_data A data.table from load_uclust()
#' @param by_family Logical, whether to group statistics by KT family (default FALSE)
#' @param by_cluster Logical, whether to include per-cluster statistics (default FALSE)
#'
#' @return A list or data.table of clustering summary statistics
#'
#' @examples
#' \dontrun{
#' clustering_data <- load_uclust("inst/extdata/kt1/kt1.id90.uc")
#'
#' # Overall clustering summary
#' cluster_summary <- summarize_clustering(clustering_data)
#' print(cluster_summary)
#'
#' # Summary by KT family
#' family_summary <- summarize_clustering(clustering_data, by_family = TRUE)
#' }
#'
#' @export
summarize_clustering <- function(uclust_data, by_family = FALSE, by_cluster = FALSE) {

  if (!inherits(uclust_data, "data.table")) {
    stop("Input must be a data.table from load_uclust()")
  }

  if (!by_family && !by_cluster) {
    # Overall summary
    total_sequences <- nrow(uclust_data)
    n_seeds <- sum(uclust_data$record_type == "S")
    n_hits <- sum(uclust_data$record_type == "H")
    n_centroids <- sum(uclust_data$record_type == "C")
    n_clusters <- data.table::uniqueN(uclust_data$cluster_id)

    # Cluster size statistics (only meaningful if we have hits)
    if (n_hits > 0) {
      cluster_sizes <- uclust_data[, .N, by = cluster_id]$N
      cluster_stats <- list(
        min = min(cluster_sizes),
        max = max(cluster_sizes),
        mean = round(mean(cluster_sizes), 1),
        median = median(cluster_sizes)
      )
    } else {
      cluster_stats <- list(min = 1, max = 1, mean = 1, median = 1)
    }

    # Identity statistics (only for hits)
    if (any(!is.na(uclust_data$percent_identity))) {
      identity_stats <- list(
        min = round(min(uclust_data$percent_identity, na.rm = TRUE), 1),
        max = round(max(uclust_data$percent_identity, na.rm = TRUE), 1),
        mean = round(mean(uclust_data$percent_identity, na.rm = TRUE), 1),
        median = round(median(uclust_data$percent_identity, na.rm = TRUE), 1)
      )
    } else {
      identity_stats <- list(min = NA, max = NA, mean = NA, median = NA)
    }

    summary_list <- list(
      total_sequences = total_sequences,
      cluster_representatives = n_seeds,
      cluster_members = n_hits,
      cluster_centroids = n_centroids,
      total_clusters = n_clusters,
      cluster_size_stats = cluster_stats,
      identity_stats = identity_stats,
      compression_ratio = round((n_seeds + n_centroids) / total_sequences, 3)
    )

    return(summary_list)

  } else if (by_family) {
    # Summary by KT family
    family_summary <- uclust_data[, .(
      total_sequences = .N,
      seeds = sum(record_type == "S"),
      hits = sum(record_type == "H"),
      centroids = sum(record_type == "C"),
      unique_clusters = data.table::uniqueN(cluster_id),
      mean_identity = round(mean(percent_identity, na.rm = TRUE), 1),
      unique_organisms = data.table::uniqueN(organism)
    ), by = kt_family][order(-total_sequences)]

    # Add compression ratio
    family_summary$compression_ratio <- round((family_summary$seeds + family_summary$centroids) / family_summary$total_sequences, 3)

    return(family_summary)

  } else if (by_cluster) {
    # Per-cluster statistics
    cluster_summary <- uclust_data[, .(
      cluster_size = .N,
      has_hits = any(record_type == "H"),
      families_in_cluster = data.table::uniqueN(kt_family),
      organisms_in_cluster = data.table::uniqueN(organism),
      representative = query_label[record_type %in% c("S", "C")][1],
      mean_identity = round(mean(percent_identity, na.rm = TRUE), 1)
    ), by = cluster_id][order(-cluster_size)]

    return(cluster_summary)
  }
}