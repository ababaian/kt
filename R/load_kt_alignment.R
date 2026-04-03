#' Load USEARCH allpairs alignment output
#'
#' Reads USEARCH allpairs_local output in BLAST6 format into a structured data.table
#' with parsed sequence identifiers and alignment statistics.
#'
#' @param file_path Path to the USEARCH alignment file (.aln or similar)
#' @param include_self Logical, whether to include self-alignments (default FALSE)
#' @param min_identity Numeric, minimum percent identity to include (default 0)
#' @param min_length Numeric, minimum alignment length to include (default 0)
#' @param max_evalue Numeric, maximum e-value to include (default 1)
#'
#' @return A data.table with columns:
#' \describe{
#'   \item{query_id}{Query sequence identifier}
#'   \item{subject_id}{Subject sequence identifier}
#'   \item{percent_identity}{Percentage of identical matches}
#'   \item{alignment_length}{Length of alignment}
#'   \item{mismatches}{Number of mismatches}
#'   \item{gap_opens}{Number of gap openings}
#'   \item{query_start}{Start of alignment in query}
#'   \item{query_end}{End of alignment in query}
#'   \item{subject_start}{Start of alignment in subject}
#'   \item{subject_end}{End of alignment in subject}
#'   \item{evalue}{Expect value}
#'   \item{bitscore}{Bit score}
#'   \item{query_search_group}{Parsed search group from query identifier}
#'   \item{query_kt_family}{Parsed KT family from query identifier}
#'   \item{query_organism}{Parsed organism from query identifier}
#'   \item{query_accession}{Parsed accession from query identifier}
#'   \item{subject_search_group}{Parsed search group from subject identifier}
#'   \item{subject_kt_family}{Parsed KT family from subject identifier}
#'   \item{subject_organism}{Parsed organism from subject identifier}
#'   \item{subject_accession}{Parsed accession from subject identifier}
#' }
#'
#' @examples
#' \dontrun{
#' # Load all alignments
#' alignments <- load_kt_alignment("inst/extdata/kt1/kt1.aln")
#'
#' # Load high-quality alignments only
#' hq_alignments <- load_kt_alignment("inst/extdata/kt1/kt1.aln",
#'                                   min_identity = 30,
#'                                   min_length = 50,
#'                                   max_evalue = 1e-5)
#'
#' # Basic exploration
#' summary(alignments$percent_identity)
#' table(alignments$query_kt_family, alignments$subject_kt_family)
#' }
#'
#' @importFrom data.table data.table fread copy
#' @export
load_kt_alignment <- function(file_path, include_self = FALSE, min_identity = 0,
                             min_length = 0, max_evalue = 1) {

  # Check if file exists
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  message("Reading USEARCH alignment file: ", basename(file_path))

  # Define BLAST6 column names
  col_names <- c("query_id", "subject_id", "percent_identity", "alignment_length",
                 "mismatches", "gap_opens", "query_start", "query_end",
                 "subject_start", "subject_end", "evalue", "bitscore")

  # Read the alignment file
  dt <- data.table::fread(file_path,
                         header = FALSE,
                         sep = "\t",
                         col.names = col_names,
                         colClasses = list(character = c(1, 2),
                                         numeric = c(3, 11, 12),
                                         integer = c(4:10)))

  message("Read ", nrow(dt), " alignment records")

  # Apply filters
  original_count <- nrow(dt)

  # Filter by identity
  if (min_identity > 0) {
    dt <- dt[percent_identity >= min_identity]
    message("Filtered by identity >= ", min_identity, "%: ", nrow(dt), " records remain")
  }

  # Filter by alignment length
  if (min_length > 0) {
    dt <- dt[alignment_length >= min_length]
    message("Filtered by length >= ", min_length, ": ", nrow(dt), " records remain")
  }

  # Filter by e-value
  if (max_evalue < 1) {
    dt <- dt[evalue <= max_evalue]
    message("Filtered by e-value <= ", max_evalue, ": ", nrow(dt), " records remain")
  }

  # Filter self-alignments if requested
  if (!include_self) {
    dt <- dt[query_id != subject_id]
    message("Removed self-alignments: ", nrow(dt), " records remain")
  }

  # Parse sequence identifiers
  message("Parsing sequence identifiers...")

  # Helper operator for NULL default values
  `%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
  }

  # Parse query identifiers
  query_parts <- strsplit(dt$query_id, ";", fixed = TRUE)
  query_search_group <- character(nrow(dt))
  query_kt_family <- character(nrow(dt))
  query_organism <- character(nrow(dt))
  query_accession <- character(nrow(dt))

  for (i in seq_along(query_parts)) {
    parts <- query_parts[[i]]
    if (length(parts) < 4) {
      parts <- c(parts, rep(NA, 4 - length(parts)))
    }
    query_search_group[i] <- parts[1] %||% NA_character_
    query_kt_family[i] <- parts[2] %||% NA_character_
    query_organism[i] <- parts[3] %||% NA_character_
    query_accession[i] <- parts[4] %||% NA_character_
  }

  # Parse subject identifiers
  subject_parts <- strsplit(dt$subject_id, ";", fixed = TRUE)
  subject_search_group <- character(nrow(dt))
  subject_kt_family <- character(nrow(dt))
  subject_organism <- character(nrow(dt))
  subject_accession <- character(nrow(dt))

  for (i in seq_along(subject_parts)) {
    parts <- subject_parts[[i]]
    if (length(parts) < 4) {
      parts <- c(parts, rep(NA, 4 - length(parts)))
    }
    subject_search_group[i] <- parts[1] %||% NA_character_
    subject_kt_family[i] <- parts[2] %||% NA_character_
    subject_organism[i] <- parts[3] %||% NA_character_
    subject_accession[i] <- parts[4] %||% NA_character_
  }

  # Add parsed columns
  dt[, ':='(
    query_search_group = query_search_group,
    query_kt_family = query_kt_family,
    query_organism = query_organism,
    query_accession = query_accession,
    subject_search_group = subject_search_group,
    subject_kt_family = subject_kt_family,
    subject_organism = subject_organism,
    subject_accession = subject_accession
  )]

  message("Summary:")
  message("  Total alignments: ", nrow(dt))
  message("  Filtered out: ", original_count - nrow(dt), " records")
  message("  Unique query families: ", length(unique(dt$query_kt_family[!is.na(dt$query_kt_family)])))
  message("  Unique subject families: ", length(unique(dt$subject_kt_family[!is.na(dt$subject_kt_family)])))
  message("  Identity range: ", round(min(dt$percent_identity), 1), "% - ",
          round(max(dt$percent_identity), 1), "%")
  message("  E-value range: ", format(min(dt$evalue), scientific = TRUE), " - ",
          format(max(dt$evalue), scientific = TRUE))

  return(dt)
}


#' Summarize alignment statistics
#'
#' Provides comprehensive statistics for USEARCH alignment results
#'
#' @param alignment_data A data.table from load_kt_alignment()
#' @param by_family Logical, whether to group by KT family pairs (default FALSE)
#' @param by_organism Logical, whether to group by organism pairs (default FALSE)
#' @param top_n Number of top categories to show for grouped summaries (default 20)
#'
#' @return A list or data.table of alignment summary statistics
#'
#' @examples
#' \dontrun{
#' alignments <- load_kt_alignment("inst/extdata/kt1/kt1.aln")
#'
#' # Overall summary
#' overall_stats <- summarize_alignment(alignments)
#' print(overall_stats)
#'
#' # Summary by family pairs
#' family_stats <- summarize_alignment(alignments, by_family = TRUE)
#' }
#'
#' @export
summarize_alignment <- function(alignment_data, by_family = FALSE, by_organism = FALSE, top_n = 20) {

  if (!inherits(alignment_data, "data.table")) {
    stop("Input must be a data.table from load_kt_alignment()")
  }

  if (!by_family && !by_organism) {
    # Overall summary
    stats <- list(
      total_alignments = nrow(alignment_data),
      unique_queries = length(unique(alignment_data$query_id)),
      unique_subjects = length(unique(alignment_data$subject_id)),

      identity_stats = list(
        min = round(min(alignment_data$percent_identity), 1),
        max = round(max(alignment_data$percent_identity), 1),
        mean = round(mean(alignment_data$percent_identity), 1),
        median = round(median(alignment_data$percent_identity), 1)
      ),

      length_stats = list(
        min = min(alignment_data$alignment_length),
        max = max(alignment_data$alignment_length),
        mean = round(mean(alignment_data$alignment_length), 1),
        median = median(alignment_data$alignment_length)
      ),

      evalue_stats = list(
        min = min(alignment_data$evalue),
        max = max(alignment_data$evalue),
        mean = mean(alignment_data$evalue),
        median = median(alignment_data$evalue)
      ),

      bitscore_stats = list(
        min = round(min(alignment_data$bitscore), 1),
        max = round(max(alignment_data$bitscore), 1),
        mean = round(mean(alignment_data$bitscore), 1),
        median = round(median(alignment_data$bitscore), 1)
      )
    )

    return(stats)

  } else if (by_family) {
    # Summary by KT family pairs
    family_summary <- alignment_data[, .(
      alignment_count = .N,
      unique_queries = data.table::uniqueN(query_id),
      unique_subjects = data.table::uniqueN(subject_id),
      mean_identity = round(mean(percent_identity), 1),
      median_identity = round(median(percent_identity), 1),
      max_identity = round(max(percent_identity), 1),
      mean_length = round(mean(alignment_length), 1),
      mean_bitscore = round(mean(bitscore), 1)
    ), by = .(query_kt_family, subject_kt_family)][order(-alignment_count)]

    # Add family pair column for easier reading
    family_summary[, family_pair := paste(query_kt_family, "vs", subject_kt_family)]

    # Return top N if specified
    if (!is.null(top_n) && nrow(family_summary) > top_n) {
      family_summary <- head(family_summary, top_n)
    }

    return(family_summary)

  } else if (by_organism) {
    # Summary by organism pairs (simplified to avoid too much detail)
    organism_summary <- alignment_data[, .(
      alignment_count = .N,
      mean_identity = round(mean(percent_identity), 1),
      max_identity = round(max(percent_identity), 1),
      mean_bitscore = round(mean(bitscore), 1)
    ), by = .(query_organism, subject_organism)][order(-alignment_count)]

    # Return top N if specified
    if (!is.null(top_n) && nrow(organism_summary) > top_n) {
      organism_summary <- head(organism_summary, top_n)
    }

    return(organism_summary)
  }
}


#' Filter alignments by various criteria
#'
#' Filter alignment data based on multiple criteria including identity, length,
#' family relationships, and statistical thresholds
#'
#' @param alignment_data A data.table from load_kt_alignment()
#' @param min_identity Numeric, minimum percent identity (default NULL)
#' @param max_identity Numeric, maximum percent identity (default NULL)
#' @param min_length Numeric, minimum alignment length (default NULL)
#' @param max_evalue Numeric, maximum e-value (default NULL)
#' @param min_bitscore Numeric, minimum bit score (default NULL)
#' @param same_family Logical, keep only same-family alignments (default NULL)
#' @param different_family Logical, keep only different-family alignments (default NULL)
#' @param same_organism Logical, keep only same-organism alignments (default NULL)
#' @param different_organism Logical, keep only different-organism alignments (default NULL)
#' @param families Character vector, keep only specified families (default NULL)
#'
#' @return Filtered data.table
#'
#' @examples
#' \dontrun{
#' alignments <- load_kt_alignment("inst/extdata/kt1/kt1.aln")
#'
#' # High-quality alignments
#' hq_alignments <- filter_alignments(alignments,
#'                                   min_identity = 50,
#'                                   min_length = 100,
#'                                   max_evalue = 1e-10)
#'
#' # Cross-family alignments
#' cross_family <- filter_alignments(alignments, different_family = TRUE)
#'
#' # Focus on specific families
#' m1_m2 <- filter_alignments(alignments, families = c("M1", "M2"))
#' }
#'
#' @export
filter_alignments <- function(alignment_data, min_identity = NULL, max_identity = NULL,
                             min_length = NULL, max_evalue = NULL, min_bitscore = NULL,
                             same_family = NULL, different_family = NULL,
                             same_organism = NULL, different_organism = NULL,
                             families = NULL) {

  if (!inherits(alignment_data, "data.table")) {
    stop("Input must be a data.table from load_kt_alignment()")
  }

  dt <- copy(alignment_data)  # Work on a copy
  original_count <- nrow(dt)

  # Apply identity filters
  if (!is.null(min_identity)) {
    dt <- dt[percent_identity >= min_identity]
    message("Filtered by min identity ", min_identity, "%: ", nrow(dt), " records remain")
  }

  if (!is.null(max_identity)) {
    dt <- dt[percent_identity <= max_identity]
    message("Filtered by max identity ", max_identity, "%: ", nrow(dt), " records remain")
  }

  # Apply length filter
  if (!is.null(min_length)) {
    dt <- dt[alignment_length >= min_length]
    message("Filtered by min length ", min_length, ": ", nrow(dt), " records remain")
  }

  # Apply e-value filter
  if (!is.null(max_evalue)) {
    dt <- dt[evalue <= max_evalue]
    message("Filtered by max e-value ", max_evalue, ": ", nrow(dt), " records remain")
  }

  # Apply bitscore filter
  if (!is.null(min_bitscore)) {
    dt <- dt[bitscore >= min_bitscore]
    message("Filtered by min bitscore ", min_bitscore, ": ", nrow(dt), " records remain")
  }

  # Apply family relationship filters
  if (!is.null(same_family) && same_family) {
    dt <- dt[query_kt_family == subject_kt_family]
    message("Filtered for same-family alignments: ", nrow(dt), " records remain")
  }

  if (!is.null(different_family) && different_family) {
    dt <- dt[query_kt_family != subject_kt_family]
    message("Filtered for different-family alignments: ", nrow(dt), " records remain")
  }

  # Apply organism relationship filters
  if (!is.null(same_organism) && same_organism) {
    dt <- dt[query_organism == subject_organism]
    message("Filtered for same-organism alignments: ", nrow(dt), " records remain")
  }

  if (!is.null(different_organism) && different_organism) {
    dt <- dt[query_organism != subject_organism]
    message("Filtered for different-organism alignments: ", nrow(dt), " records remain")
  }

  # Apply family-specific filter
  if (!is.null(families)) {
    dt <- dt[query_kt_family %in% families | subject_kt_family %in% families]
    message("Filtered for families [", paste(families, collapse = ", "), "]: ",
            nrow(dt), " records remain")
  }

  message("Total filtered: ", original_count - nrow(dt), " records removed")

  return(dt)
}


#' Find best reciprocal alignments
#'
#' Identify reciprocal best hits (RBH) and best alignments for each query-subject pair
#'
#' @param alignment_data A data.table from load_kt_alignment()
#' @param method Character, method for determining best hits: "bitscore", "evalue", or "identity" (default "bitscore")
#' @param reciprocal Logical, whether to find reciprocal best hits (default TRUE)
#'
#' @return A data.table with best alignments, optionally filtered for reciprocal hits
#'
#' @examples
#' \dontrun{
#' alignments <- load_kt_alignment("inst/extdata/kt1/kt1.aln")
#'
#' # Find reciprocal best hits by bitscore
#' rbh <- find_best_alignments(alignments, reciprocal = TRUE)
#'
#' # Find best hits by identity (not necessarily reciprocal)
#' best_identity <- find_best_alignments(alignments, method = "identity", reciprocal = FALSE)
#' }
#'
#' @export
find_best_alignments <- function(alignment_data, method = "bitscore", reciprocal = TRUE) {

  if (!inherits(alignment_data, "data.table")) {
    stop("Input must be a data.table from load_kt_alignment()")
  }

  valid_methods <- c("bitscore", "evalue", "identity")
  if (!method %in% valid_methods) {
    stop("method must be one of: ", paste(valid_methods, collapse = ", "))
  }

  dt <- copy(alignment_data)

  # Determine sorting criteria based on method
  if (method == "bitscore") {
    best_forward <- dt[order(query_id, -bitscore), .SD[1], by = query_id]
    best_reverse <- dt[order(subject_id, -bitscore), .SD[1], by = subject_id]
  } else if (method == "evalue") {
    best_forward <- dt[order(query_id, evalue), .SD[1], by = query_id]
    best_reverse <- dt[order(subject_id, evalue), .SD[1], by = subject_id]
  } else if (method == "identity") {
    best_forward <- dt[order(query_id, -percent_identity), .SD[1], by = query_id]
    best_reverse <- dt[order(subject_id, -percent_identity), .SD[1], by = subject_id]
  }

  if (!reciprocal) {
    message("Found ", nrow(best_forward), " best forward alignments by ", method)
    return(best_forward)
  }

  # Find reciprocal best hits
  # Create standardized query-subject pairs for both directions
  best_forward[, pair := paste(pmin(query_id, subject_id), pmax(query_id, subject_id), sep = "___")]
  best_reverse[, pair := paste(pmin(query_id, subject_id), pmax(query_id, subject_id), sep = "___")]

  # Find pairs that appear in both directions
  reciprocal_pairs <- intersect(best_forward$pair, best_reverse$pair)

  # Get the reciprocal best hits
  rbh <- best_forward[pair %in% reciprocal_pairs]

  # Clean up temporary column
  rbh[, pair := NULL]

  message("Found ", length(reciprocal_pairs), " reciprocal best hit pairs by ", method)
  message("Total RBH alignments: ", nrow(rbh))

  return(rbh)
}


#' Export alignment data to various formats
#'
#' Export processed alignment data to CSV, TSV, or other formats
#'
#' @param alignment_data A data.table from load_kt_alignment()
#' @param file_path Output file path
#' @param format Output format: "csv", "tsv", or "rds" (default "csv")
#' @param include_parsed Logical, whether to include parsed identifier columns (default TRUE)
#'
#' @examples
#' \dontrun{
#' alignments <- load_kt_alignment("inst/extdata/kt1/kt1.aln")
#'
#' # Export to CSV
#' export_alignment_data(alignments, "kt_alignments.csv")
#'
#' # Export core alignment data only
#' export_alignment_data(alignments, "kt_alignments_core.tsv",
#'                      format = "tsv", include_parsed = FALSE)
#' }
#'
#' @export
export_alignment_data <- function(alignment_data, file_path, format = "csv", include_parsed = TRUE) {

  if (!inherits(alignment_data, "data.table")) {
    stop("Input must be a data.table from load_kt_alignment()")
  }

  # Prepare data for export
  if (include_parsed) {
    export_data <- alignment_data
  } else {
    # Keep only core BLAST6 columns
    core_cols <- c("query_id", "subject_id", "percent_identity", "alignment_length",
                   "mismatches", "gap_opens", "query_start", "query_end",
                   "subject_start", "subject_end", "evalue", "bitscore")
    export_data <- alignment_data[, ..core_cols]
  }

  # Export based on format
  switch(tolower(format),
    "csv" = {
      data.table::fwrite(export_data, file = file_path)
      message("Alignment data exported to CSV: ", file_path)
    },
    "tsv" = {
      data.table::fwrite(export_data, file = file_path, sep = "\t")
      message("Alignment data exported to TSV: ", file_path)
    },
    "rds" = {
      saveRDS(export_data, file = file_path)
      message("Alignment data exported to RDS: ", file_path)
    },
    stop("Unsupported format: ", format, ". Use 'csv', 'tsv', or 'rds'")
  )

  invisible(file_path)
}