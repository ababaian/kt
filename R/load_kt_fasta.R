#' Load and process KT FASTA sequences
#'
#' Generic function to read KT FASTA files (kt0, kt1, or similar formats) and
#' process them into a structured data.table with parsed header information.
#'
#' @param file_path Path to the KT FASTA file. If NULL, uses the default package path.
#' @param dataset_name Name of the dataset for reference ("kt0", "kt1", etc.).
#'   If NULL, inferred from filename.
#' @param package_data Logical, whether to use package data path (default FALSE)
#'
#' @return A data.table with columns:
#' \describe{
#'   \item{search_group}{Search group identifier (typically "kt")}
#'   \item{kt_family}{KT protein family (e.g., "M1", "M2", "SMK", "Mlus")}
#'   \item{organism}{Organism/species name}
#'   \item{accession}{Sequence accession number}
#'   \item{sequence}{Amino acid sequence}
#'   \item{seq_length}{Length of the sequence in amino acids}
#' }
#'
#' @examples
#' \dontrun{
#' # Load kt1 data
#' kt1_data <- load_kt_fasta("inst/extdata/kt1/kt1.fa", dataset_name = "kt1")
#'
#' # Load kt0 data
#' kt0_data <- load_kt_fasta("inst/extdata/kt0/kt0.preclust.fa", dataset_name = "kt0")
#'
#' # Basic exploration
#' table(kt1_data$kt_family)
#' summary(kt0_data$seq_length)
#' }
#'
#' @importFrom data.table data.table setDT
#' @export
load_kt_fasta <- function(file_path, dataset_name = NULL, package_data = FALSE) {

  # Determine file path and dataset name
  if (is.null(file_path)) {
    stop("file_path must be provided")
  }

  # Check if file exists
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  # Infer dataset name from filename if not provided
  if (is.null(dataset_name)) {
    dataset_name <- gsub("\\.(fa|fasta)$", "", basename(file_path))
    dataset_name <- gsub("^kt", "kt", dataset_name)  # Standardize kt prefix
  }

  # Read the FASTA file
  message("Reading FASTA file: ", basename(file_path))
  message("Dataset: ", dataset_name)

  # Read file line by line for processing
  lines <- readLines(file_path, warn = FALSE)

  # Initialize vectors to store data
  headers <- character()
  sequences <- character()
  current_seq <- character()

  # Parse FASTA format
  for (line in lines) {
    if (startsWith(line, ">")) {
      # If we have a current sequence, save it
      if (length(current_seq) > 0) {
        sequences <- c(sequences, paste(current_seq, collapse = ""))
        current_seq <- character()
      }
      # Store header (remove the ">" character)
      headers <- c(headers, substr(line, 2, nchar(line)))
    } else {
      # Accumulate sequence lines (remove any whitespace)
      clean_line <- gsub("[[:space:]]", "", line)
      if (nchar(clean_line) > 0) {  # Skip empty lines
        current_seq <- c(current_seq, clean_line)
      }
    }
  }

  # Don't forget the last sequence
  if (length(current_seq) > 0) {
    sequences <- c(sequences, paste(current_seq, collapse = ""))
  }

  # Check that we have equal numbers of headers and sequences
  if (length(headers) != length(sequences)) {
    stop("Mismatch between number of headers (", length(headers),
         ") and sequences (", length(sequences), ")")
  }

  message("Parsed ", length(headers), " sequences")

  # Parse headers (expected format: kt;family;organism;accession)
  header_parts <- strsplit(headers, ";", fixed = TRUE)

  # Check header format and extract components
  search_group <- character(length(header_parts))
  kt_family <- character(length(header_parts))
  organism <- character(length(header_parts))
  accession <- character(length(header_parts))

  for (i in seq_along(header_parts)) {
    parts <- header_parts[[i]]

    if (length(parts) < 4) {
      warning("Header ", i, " has fewer than 4 parts: ", headers[i])
      # Fill missing parts with NA
      parts <- c(parts, rep(NA, 4 - length(parts)))
    }

    search_group[i] <- parts[1] %||% NA_character_
    kt_family[i] <- parts[2] %||% NA_character_
    organism[i] <- parts[3] %||% NA_character_
    accession[i] <- parts[4] %||% NA_character_
  }

  # Create data.table
  dt <- data.table::data.table(
    search_group = search_group,
    kt_family = kt_family,
    organism = organism,
    accession = accession,
    sequence = sequences,
    seq_length = nchar(sequences)
  )

  message("Created data.table with ", nrow(dt), " rows and ", ncol(dt), " columns")

  # Print summary information
  message("Summary:")
  message("  Dataset: ", dataset_name)
  message("  Unique families: ", length(unique(dt$kt_family[!is.na(dt$kt_family)])))
  message("  Unique organisms: ", length(unique(dt$organism[!is.na(dt$organism)])))
  message("  Sequence length range: ", min(dt$seq_length), " - ", max(dt$seq_length), " aa")

  return(dt)
}


#' Helper operator for NULL default values
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}


#' Get summary statistics for KT FASTA datasets
#'
#' Provides summary statistics for loaded KT sequence data
#'
#' @param kt_data A data.table from load_kt_fasta()
#' @param group_by Character vector of columns to group by (e.g., "kt_family", "organism")
#' @param detailed Logical, whether to return detailed statistics (default FALSE)
#'
#' @return A list or data.table of summary statistics
#'
#' @examples
#' \dontrun{
#' kt_data <- load_kt_fasta("inst/extdata/kt1/kt1.fa")
#' summary_stats <- summarize_kt_data(kt_data)
#' print(summary_stats)
#'
#' # Group by family
#' family_stats <- summarize_kt_data(kt_data, group_by = "kt_family")
#' }
#'
#' @export
summarize_kt_data <- function(kt_data, group_by = NULL, detailed = FALSE) {

  if (!inherits(kt_data, "data.table")) {
    stop("Input must be a data.table from load_kt_fasta()")
  }

  if (is.null(group_by)) {
    # Overall summary
    stats <- list(
      total_sequences = nrow(kt_data),
      unique_families = length(unique(kt_data$kt_family[!is.na(kt_data$kt_family)])),
      unique_organisms = length(unique(kt_data$organism[!is.na(kt_data$organism)])),
      unique_accessions = length(unique(kt_data$accession[!is.na(kt_data$accession)])),
      sequence_length = list(
        min = min(kt_data$seq_length),
        max = max(kt_data$seq_length),
        mean = round(mean(kt_data$seq_length), 1),
        median = median(kt_data$seq_length)
      )
    )

    if (detailed) {
      # Family and organism distributions
      stats$family_distribution <- sort(table(kt_data$kt_family), decreasing = TRUE)
      stats$organism_distribution <- sort(table(kt_data$organism), decreasing = TRUE)[1:20]  # Top 20
    }

    return(stats)

  } else {
    # Grouped summary
    summary_cols <- c("sequence_count = .N",
                     "unique_organisms = data.table::uniqueN(organism)",
                     "mean_length = round(mean(seq_length), 1)",
                     "median_length = median(seq_length)",
                     "min_length = min(seq_length)",
                     "max_length = max(seq_length)")

    summary_expr <- paste(summary_cols, collapse = ", ")
    group_by_expr <- paste(group_by, collapse = ", ")

    # Use eval to build the expression
    cmd <- paste0("kt_data[, .(", summary_expr, "), by = .(", group_by_expr, ")][order(-sequence_count)]")
    result <- eval(parse(text = cmd))

    return(result)
  }
}


#' Compare two KT datasets
#'
#' Compare two KT datasets by family, organism, or sequence characteristics
#'
#' @param kt_data1 First data.table from load_kt_fasta()
#' @param kt_data2 Second data.table from load_kt_fasta()
#' @param dataset1_name Name for first dataset (default "Dataset1")
#' @param dataset2_name Name for second dataset (default "Dataset2")
#' @param compare_by What to compare: "family", "organism", "length", or "summary"
#' @param top_n Number of top categories to show (default 20)
#'
#' @return A data.table with comparison results
#'
#' @examples
#' \dontrun{
#' # Load datasets
#' kt0_data <- load_kt_fasta("inst/extdata/kt0/kt0.preclust.fa", "kt0")
#' kt1_data <- load_kt_fasta("inst/extdata/kt1/kt1.fa", "kt1")
#'
#' # Compare family distributions
#' family_comparison <- compare_kt_data(kt0_data, kt1_data, "kt0", "kt1", "family")
#' organism_comparison <- compare_kt_data(kt0_data, kt1_data, "kt0", "kt1", "organism")
#' }
#'
#' @export
compare_kt_data <- function(kt_data1, kt_data2, dataset1_name = "Dataset1", dataset2_name = "Dataset2",
                           compare_by = "family", top_n = 20) {

  if (!inherits(kt_data1, "data.table") || !inherits(kt_data2, "data.table")) {
    stop("Both inputs must be data.tables from load_kt_fasta()")
  }

  valid_comparisons <- c("family", "organism", "length", "summary")
  if (!compare_by %in% valid_comparisons) {
    stop("compare_by must be one of: ", paste(valid_comparisons, collapse = ", "))
  }

  switch(compare_by,
    "family" = {
      # Get family counts for each dataset
      families1 <- kt_data1[, .N, by = kt_family]
      families2 <- kt_data2[, .N, by = kt_family]

      # Get top families across both datasets if top_n specified
      if (!is.null(top_n)) {
        combined_families <- rbind(
          families1[, .(kt_family, total_n = N)],
          families2[, .(kt_family, total_n = N)]
        )[, .(total_n = sum(total_n)), by = kt_family][order(-total_n)]
        top_families <- head(combined_families$kt_family, top_n)
        families1 <- families1[kt_family %in% top_families]
        families2 <- families2[kt_family %in% top_families]
      }

      # Merge and format results
      result <- merge(families1, families2, by = "kt_family", all = TRUE, suffixes = c(paste0("_", dataset1_name), paste0("_", dataset2_name)))
      result[is.na(result)] <- 0
      result$total <- result[[paste0("N_", dataset1_name)]] + result[[paste0("N_", dataset2_name)]]
      result <- result[order(-total)]
    },
    "organism" = {
      # Get organism counts for each dataset
      organisms1 <- kt_data1[, .N, by = organism]
      organisms2 <- kt_data2[, .N, by = organism]

      # Get top organisms across both datasets if top_n specified
      if (!is.null(top_n)) {
        combined_organisms <- rbind(
          organisms1[, .(organism, total_n = N)],
          organisms2[, .(organism, total_n = N)]
        )[, .(total_n = sum(total_n)), by = organism][order(-total_n)]
        top_organisms <- head(combined_organisms$organism, top_n)
        organisms1 <- organisms1[organism %in% top_organisms]
        organisms2 <- organisms2[organism %in% top_organisms]
      }

      # Merge and format results
      result <- merge(organisms1, organisms2, by = "organism", all = TRUE, suffixes = c(paste0("_", dataset1_name), paste0("_", dataset2_name)))
      result[is.na(result)] <- 0
      result$total <- result[[paste0("N_", dataset1_name)]] + result[[paste0("N_", dataset2_name)]]
      result <- result[order(-total)]
    },
    "length" = {
      result <- data.table::data.table(
        dataset = c(dataset1_name, dataset2_name),
        sequences = c(nrow(kt_data1), nrow(kt_data2)),
        mean_length = c(round(mean(kt_data1$seq_length), 1), round(mean(kt_data2$seq_length), 1)),
        median_length = c(median(kt_data1$seq_length), median(kt_data2$seq_length)),
        min_length = c(min(kt_data1$seq_length), min(kt_data2$seq_length)),
        max_length = c(max(kt_data1$seq_length), max(kt_data2$seq_length))
      )
    },
    "summary" = {
      result <- data.table::data.table(
        dataset = c(dataset1_name, dataset2_name),
        sequences = c(nrow(kt_data1), nrow(kt_data2)),
        unique_families = c(data.table::uniqueN(kt_data1$kt_family), data.table::uniqueN(kt_data2$kt_family)),
        unique_organisms = c(data.table::uniqueN(kt_data1$organism), data.table::uniqueN(kt_data2$organism)),
        mean_length = c(round(mean(kt_data1$seq_length), 1), round(mean(kt_data2$seq_length), 1))
      )
    }
  )

  return(result)
}


#' Export KT data to various formats
#'
#' Export processed KT data to CSV, TSV, or other formats
#'
#' @param kt_data A data.table from load_kt_fasta()
#' @param file_path Output file path
#' @param format Output format: "csv", "tsv", or "rds" (default "csv")
#' @param include_sequence Logical, whether to include the full sequence column (default TRUE)
#'
#' @examples
#' \dontrun{
#' kt_data <- load_kt_fasta("inst/extdata/kt1/kt1.fa")
#'
#' # Export to CSV
#' export_kt_data(kt_data, "kt_processed.csv")
#'
#' # Export metadata only
#' export_kt_data(kt_data, "kt_metadata.tsv", format = "tsv", include_sequence = FALSE)
#' }
#'
#' @export
export_kt_data <- function(kt_data, file_path, format = "csv", include_sequence = TRUE) {

  if (!inherits(kt_data, "data.table")) {
    stop("Input must be a data.table from load_kt_fasta()")
  }

  # Prepare data for export
  if (include_sequence) {
    export_data <- kt_data
  } else {
    seq_col <- which(names(kt_data) == "sequence")
    if (length(seq_col) > 0) {
      export_data <- kt_data[, -seq_col, with = FALSE]
    } else {
      export_data <- kt_data
    }
  }

  # Export based on format
  switch(tolower(format),
    "csv" = {
      data.table::fwrite(export_data, file = file_path)
      message("Data exported to CSV: ", file_path)
    },
    "tsv" = {
      data.table::fwrite(export_data, file = file_path, sep = "\t")
      message("Data exported to TSV: ", file_path)
    },
    "rds" = {
      saveRDS(export_data, file = file_path)
      message("Data exported to RDS: ", file_path)
    },
    stop("Unsupported format: ", format, ". Use 'csv', 'tsv', or 'rds'")
  )

  invisible(file_path)
}