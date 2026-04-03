#' Plot Family Count Comparison
#'
#' Creates a bar graph comparing the total number of sequences for each
#' KT family across multiple datasets (e.g., KT1 vs KT0).
#'
#' @param ... Named data.frames from load_kt_fasta(). Names will be used as dataset labels.
#' @param title Character, plot title. Default: "Family Count Comparison"
#' @param top_n Integer, number of top families to show. Default: 15
#' @param color_palette Character vector or function for colors. Default: viridis colors
#' @param show_counts Logical, whether to show count labels on bars. Default: TRUE
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' kt1_data <- load_kt_fasta("kt1.fa", "kt1")
#' kt0_data <- load_kt_fasta("kt0.fa", "kt0")
#' plot_family_count(KT1 = kt1_data, KT0 = kt0_data)
#'
#' # Custom styling
#' plot_family_count(KT1 = kt1_data, KT0 = kt0_data,
#'                  title = "KT Family Distribution",
#'                  top_n = 10,
#'                  color_palette = c("#1f77b4", "#ff7f0e"))
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_col facet_wrap theme_minimal
#' @importFrom ggplot2 labs scale_fill_manual geom_text position_dodge
#' @importFrom ggplot2 theme element_text
#' @importFrom dplyr bind_rows count arrange desc slice_head
#' @importFrom viridis viridis
#' @export
plot_family_count <- function(...,
                             title = "Family Count Comparison",
                             top_n = 15,
                             color_palette = NULL,
                             show_counts = TRUE) {

  # Get datasets from ...
  datasets <- list(...)

  if (length(datasets) == 0) {
    stop("Please provide at least one dataset")
  }

  # Extract dataset names
  dataset_names <- names(datasets)
  if (is.null(dataset_names) || any(dataset_names == "")) {
    dataset_names <- paste0("Dataset_", seq_along(datasets))
    names(datasets) <- dataset_names
  }

  # Combine all datasets with dataset labels
  combined_data <- bind_rows(datasets, .id = "dataset")

  # Count families per dataset
  family_counts <- combined_data %>%
    count(dataset, kt_family, name = "count") %>%
    arrange(desc(count))

  # Get top families overall
  top_families <- family_counts %>%
    count(kt_family, wt = count, name = "total_count") %>%
    arrange(desc(total_count)) %>%
    slice_head(n = top_n) %>%
    pull(kt_family)

  # Filter to top families
  plot_data <- family_counts %>%
    filter(kt_family %in% top_families) %>%
    # Ensure factor levels are ordered by total count
    mutate(kt_family = factor(kt_family, levels = rev(top_families)))

  # Set up colors
  if (is.null(color_palette)) {
    colors <- viridis(length(dataset_names))
  } else if (is.function(color_palette)) {
    colors <- color_palette(length(dataset_names))
  } else {
    colors <- rep(color_palette, length.out = length(dataset_names))
  }
  names(colors) <- dataset_names

  # Create plot
  p <- ggplot(plot_data, aes(x = kt_family, y = count, fill = dataset)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = colors) +
    labs(
      title = title,
      x = "KT Family",
      y = "Number of Sequences",
      fill = "Dataset"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    coord_flip()

  # Add count labels if requested
  if (show_counts) {
    p <- p + geom_text(
      aes(label = count),
      position = position_dodge(width = 0.9),
      hjust = -0.1,
      size = 3
    )
  }

  return(p)
}

#' Plot Organism Count Comparison
#'
#' Creates a bar graph comparing the total number of sequences for each
#' organism across multiple datasets (e.g., KT1 vs KT0).
#'
#' @param ... Named data.frames from load_kt_fasta(). Names will be used as dataset labels.
#' @param title Character, plot title. Default: "Organism Count Comparison"
#' @param top_n Integer, number of top organisms to show. Default: 15
#' @param color_palette Character vector or function for colors. Default: viridis colors
#' @param show_counts Logical, whether to show count labels on bars. Default: TRUE
#' @param shorten_names Logical, whether to shorten organism names. Default: TRUE
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' kt1_data <- load_kt_fasta("kt1.fa", "kt1")
#' kt0_data <- load_kt_fasta("kt0.fa", "kt0")
#' plot_organism_count(KT1 = kt1_data, KT0 = kt0_data)
#'
#' # Show full organism names
#' plot_organism_count(KT1 = kt1_data, KT0 = kt0_data,
#'                    shorten_names = FALSE,
#'                    top_n = 10)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_col theme_minimal labs scale_fill_manual
#' @importFrom ggplot2 geom_text position_dodge theme element_text coord_flip
#' @importFrom dplyr bind_rows count arrange desc slice_head mutate
#' @importFrom stringr str_trunc str_replace_all
#' @importFrom viridis viridis
#' @export
plot_organism_count <- function(...,
                               title = "Organism Count Comparison",
                               top_n = 15,
                               color_palette = NULL,
                               show_counts = TRUE,
                               shorten_names = TRUE) {

  # Get datasets from ...
  datasets <- list(...)

  if (length(datasets) == 0) {
    stop("Please provide at least one dataset")
  }

  # Extract dataset names
  dataset_names <- names(datasets)
  if (is.null(dataset_names) || any(dataset_names == "")) {
    dataset_names <- paste0("Dataset_", seq_along(datasets))
    names(datasets) <- dataset_names
  }

  # Combine all datasets with dataset labels
  combined_data <- bind_rows(datasets, .id = "dataset")

  # Count organisms per dataset
  organism_counts <- combined_data %>%
    count(dataset, organism, name = "count") %>%
    arrange(desc(count))

  # Get top organisms overall
  top_organisms <- organism_counts %>%
    count(organism, wt = count, name = "total_count") %>%
    arrange(desc(total_count)) %>%
    slice_head(n = top_n) %>%
    pull(organism)

  # Filter to top organisms and optionally shorten names
  plot_data <- organism_counts %>%
    filter(organism %in% top_organisms) %>%
    mutate(
      organism_display = if (shorten_names) {
        # Shorten organism names for better display
        str_replace_all(organism, "_", " ") %>%
          str_trunc(25, "right", "...")
      } else {
        organism
      }
    ) %>%
    # Order by total count
    mutate(organism_display = factor(organism_display,
                                   levels = rev(unique(organism_display[order(match(organism, top_organisms))]))))

  # Set up colors
  if (is.null(color_palette)) {
    colors <- viridis(length(dataset_names))
  } else if (is.function(color_palette)) {
    colors <- color_palette(length(dataset_names))
  } else {
    colors <- rep(color_palette, length.out = length(dataset_names))
  }
  names(colors) <- dataset_names

  # Create plot
  p <- ggplot(plot_data, aes(x = organism_display, y = count, fill = dataset)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = colors) +
    labs(
      title = title,
      x = "Organism",
      y = "Number of Sequences",
      fill = "Dataset"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    coord_flip()

  # Add count labels if requested
  if (show_counts) {
    p <- p + geom_text(
      aes(label = count),
      position = position_dodge(width = 0.9),
      hjust = -0.1,
      size = 3
    )
  }

  return(p)
}

#' Plot Family Length Distribution
#'
#' Creates a histogram showing sequence length distribution with options
#' to display by family or all families combined. Supports interactive
#' family selection for detailed analysis.
#'
#' @param data A data.frame from load_kt_fasta() containing seq_length column
#' @param family Character, specific family to plot. If NULL, shows all families.
#'   Special value "All" shows combined histogram. Default: NULL
#' @param bins Integer, number of histogram bins. Default: 30
#' @param title Character, plot title. Default: auto-generated
#' @param color Character or function, color for histogram. Default: "steelblue"
#' @param alpha Numeric, transparency for histogram bars. Default: 0.7
#' @param show_stats Logical, whether to show summary statistics. Default: TRUE
#' @param facet_families Logical, whether to create separate panels for each family. Default: FALSE
#' @param interactive Logical, whether to create interactive plot with plotly. Default: FALSE
#'
#' @return A ggplot object (or plotly object if interactive = TRUE)
#'
#' @examples
#' \dontrun{
#' kt1_data <- load_kt_fasta("kt1.fa", "kt1")
#'
#' # All families combined
#' plot_family_length(kt1_data, family = "All")
#'
#' # Specific family
#' plot_family_length(kt1_data, family = "M1")
#'
#' # Separate panels for each family
#' plot_family_length(kt1_data, facet_families = TRUE)
#'
#' # Interactive version
#' plot_family_length(kt1_data, interactive = TRUE)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_histogram facet_wrap theme_minimal
#' @importFrom ggplot2 labs theme element_text geom_vline annotate
#' @importFrom dplyr filter summarise group_by
#' @importFrom stats median
#' @export
plot_family_length <- function(data,
                              family = NULL,
                              bins = 30,
                              title = NULL,
                              color = "steelblue",
                              alpha = 0.7,
                              show_stats = TRUE,
                              facet_families = FALSE,
                              interactive = FALSE) {

  if (!"seq_length" %in% names(data)) {
    stop("Data must contain 'seq_length' column. Use load_kt_fasta() to load data properly.")
  }

  # Filter data based on family selection
  plot_data <- data

  if (!is.null(family) && family != "All") {
    if (!family %in% data$kt_family) {
      stop(paste("Family", family, "not found in data. Available families:",
                 paste(unique(data$kt_family), collapse = ", ")))
    }
    plot_data <- filter(data, kt_family == family)
    plot_title <- title %||% paste("Sequence Length Distribution -", family, "Family")
  } else {
    plot_title <- title %||% "Sequence Length Distribution - All Families"
  }

  # Create base plot
  p <- ggplot(plot_data, aes(x = seq_length))

  # Add histogram
  if (facet_families && (is.null(family) || family == "All")) {
    # Separate panels for each family
    p <- p +
      geom_histogram(bins = bins, fill = color, alpha = alpha, color = "white") +
      facet_wrap(~kt_family, scales = "free_y")
    plot_title <- title %||% "Sequence Length Distribution by Family"
  } else {
    # Single histogram
    p <- p + geom_histogram(bins = bins, fill = color, alpha = alpha, color = "white")
  }

  # Add statistics if requested
  if (show_stats) {
    stats_data <- plot_data %>%
      summarise(
        mean_length = mean(seq_length, na.rm = TRUE),
        median_length = median(seq_length, na.rm = TRUE),
        .groups = "drop"
      )

    p <- p +
      geom_vline(xintercept = stats_data$mean_length,
                 color = "red", linetype = "dashed", size = 1) +
      geom_vline(xintercept = stats_data$median_length,
                 color = "blue", linetype = "dashed", size = 1) +
      annotate("text",
               x = Inf, y = Inf,
               label = paste0("Mean: ", round(stats_data$mean_length, 1), "\n",
                             "Median: ", round(stats_data$median_length, 1)),
               hjust = 1.1, vjust = 1.1,
               color = "black", size = 3.5,
               fontface = "bold")
  }

  # Apply theme and labels
  p <- p +
    labs(
      title = plot_title,
      x = "Sequence Length (amino acids)",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      strip.text = element_text(size = 10, face = "bold")
    )

  # Make interactive if requested
  if (interactive) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("plotly package not available. Returning static plot.")
      return(p)
    } else {
      return(plotly::ggplotly(p))
    }
  }

  return(p)
}

#' Create Interactive Family Length Explorer
#'
#' Creates a Shiny app for interactive exploration of sequence length
#' distributions across different KT families.
#'
#' @param data A data.frame from load_kt_fasta()
#' @param launch Logical, whether to launch the app immediately. Default: TRUE
#'
#' @return A Shiny app object (if launch = FALSE)
#'
#' @examples
#' \dontrun{
#' kt1_data <- load_kt_fasta("kt1.fa", "kt1")
#' explore_family_lengths(kt1_data)
#' }
#'
#' @export
explore_family_lengths <- function(data, launch = TRUE) {

  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("shiny package required for interactive explorer. Install with install.packages('shiny')")
  }

  # Get available families
  families <- c("All", sort(unique(data$kt_family)))

  # Define UI
  ui <- shiny::fluidPage(
    shiny::titlePanel("KT Protein Family Length Explorer"),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::selectInput("family", "Select Family:",
                          choices = families,
                          selected = "All"),
        shiny::sliderInput("bins", "Number of Bins:",
                          min = 10, max = 100, value = 30),
        shiny::checkboxInput("show_stats", "Show Statistics", TRUE),
        shiny::checkboxInput("facet_families", "Separate Panels per Family", FALSE),
        shiny::hr(),
        shiny::h4("Dataset Summary"),
        shiny::verbatimTextOutput("summary")
      ),

      shiny::mainPanel(
        shiny::plotOutput("length_plot", height = "600px")
      )
    )
  )

  # Define Server
  server <- function(input, output, session) {

    output$length_plot <- shiny::renderPlot({
      plot_family_length(
        data = data,
        family = input$family,
        bins = input$bins,
        show_stats = input$show_stats,
        facet_families = input$facet_families
      )
    })

    output$summary <- shiny::renderText({
      if (input$family == "All") {
        subset_data <- data
      } else {
        subset_data <- dplyr::filter(data, kt_family == input$family)
      }

      paste(
        paste("Sequences:", nrow(subset_data)),
        paste("Families:", length(unique(subset_data$kt_family))),
        paste("Length range:", min(subset_data$seq_length, na.rm = TRUE),
              "-", max(subset_data$seq_length, na.rm = TRUE), "aa"),
        sep = "\n"
      )
    })
  }

  # Create and optionally launch app
  app <- shiny::shinyApp(ui = ui, server = server)

  if (launch) {
    shiny::runApp(app)
  } else {
    return(app)
  }
}