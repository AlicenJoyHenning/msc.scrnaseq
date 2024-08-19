#' @name preprocess
#'
#' @title Preprocess Full Function
#'
#' @description
#'    Function responsible for heavy lifting preprocessing work.
#'    This includes ambient RNA correction with SoupX, removing red blood cells,
#'    isolating immune cells, filtering damaged cells detected by DropletQC and
#'    limiric, and filtering doublets detected by DoubletFinder. Also provided
#'    is preliminary check of sample quality by immune cell marker visualisation (UMAP).
#'
#'
#' @param project_name String with project or sample name.
#' @param organism String with one of the following Mmul, Hsap, Mmus.
#' @param ptprc_threshold Percentage CD45 over which cells will be retained. E.g: If 0, all cells must express CD45
#' @param filtered_path Directory of filtered alignment output, matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz.
#' @param raw_path Directory of raw alignment output, matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz.
#' @param velocyto_path Directory of Velocyto filtered alignment output, spliced.mtx.gz, unspliced.mtx.gz, ambiguous.mtx.gz, features.tsv.gz, barcodes.tsv.gz.
#' @param filter_damaged Should output contain no damaged cells. Default is TRUE.
#' @param output_dir Directory where preprocess_core output should be generated.
#' @param sample_list Input multiple samples in list. Default is NULL
#'
#' @return Seurat object
#'
#' @import Seurat
#' @import SoupX
#' @import Matrix
#' @import ggplot2
#' @import dplyr
#' @import DoubletFinder
#' @importFrom cowplot draw_label ggdraw
#' @importFrom patchwork align_plots
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Function can run multiple samples or a single sample
#'
#'   # Run a single sample
#'   SRR1234567 <- preprocess_core(project_name = "SRR1234567",
#'                                 filtered_path = "/home/username/alignment/SRR1234567/filtered",
#'                                 raw_path = "/home/username/alignment/SRR1234567/raw",
#'                                 velocyto_path = "/home/username/alignment/SRR1234567/velocyto",
#'                                 output_dir = "/home/username/postalignemnt/SRR1234567/preprocess")
#'
#'  # Run multiple samples
#'  samples <- list(
#'  list(project_name = "SRR1234567",
#'      filtered_path = "/home/username/alignment/SRR1234567/filtered",
#'      raw_path = "/home/username/alignment/SRR1234567/raw",
#'      velocyto_path = "/home/username/alignment/SRR1234567/velocyto",
#'      output_dir = "/home/username/postalignemnt/SRR1234567/preprocess"),
#'
#'  list(project_name = "SRR1234568",
#'      filtered_path = "/home/username/alignment/SRR1234568/filtered",
#'      raw_path = "/home/username/alignment/SRR1234568/raw",
#'      velocyto_path = "/home/username/alignment/SRR1234568/velocyto",
#'      output_dir = "/home/username/postalignemnt/SRR1234568/preprocess"))
#'
#' GSE1234567 <- preprocess(sample_list = samples)
#'
#' }

# Defining variables that exist within the preprocess_core function space
utils::globalVariables(c(
  "current_plots", "end_index", "IMCQC", "IMCQC_path", "nrow", "num_current_plots",
  "num_plots", "page_num", "plot_grid", "plots", "plots_per_page", "png_path",
  "RBCQC", "RBCQC_path", "rds_dir", "rds_files"
))

preprocess <- function(
    organism        = NULL,
    ptprc_threshold = NULL,
    filtered_path   = NULL,
    raw_path        = NULL,
    velocyto_path   = NULL,
    project_name    = NULL,
    filter_damaged  = NULL,
    output_dir      = NULL,
    sample_list     = NULL
){

  # First, handle two possible cases, either one sample as input or multiple in the form of a list

  if (is.null(sample_list)) {

    # Run single sample (default) using the core preprocess function
    preprocess_core(
      organism        = organism,
      ptprc_threshold = ptprc_threshold,
      filtered_path   = filtered_path,
      raw_path        = raw_path,
      velocyto_path   = velocyto_path,
      project_name    = project_name,
      filter_damaged  = filter_damaged,
      output_dir      = output_dir
    )

    # Clean output directories: convert QC rds to png plots

    # Red blood cell QC
    RBCQC_path <- file.path(output_dir, "RBCQC", paste0(project_name, "_RBCQC", ".rds"))
    RBCQC      <- readRDS(RBCQC_path)
    ggsave(file.path(output_dir, "RBCQC", paste0(project_name, "_RBCQC", ".png")), plot = RBCQC, width = 5, height = 5, dpi = 300)
    file.remove(RBCQC_path)

    # And immune QC
    IMCQC_path <- file.path(output_dir, "IMCQC", paste0(project_name, "_IMCQC", ".rds"))
    IMCQC      <- readRDS(IMCQC_path)
    ggsave(file.path(output_dir, "IMCQC", paste0(project_name, "_IMCQC", ".png")), plot = IMCQC, width = 5, height = 5, dpi = 300)
    file.remove(IMCQC_path)

  }

  # If there is a sample list, extract elements from each list item and run with core preprocess function
  else

  {

    # Initialize the output list where resulting Seurat objects will be stored
    results <- list()

    # Run preprocess function for each sample
    for (i in seq_along(sample_list)) {

      # Define sample
      sample <- sample_list[[i]]

      # Define the inputs from the list
      organism        <- sample$organism
      ptprc_threshold <- sample$ptprc_threshold
      filtered_path   <- sample$filtered_path
      raw_path        <- sample$raw_path
      velocyto_path   <- sample$velocyto_path
      project_name    <- sample$project_name
      filter_damaged  <- sample$filter_damaged
      output_dir      <- sample$output_dir

      # Call the core function with error handling
      tryCatch({

        temp_result <- preprocess_core(
          organism        = organism,
          ptprc_threshold = ptprc_threshold,
          filtered_path   = filtered_path,
          raw_path        = raw_path,
          velocyto_path   = velocyto_path,
          project_name    = project_name,
          filter_damaged  = filter_damaged,
          output_dir      = output_dir
        )

        results[[project_name]] <- temp_result$seuratObject

      },

      # If error in one sample, print message and continue with next sample
      error = function(e) {
        cat(paste(e$message, "\n", "Error in processing", project_name, "\n", "Likely low cell number\n\n"))

      })
    }

    # Collate output QC plots ####
    # RBC
    rds_dir   <- file.path(output_dir, "RBCQC") # Define the directory containing the RDS files
    rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)
    plots     <- lapply(rds_files, readRDS) # Read all RDS files into a list of ggplot objects

    # Function to save a grid of plots(5 in each row and max 3 rows)
    save_plot_grid <- function(plots, file_path, ncol = 5, nrow = 3) {
      plot_grid <- plot_grid(plotlist = plots, ncol = ncol, nrow = nrow)
      ggsave(file_path, plot = plot_grid, width = ncol * 5, height = nrow * 5, dpi = 300)
    }

    # Loop through the plots & calculate # plot pages nedded (max 15 plots each pg)
    num_plots <- length(plots)
    plots_per_page <- 15
    page_num <- 1

    for (i in seq(1, num_plots, by = plots_per_page)) {
      end_index <- min(i + plots_per_page - 1, num_plots)
      current_plots <- plots[i:end_index]

      # Calculate the number of rows needed for the current set of plots
      num_current_plots <- length(current_plots)
      nrow <- ceiling(num_current_plots / 5)

      # Save the current set of plots as a PNG
      png_path <- file.path(output_dir, "RBCQC", paste0(project_name, "_RBCQC_", page_num, ".png"))
      save_plot_grid(current_plots, png_path, ncol = 5, nrow = nrow)

      page_num <- page_num + 1
    }

    # Delete all rds files (only want output image)
    file.remove(rds_files)

    # Collate output IMC QC plots ####
    # IMC
    rds_dir   <- file.path(output_dir, "IMCQC") # Define the directory containing the RDS files
    rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)
    plots     <- lapply(rds_files, readRDS) # Read all RDS files into a list of ggplot objects

    # Function to save a grid of plots(5 in each row and max 3 rows)
    save_plot_grid <- function(plots, file_path, ncol = 5, nrow = 3) {
      plot_grid <- plot_grid(plotlist = plots, ncol = ncol, nrow = nrow)
      ggsave(file_path, plot = plot_grid, width = ncol * 5, height = nrow * 5, dpi = 300)
    }

    # Loop through the plots & calculate # plot pages
    num_plots <- length(plots)
    plots_per_page <- 15
    page_num <- 1

    for (i in seq(1, num_plots, by = plots_per_page)) {
      end_index <- min(i + plots_per_page - 1, num_plots)
      current_plots <- plots[i:end_index]

      # Calculate the number of rows needed for the current set of plots
      num_current_plots <- length(current_plots)
      nrow <- ceiling(num_current_plots / 5)

      # Save the current set of plots as a PNG
      png_path <- file.path(output_dir, "IMCQC", paste0(project_name, "_IMCQC_", page_num, ".png"))
      save_plot_grid(current_plots, png_path, ncol = 5, nrow = nrow)

      page_num <- page_num + 1
    }

    # Delete all rds files (only want output image)
    file.remove(rds_files)


    return(results)


  }

}
