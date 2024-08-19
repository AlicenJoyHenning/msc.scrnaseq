#' @name homologize_core
#'
#' @title Homologize Core Function
#'
#' @description
#'    Helper function responsible for converting genes in a Seurat object
#'    to homologous genes from another organism.
#'
#' @param project_name String with project or sample name.
#' @param seurat_dir Path to the Seurat object used as input.
#' @param organism_in Organism of the input Seurat object (Hsap, Mmul, Mmus)
#' @param organism_out Organism for which gene homologs must be found
#' @param check_umap Whether or not marker umaps shoudl be generated. Default is TRUE.
#' @param output_dir Directory where homologize_core output should be generated.
#' @param sample_list Input multiple samples in list. Default is NULL
#'
#' @return Seurat object
#'
#' @import Seurat
#' @importFrom dplyr group_by summarise mutate arrange slice pull
#' @importFrom ggplot2 theme element_text element_blank element_rect labs scale_color_manual scale_fill_manual
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom utils globalVariables
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   seurat_obj <- homologize(
#'     project_name = "example_project",
#'     seurat_dir   = "path/to/seurat_object.rds",
#'     organism_in  = "human",
#'     organism_out = "mouse",
#'     check_umap   = TRUE,
#'     output_dir   = "path/to/output"
#'   )
#' }

homologize <- function(

  project_name    = NULL,
  seurat_dir      = NULL,
  organism_in     = NULL,
  organism_out    = NULL,
  check_umap      = NULL,
  output_dir      = NULL,
  sample_list     = NULL

){

  # Two cases, either one sample as input or multiple in the form of a list

  if (is.null(sample_list)) {

    # Run single sample (default) using the core Homologize function
    temp <- homologize_core(

      project_name = project_name,
      seurat_dir   = seurat_dir,
      organism_in  = organism_in,
      organism_out = organism_out,
      check_umap   = check_umap,
      output_dir   =  output_dir

    )

    return(temp)

  }

  else

  {

    # Initialize the output list where resulting Seurat objects will be stored
    results <- list()

    # Run Homologize function for each sample
    for (i in seq_along(sample_list)) {

      # Define sample
      sample <- sample_list[[i]]

      # Define the inputs from the list
      project_name <- sample$project_name
      seurat_dir   <- sample$seurat_dir
      organism_in  <- sample$organism_in
      organism_out <- sample$organism_out
      check_umap   <- sample$check_umap
      output_dir   <- sample$output_dir

      # Account for defaults
      if (is.null(check_umap)) { check_umap = TRUE}
      # if (is.null(output_dir)) { output_dir = "./"}

      # Call the core function with error handling
      tryCatch({

        temp <- homologize_core(

          project_name = project_name,
          seurat_dir   = seurat_dir,
          organism_in  = organism_in,
          organism_out = organism_out,
          check_umap   = check_umap,
          output_dir   =  output_dir

        )

        results[[project_name]] <- temp

      },

      # If error in one sample, print message and continue with next sample
      error = function(e) {
        cat(paste(e$message, "\n", "Error in processing", project_name, "\n", "Likely low cell number\n\n"))

      })

    }

    return(results)
  }

}
