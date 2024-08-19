#' @name integrate_species
#'
#' @title Integrate Samples of Different Species Origin
#'
#' @param project_name String with project name.
#' @param organism Organism that all samples gene annotations conform to.
#' @param sample_list Required list of directories for input Seurat objects.
#' @param method Integration method to use, either CCA or harmony, default is harmony.
#' @param output_dir Directory where integrate_species output should be generated.
#'
#' @return Single, integrated Seurat object
#'
#' @import dplyr
#' @import Seurat
#' @import harmony
#' @import ggplot2
#' @import patchwork
#' @import stringr
#' @import purrr
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  # Define input list, here all are defined in terms of mouse genes
#'   samples <- list(
#'     "/home/user/postalignemnt/Homologize/R_objects/human_1_mouse.rds",
#'     "/home/user/postalignemnt/Homologize/R_objects/human_2_mouse.rds",
#'     "/home/user/postalignemnt/Homologize/R_objects/mouse_1_mouse.rds",
#'     "/home/user/postalignemnt/Homologize/R_objects/mouse_2_mouse.rds",
#'     "/home/user/postalignemnt/Homologize/R_objects/macaque_1_mouse.rds",
#'     "/home/user/postalignemnt/Homologize/R_objects/macaque_2_mouse.rds")
#'
#'   # Run the integrate species function
#'   output <- integrate_species(
#'     project_name = "mouse_genes_harmony",
#'     organism = "mouse",
#'     sample_list = samples,
#'     method = "harmony",
#'     output_dir = "/home/user/postalignment/integrate_species/"
#'   )
#'
#' }

# Defining variables that exist within the integrate_species function space
utils::globalVariables(c("cell_ids", "colours_full", "colours_split",
                         "FeaturesLabels", "markers_plot", "merge_input",
                         "method", "organism", "output_dir", "plot_clusters",
                         "rank", "sample_list", "samples", "seurat",
                         "seurat_clusters", "seurat_merge", "SpeciesFeatures",
                         "species_split", "UMAP_markers"))

integrate_species <- function(project_name,
                              organism,
                              sample_list,
                              method,
                              output_dir
){

  # Read in data --------------
  cat("Reading in files...\n")

  # Initialize an empty list to store the Seurat objects
  samples <- list()

  # Read in the Seurat objects and store in list
  for (i in seq_along(sample_list)){
    samples[[paste0("sample_", i)]] <- readRDS(sample_list[[i]])
  }


  # Merge data --------------

  cat("Merging the Seurat objects ...\n")

  # Prepare arguments for merge

  # Create merge input using input number of list entries
  merge_input <- list()
  cell_ids <- c()

  for (i in seq_along(samples)) {
    merge_input <- c(merge_input, samples[[i]])
    cell_ids <- c(cell_ids, paste0("sample_", i))
  }

  # Remove the first item from merge_input as it will be specified separately in the merge function
  merge_input <- merge_input[-1]

  # Dynamically pass the arguments to the merge function
  seurat_merge <- do.call(merge, c(list(x = samples[[1]], y = merge_input, add.cell.ids = cell_ids)))


  # Prepare for integration --------------

  cat("Preparing merged object for integration ...\n")


  seurat_merge <- NormalizeData(seurat_merge) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA()


  # Integrate --------------

  cat("Integrating ...\n")

  # Integrate the layers of the merged object

  if (method == "harmony"){

    seurat <- IntegrateLayers(
      object = seurat_merge,
      method = HarmonyIntegration,
      theta = 3,
      layers = c("data"), # data
      assay = "RNA",
      orig.reduction = "pca",
      new.reduction = "integrated")

  }

  if (method == "cca"){

    seurat <- IntegrateLayers(
      object = seurat_merge,
      method = CCAIntegration,
      layers = c("data"), # data
      assay = "RNA",
      orig.reduction = "pca",
      new.reduction = "integrated")

  }


  saveRDS(seurat, file.path(output_dir, paste0(project_name, "_temp.rds"))) # incase error later on to save redoing this


  # Dimensionality reduction --------------
  # Run UMAP reduction

  cat("Reducing dimensions ...\n")

  seurat  <- FindNeighbors(seurat, reduction = "integrated", dims = 1:30) %>%

    FindClusters(resolution = 0.2) %>%

    RunUMAP(reduction = "integrated", dims = 1:30)


  cells <- length(Cells(seurat))

  saveRDS(seurat, file.path(output_dir, paste0(project_name, ".rds")))
  file.remove(file.path(output_dir, paste0(project_name, "_temp.rds")))


  # Visualising --------------

  cat("Visualising samples...\n")

  # See the split of different species

  colours_split <- c("#6C79F0", "#A97ECC", "#55AFA6")

  species_split <- DimPlot(seurat,
                           pt.size = 0.7,
                           group.by = "organism",
                           split.by = "organism"
                           #label = TRUE,
                           #label.box = TRUE,
                           #label.size = 5,
                           #label.color = "white",
                           #repel = TRUE
  ) +
    scale_color_manual(values = colours_split) +
    scale_fill_manual(values = colours_split) +
    labs(title = paste0(project_name, " integrated sample clustering"), subtitle = paste0("Total cells: ", cells)) +
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme(plot.title = element_text(hjust = 0.5, size = rel(1), face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, vjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, linewidth =1))

  ggsave(filename = file.path(output_dir, paste0(project_name, "_split.png")), plot = species_split, width = 18, height = 6, dpi = 300)


  if (organism == "human") {

    SpeciesFeatures <- c("PTPRC", "MS4A1", "CD3E", "FOXP3", "PRF1", "S100A8", "CXCR2", "MARCO", "MKI67", "CPA3", "CPVL", "CLEC4C")
    FeaturesLabels <- c("Immune cells", "B cells", "T cells", "T regulatory", "Natural Killer", "Monocytes", "Neutrophils", "Macrophages", "Proliferative", "Mast cells", "Dendritic cells", "P. Dendritic cells")

  }

  if (organism == "macaque") {

    SpeciesFeatures <- c("PTPRC", "MS4A1", "CD3E", "FOXP3", "NKG7", "CD163", "IL1R2", "MRC1", "MKI67", "CPA3", "CPVL", "CDH15")
    FeaturesLabels <- c("Immune cells", "B cells", "T cells", "T regulatory", "Natural Killer", "Monocytes", "Neutrophils", "Macrophages", "Proliferative", "Mast cells", "Dendritic cells", "P. Dendritic cells")

  }

  if (organism == "mouse") {

    # Define marker gene sets for cell types difficult to distinguish globally by one genes
    SpeciesFeatures <- c("Ptprc", "Ms4a1", "Cd3e", "Foxp3", "Nkg7", "Ccr2", "S100a8", "Marco", "Mki67", "Fcer1a", "Clec9a", "Siglech")
    FeaturesLabels  <- c("Immune cells", "B cells", "T cells", "T regulatory", "Natural Killer", "Monocytes", "Neutrophils", "Macrophages", "Proliferative", "Mast cells", "Dendritic cells", "P. Dendritic cells")

  }


  # UMAP plots shaded by gene expression

  DefaultAssay(seurat) <- "RNA"

  markers_plot <- suppressWarnings(FeaturePlot(
    seurat,
    features = SpeciesFeatures,
    cols = c("#F1F1F1", "#3E77B6")) +
      theme(
        axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(),
        plot.caption = element_blank()
      ))

  for (plot in 1:length(markers_plot)) {
    markers_plot[[plot]] <- markers_plot[[plot]] +
      NoLegend() +
      NoAxes() +
      labs(title = FeaturesLabels[plot]) +
      theme(panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
            plot.title = element_text(size = 12))
  }

  # UMAP of clusters
  colours_full <- c("#9FCBE1","#70B1D2", "#A1A9F5","#6C79F0", "#C5E7A7","#ABDC7E", "#ADDBB6", "#72BC7B", "#C9AFDF","#9D71B3",
                    "#9FCBE1","#70B1D2", "#A1A9F5","#6C79F0", "#C5E7A7","#ABDC7E", "#ADDBB6", "#72BC7B", "#C9AFDF","#9D71B3",
                    "#9FCBE1","#70B1D2", "#A1A9F5","#6C79F0", "#C5E7A7","#ABDC7E", "#ADDBB6", "#72BC7B", "#C9AFDF","#9D71B3",
                    "#9FCBE1","#70B1D2", "#A1A9F5","#6C79F0", "#C5E7A7","#ABDC7E", "#ADDBB6", "#72BC7B", "#C9AFDF","#9D71B3",
                    "#9FCBE1","#70B1D2", "#A1A9F5","#6C79F0", "#C5E7A7","#ABDC7E", "#ADDBB6", "#72BC7B", "#C9AFDF","#9D71B3",
                    "#9FCBE1","#70B1D2", "#A1A9F5","#6C79F0", "#C5E7A7","#ABDC7E", "#ADDBB6", "#72BC7B", "#C9AFDF","#9D71B3",
                    "#9FCBE1","#70B1D2", "#A1A9F5","#6C79F0", "#C5E7A7","#ABDC7E", "#ADDBB6", "#72BC7B", "#C9AFDF","#9D71B3")

  # Open a connection to a temporary file for writing (output verbose )
  tmp_conn <- file(tempfile(), open = "wt")
  sink(tmp_conn) # Redirect standard output and messages to the temporary file
  sink(tmp_conn, type = "message")

  # For label consistency with PreProcess plots
  cells <- length(Cells(seurat))

  plot_clusters <- DimPlot(seurat,
                           pt.size = 1,
                           group.by = "orig.ident",
                           label = TRUE,
                           label.box = TRUE,
                           label.size = 5,
                           label.color = "white",
                           repel = TRUE) +
    scale_color_manual(values = colours_full) +
    scale_fill_manual(values = colours_full) +
    labs(title = paste0(project_name, " integrated sample clustering and marker visualisation"), subtitle = paste0("Total cells: ", cells)) +
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme(plot.title = element_text(hjust = 0.5, size = rel(1), face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, vjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, linewidth =1))

  sink(NULL)
  sink(NULL, type = "message")
  close(tmp_conn) # Close connection

  UMAP_markers <- plot_clusters + wrap_plots(markers_plot) + plot_layout(ncol = 2)

  ggsave(filename = file.path(output_dir, paste0(project_name, ".png")), plot = UMAP_markers, width = 20, height = 9, dpi = 300)

  cat("Complete!")

}
