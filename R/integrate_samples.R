#' @name integrate_samples
#'
#' @title Integrate Samples of Same Species Origin
#'
#' @param project_name String with project name.
#' @param organism Which organism do all inputs belong to, either Hsap, Mmul, or Mmus.
#' @param sample_list Required list of directories for input 'Seurat' objects.
#' @param output_dir Directory where integrate_samples output should be generated.
#'
#' @return Single, integrated Seurat object
#'
#' @import Seurat
#' @import ggplot2
#' @import patchwork
#' @import dplyr
#' @import stringr
#' @import purrr
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#'  # Define input list, here all are defined in terms of mouse genes
#'   samples <- list(
#'     "/home/user/postalignemnt/Homologize/R_objects/human_1_mouse.rds",
#'     "/home/user/postalignemnt/Homologize/R_objects/human_2_mouse.rds",
#'     "/home/user/postalignemnt/Homologize/R_objects/human_3_mouse.rds",
#'     "/home/user/postalignemnt/Homologize/R_objects/human_4_mouse.rds"
#'     )
#'
#'   # Run the integrate samples function
#'   output <- integrate_samples(
#'     project_name = "mouse_samples",
#'     organism = "mouse",
#'     sample_list = samples,
#'     output_dir = "/home/user/postalignment/integrate_samples/"
#'   )
#'
#' }

# Defining variables that exist within the integrate_samples function space
utils::globalVariables(c(
  "annoying", "cell_ids", "colours_full", "DC_Mon_genes", "FeaturesLabels", "mac_genes",
  "markers_plot", "merge_input", "merged_seurat", "method", "Mon_2_genes", "Mon_genes", "neu_genes",
  "organism", "output_dir", "plot_clusters", "project_name", "rank", "sample_list",
  "samples", "seurat", "seurat_check", "seurat_clusters", "seurat_merge", "SpeciesFeatures",
  "UMAP_markers"
))

integrate_samples <- function(

  project_name,
  organism, # human mouse macaque
  sample_list,
  output_dir

){
  # READ IN DATA ####
  cat("Reading in files...\n")

  # Initialize an empty list to store the Seurat objects
  samples <- list()

  # Read in the Seurat objects and store in list
  for (i in seq_along(sample_list)){
    samples[[paste0("sample_", i)]] <- readRDS(sample_list[[i]])
  }


  # MERGE DATA ####

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


  # PREPARE FOR INTEGRATION ####

  cat("Preparing merged object for integration ...\n")


  # Remove annoying columns that somehow managed to stay alive

  annoying <- c("ptprc.percent", "quality", "IMC", "complexity")

  for (annoy in annoying) {
    if (annoy %in% colnames(seurat_merge@meta.data)) {
      seurat_merge[[annoy]] <- NULL
    }
  }

  seurat_merge <- NormalizeData(seurat_merge) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA()


  # INTEGRATE ####

  cat("Integrating ...\n")

  # Integrate the layers of the merged object

  seurat <- IntegrateLayers(
    object = seurat_merge,
    method = CCAIntegration,
    layers = c("counts", "data", "scale.data"),
    scale.layer = "scale.data",
    assay = "RNA",
    orig.reduction = "pca",
    new.reduction = "integrated")

  saveRDS(seurat, file.path(output_dir, paste0(project_name, "_temp.rds"))) # incase error later on to save redoing this


  # DIMENSIONALITY REDUCTION ####
  # Run UMAP reduction

  cat("Reducing dimensions ...\n")

  seurat  <- FindNeighbors(merged_seurat,
                           reduction = "integrated",
                           dims = 1:30) %>%

    FindClusters(resolution = 0.6,
                 group.singletons = TRUE,
                 cluster.name = "cca_clusters") %>%

    RunUMAP(reduction = "integrated",
            dims = 1:30)


  saveRDS(seurat, file.path(output_dir, paste0(project_name, ".rds")))
  file.remove(file.path(output_dir, paste0(project_name, "_temp.rds")))


  # Visualising ####

  if (organism == "human") {

    # Highlight neutrophils but exlude NK cells that may express CXCR2
    neu_genes <- c("CXCR2")
    seurat_check[['neu.percent']] <- PercentageFeatureSet(
      object   = seurat_check,
      features = intersect(neu_genes, rownames(seurat_check@assays$RNA)),
      assay    = "RNA"
    )

    nk_genes <- c("NKG7")
    seurat_check[['nk.percent']] <- PercentageFeatureSet(
      object   = seurat_check,
      features = intersect(nk_genes, rownames(seurat_check@assays$RNA)),
      assay    = "RNA"
    )
    seurat_check[['neu.percent']] <- ifelse(seurat_check$nk.percent < seurat_check$neu.percent, seurat_check$neu.percent, 0)


    # macrophages, express MARCO but some monocytes do aswell
    mac_genes <- c("MARCO")
    seurat_check[['mac.percent']] <- PercentageFeatureSet(
      object = seurat_check,
      features = intersect(mac_genes, rownames(seurat_check@assays$RNA)),
      assay = "RNA"
    )

    mon_genes <- c("S100A8")
    seurat_check[['mon.percent']] <- PercentageFeatureSet(
      object = seurat_check,
      features = intersect(mon_genes, rownames(seurat_check@assays$RNA)),
      assay = "RNA"
    )
    seurat_check[['mac.percent']] <- ifelse(seurat_check$mac.percent > seurat_check$mon.percent, seurat_check$mac.percent, 0)

    # Isolating DC by those that express CPVL more than they express CD14
    DC_Mon_genes <- c("CPVL")
    seurat_check[['DC.percent']] <- PercentageFeatureSet(
      object = seurat_check,
      features = intersect(DC_Mon_genes, rownames(seurat_check@assays$RNA)),
      assay = "RNA"
    )

    Mon_genes <- c("CD14")
    seurat_check[['Mon.percent']] <- PercentageFeatureSet(
      object = seurat_check,
      features = intersect(Mon_genes, rownames(seurat_check@assays$RNA)),
      assay = "RNA"
    )
    seurat_check[['DC.percent']] <- ifelse((seurat_check$DC.percent > seurat_check$Mon.percent) & (seurat_check$DC.percent > seurat_check$mac.percent), seurat_check$DC.percent, 0)
    # ensure both monocytes & macrophages with CPVL expression are silenced

    SpeciesFeatures <- c("PTPRC", "MS4A1", "CD3E", "FOXP3", "PRF1", "S100A8", "neu.percent", "mac.percent", "MKI67", "CPA3", "DC.percent", "CLEC4C")
    FeaturesLabels <- c("Immune cells", "B cells", "T cells", "T regulatory", "Natural Killer", "Monocytes", "Neutrophils", "Macrophages", "Proliferative", "Mast cells", "Dendritic cells", "P. Dendritic cells")

  }

  if (organism == "macaque") {

    nk_genes <- c("NKG7")
    t_genes <- c("CD3E")
    seurat_check[['nk.percent']] <- PercentageFeatureSet(
      object   = seurat_check,
      features = intersect(nk_genes, rownames(seurat_check@assays$RNA)),
      assay    = "RNA"
    )
    seurat_check[['t.percent']] <- PercentageFeatureSet(
      object   = seurat_check,
      features = intersect(t_genes, rownames(seurat_check@assays$RNA)),
      assay    = "RNA"
    )
    seurat_check[['nk.percent']] <- ifelse(seurat_check$nk.percent < seurat_check$t.percent, 0, seurat_check$nk.percent)


    # Both mono & DCs express CPVL but DC do so at a higer level, isolating DCs by isolating cells where CPVL expression is higher than CD14 expression
    DC_Mon_genes <- c("CPVL")
    seurat_check[['DC_Mon.percent']] <- PercentageFeatureSet(
      object   = seurat_check,
      features = intersect(DC_Mon_genes, rownames(seurat_check@assays$RNA)),
      assay    = "RNA"
    )

    Mon_genes <- c("CD14")
    Mon_2_genes <- c("CD163")
    seurat_check[['Mon.percent']] <- PercentageFeatureSet(
      object   = seurat_check,
      features = intersect(Mon_genes, rownames(seurat_check@assays$RNA)),
      assay    = "RNA"
    ) +
      PercentageFeatureSet(
        object   = seurat_check,
        features = intersect(Mon_2_genes, rownames(seurat_check@assays$RNA)),
        assay    = "RNA"
      )

    seurat_check[['DC.percent']] <- ifelse(seurat_check$DC_Mon.percent >= seurat_check$Mon.percent,  (seurat_check$DC_Mon.percent - seurat_check$Mon.percent), seurat_check$DC_Mon.percent)

    # Macrophages : all types
    Marco_genes <- c("MARCO")
    seurat_check[['marco.percent']] <- PercentageFeatureSet(
      object   = seurat_check,
      features = intersect(Marco_genes, rownames(seurat_check@assays$RNA)),
      assay    = "RNA"
    )
    Mrc1_genes <- c("MRC1")
    seurat_check[['mrc1.percent']] <- PercentageFeatureSet(
      object   = seurat_check,
      features = intersect(Mrc1_genes, rownames(seurat_check@assays$RNA)),
      assay    = "RNA"
    )
    seurat_check[['mac.percent']] <- seurat_check$marco.percent + seurat_check$mrc1.percent

    # Cheating a bit to visualise pDCs :
    pdc_genes <- c("CDH15") # alt: CDK5R1,
    seurat_check[['pdc.percent']] <- PercentageFeatureSet(
      object   = seurat_check,
      features = intersect(pdc_genes, rownames(seurat_check@assays$RNA)),
      assay    = "RNA"
    )

    # Calculate the average CDH15 for each cluster, add ranks, and identify the best cluster
    pDC_cluster <- seurat_check@meta.data %>%
      group_by(seurat_clusters) %>%
      summarise(
        avg_pdc_percent = mean(pdc.percent, na.rm = TRUE)
      ) %>%
      mutate(
        rank = rank(-avg_pdc_percent)  # Negative for descending order (highest mt to be #1)
      ) %>%
      arrange(rank) %>%
      slice(1) %>%
      pull(seurat_clusters) %>%
      as.character()

    # Highlight the cluster where CDH15 is highest
    seurat_check[['pdc.percent']] <- ifelse(seurat_check$seurat_clusters == pDC_cluster, 1, 0)

    SpeciesFeatures <- c("PTPRC", "MS4A1", "CD3E", "FOXP3", "nk.percent", "Mon.percent", "IL1R2", "mac.percent", "MKI67", "CPA3", "DC.percent", "pdc.percent")
    FeaturesLabels <- c("Immune cells", "B cells", "T cells", "T regulatory", "Natural Killer", "Monocytes", "Neutrophils", "Macrophages", "Proliferative", "Mast cells", "Dendritic cells", "P. Dendritic cells")

  }

  if (organism == "mouse") {

    # Define marker gene sets for cell types difficult to distinguish globally by one gene
    neu_genes <- c("S100a8", "Cxcr2", "Ly6g")
    seurat_check[['neu.percent']] <- PercentageFeatureSet(
      object   = seurat_check,
      features = intersect(neu_genes, rownames(seurat_check@assays$RNA)),
      assay    = "RNA"
    )

    mac_genes <- c("Mrc1", "Marco")
    seurat_check[['mac.percent']] <- PercentageFeatureSet(
      object   = seurat_check,
      features = intersect(mac_genes, rownames(seurat_check@assays$RNA)),
      assay    = "RNA"
    )

    mon_genes <- c("Ccr2", "Itgam")
    seurat_check[['mon.percent']] <- PercentageFeatureSet(
      object   = seurat_check,
      features = intersect(mon_genes, rownames(seurat_check@assays$RNA)),
      assay    = "RNA"
    )

    den_genes <- c("Clec9a", "Flt3")
    seurat_check[['den.percent']] <- PercentageFeatureSet(
      object   = seurat_check,
      features = intersect(den_genes, rownames(seurat_check@assays$RNA)),
      assay    = "RNA"
    )

    nk_genes <- c("Klrb1c", "Prf1")
    seurat_check[['nk.percent']] <- PercentageFeatureSet(
      object   = seurat_check,
      features = intersect(nk_genes, rownames(seurat_check@assays$RNA)),
      assay    = "RNA"
    )

    mast_genes <- c("Fcer1a", "Cd200r3")
    seurat_check[['mast.percent']] <- PercentageFeatureSet(
      object   = seurat_check,
      features = intersect(mast_genes, rownames(seurat_check@assays$RNA)),
      assay    = "RNA"
    )

    SpeciesFeatures <- c("Ptprc", "Ms4a1", "Cd3e", "Foxp3", "nk.percent", "mon.percent", "neu.percent", "mac.percent", "Mki67", "mast.percent", "den.percent", "Siglech")
    FeaturesLabels  <- c("Immune cells", "B cells", "T cells", "T regulatory", "Natural Killer", "Monocytes", "Neutrophils", "Macrophages", "Proliferative", "Mast cells", "Dendritic cells", "P. Dendritic cells")
  }


  # UMAP plots shaded by gene expression

  DefaultAssay(seurat) <- "RNA"

  markers_plot <- suppressWarnings(FeaturePlot(
    seurat,
    features = SpeciesFeatures,
    cols = c("#F1F1F1", "#46948D")) +
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
  # Colours drawn from RColorBrewer::brewer.pal(12, 'Paired') with change to yellow
  colours_full <- c( "#70B1D2", "#185D8C", "#83CB41", "#277921", "#9D71B3", "#482A68",
                     "#70B1D2", "#185D8C", "#83CB41", "#277921", "#9D71B3", "#482A68",
                     "#70B1D2", "#185D8C", "#83CB41", "#277921", "#9D71B3", "#482A68",
                     "#70B1D2", "#185D8C", "#83CB41", "#277921", "#9D71B3", "#482A68",
                     "#70B1D2", "#185D8C", "#83CB41", "#277921", "#9D71B3", "#482A68",
                     "#70B1D2", "#185D8C", "#83CB41", "#277921", "#9D71B3", "#482A68")

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
    labs(title = paste0(project_name, " integrated sample clustering and marker visualisation"), subtitle = paste0("Total cells:", cells)) +
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
