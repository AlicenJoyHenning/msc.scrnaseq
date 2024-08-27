#' @name homologize_core
#'
#' @title Homologize Core Function
#'
#' @description
#'    Helper function responsible for converting genes in a 'Seurat' object
#'    to homologous genes from another organism.
#'
#' @param project_name String with project or sample name.
#' @param seurat_dir Path to the 'Seurat' object used as input
#' @param organism_in Organism of the input 'Seurat' object (Hsap, Mmul, Mmus)
#' @param organism_out Organism for which gene homologs must be found
#' @param check_umap Whether or not marker umaps shoudl be generated. Default is TRUE.
#' @param output_dir Directory where homologize_core output should be generated.
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
#'   seurat_obj <- homologize_core(
#'     project_name  = "example_project",
#'     seurat_dir    = "path/to/seurat_object.rds",
#'     organism_in   = "human",
#'     organism_out  = "mouse",
#'     check_umap    = TRUE,
#'     output_dir    = "path/to/output"
#'   )
#' }

# Defining variables that exist within the homologize_core function space
utils::globalVariables(c("seurat_clusters", "avg_pdc_percent", "rank"))

homologize_core <- function(

  project_name,          # String with the name of the project
  seurat_dir,            # Input path to Seurat object
  organism_in,           # The organism that your sample is originally (human, mouse, macaque)
  organism_out,          # The gene homologs you want your sample to have (human, mouse, macaque)
  check_umap = TRUE,     # Return the marker visualizations
  output_dir =  "./"     # Where to store the output

){
  # Create output directory structure ####
  sub_dirs <- c("sample_clustering", "R_objects")

  for (sub_dir in sub_dirs) {

    # Generate the names of the directories we expect
    full_path <- file.path(output_dir, sub_dir)

    # Create directories if they don't exist
    if (!dir.exists(full_path)) {
      dir.create(full_path, recursive = FALSE)
    }

  }


  cat("\nBeginning ", organism_in, "~", organism_out, "homology mapping for", project_name, " ...\n")

  # Reading in files ####
  cat("Reading in files ...\n")

  # Read in annotation file of organism we have (first, distinct) and the organism we want it to be
  mapping_main_dir <- c("/home/alicen/Projects/interspecies/Homologize/gene_annotations")
  mapping_path <- file.path(mapping_main_dir, paste0(organism_in, "_", organism_out, "_mapping.rds"))
  mapping <- readRDS(mapping_path)

  # Read in Seurat (hopefully from Preprocessing)
  seurat <- readRDS(seurat_dir)


  # Accounting for Preprocess error ####
  col <- c("complexity")
  if (col %in% colnames(seurat@meta.data)) {
    seurat[[col]] <- NULL }


  # Perform the mapping ####
  cat("Performing mapping ...\n")

  # Extract count matrix
  seurat_matrix <- GetAssayData(seurat, layer = 'counts')

  # From count matrix, keep genes that are present in homologous df
  seurat_matrix <- seurat_matrix[rownames(seurat_matrix) %in% mapping[[organism_in]], ] # matching gene in count matrix (mouse) to mouse genes in mapping file (aka if any genes don't map they will be excluded)

  # Find exact output organism genes that match the input organism genes remaining in the count matrix
  replacementGenes <- match(rownames(seurat_matrix), mapping[[organism_in]])
  convertedGenes <- mapping[replacementGenes, organism_out]

  # Replace the gene names with homologs
  rownames(seurat_matrix) <- convertedGenes[[organism_out]]

  # Ensure row names of mouse_matrix are unique
  rownames(seurat_matrix) <- make.unique(rownames(seurat_matrix))

  # Create Seurat object
  seurat <- suppressWarnings(CreateSeuratObject(counts = seurat_matrix))
  seurat$organism <- paste0(organism_in, "_", organism_out)

  # Save output
  saveRDS(seurat, file.path(output_dir, "R_objects", paste0(project_name, "_", organism_in, "_", organism_out, ".rds")))


  if (check_umap) {

    cat("Creating UMAP visualisations ...\n")


    # UMAP check to visualize markers #####
    # No longer altering the same seurat object, this is just for visualisation of each sample
    seurat_check <- NormalizeData(seurat, verbose = FALSE) %>%
      FindVariableFeatures(verbose = FALSE) %>%
      ScaleData(verbose = FALSE) %>%
      RunPCA(verbose = FALSE) %>%
      FindNeighbors(dims = 1:10, verbose = FALSE) %>%
      FindClusters(verbose = FALSE) %>%
      RunUMAP(dims = 1:10, verbose = FALSE)


    # Define marker genes for organism out genes

    # Use the specified organism (Mmul, Hsap, Mmus) for correct gene nomenclature

    if (organism_out == "human") {

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

    if (organism_out == "macaque") {

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

    if (organism_out == "mouse") {

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

    markers_plot <- suppressWarnings(FeaturePlot(
      seurat_check,
      features = SpeciesFeatures,
      cols = c("#F1F1F1", "#9D71B3")) +
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
    colours_full <- c("#9FCBE1","#70B1D2", "#A1A9F5","#6C79F0", "#C5E7A7","#ABDC7E", "#ADDBB6", "#72BC7B", "#C9AFDF","#9D71B3",
                      "#9FCBE1","#70B1D2", "#A1A9F5","#6C79F0", "#C5E7A7","#ABDC7E", "#ADDBB6", "#72BC7B", "#C9AFDF","#9D71B3",
                      "#9FCBE1","#70B1D2", "#A1A9F5","#6C79F0", "#C5E7A7","#ABDC7E", "#ADDBB6", "#72BC7B", "#C9AFDF","#9D71B3")

    # Open a connection to a temporary file for writing (output verbose )
    tmp_conn <- file(tempfile(), open = "wt")
    sink(tmp_conn) # Redirect standard output and messages to the temporary file
    sink(tmp_conn, type = "message")

    # For label consistency with PreProcess plots
    cells <- length(Cells(seurat))

    plot_clusters <- DimPlot(seurat_check,
                             pt.size = 1,
                             label = TRUE,
                             label.box = TRUE,
                             label.size = 5,
                             label.color = "white",
                             repel = TRUE) +
      scale_color_manual(values = colours_full) +
      scale_fill_manual(values = colours_full) +
      labs(title = paste0(project_name, " sample clustering and marker visualisation"), subtitle = paste0(organism_in, " ~ ", organism_out, " gene homologs expressed in ", cells, " cells")) +
      xlab("UMAP 1") + ylab("UMAP 2") +
      theme(plot.title = element_text(hjust = 0.5, size = rel(1), face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, vjust = 1),
            panel.border = element_rect(colour = "black", fill=NA, linewidth =1))

    sink(NULL)
    sink(NULL, type = "message")
    close(tmp_conn) # Close connection

    UMAP_markers <- plot_clusters + wrap_plots(markers_plot) + plot_layout(ncol = 2)

    ggsave(filename = file.path(output_dir, "sample_clustering", paste0(project_name, "_", organism_in, "_", organism_out, ".png")), plot = UMAP_markers, width = 22, height = 9, dpi = 300)

  }

  cat("Complete!\n")

  return(seurat)

}
