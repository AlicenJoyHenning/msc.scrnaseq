#' @name preprocess_core
#'
#' @title Preprocess Core Function
#'
#' @description
#'    Helper function responsible for heavy lifting preprocessing work.
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
#'
#' @return Seurat object
#'
#' @import Seurat
#' @import SoupX
#' @import Matrix
#' @import ggplot2
#' @import DropletQC
#' @import dplyr
#' @import DoubletFinder
#' @importFrom stats setNames
#' @importFrom cowplot draw_label ggdraw
#' @importFrom patchwork align_plots
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Note that this function runs a single sample
#'   SRR1234567 <- preprocess_core(project_name = "SRR1234567",
#'                                 filtered_path = "/home/username/alignment/SRR1234567/filtered",
#'                                 raw_path = "/home/username/alignment/SRR1234567/raw",
#'                                 velocyto_path = "/home/username/alignment/SRR1234567/velocyto",
#'                                 output_dir = "/home/username/postalignemnt/SRR1234567/preprocess")
#'
#' }

# Defining variables that exist within the preprocess_core function space
utils::globalVariables(c(
  "annotations", "avg_mt_percent", "avg_pdc_percent", "avg_rb_percent", "BCmetric",
  "best_cluster", "Cells_before_doublets", "check", "clusters", "colours_full",
  "columns", "combined_rank", "complexity_metric", "complexity_plot", "Damaged_cells",
  "Damaged_percent", "df_col_name", "DF", "doublet_classifications", "Doublet_cells",
  "end_time", "exon_sum", "full_path", "gene_name", "hbb.percent", "hemo_feature_QC",
  "IMC_percent", "imcDf", "immune_feature_QC", "initial_cells", "matrix", "markers_plot",
  "mt.percent", "mt_plot", "mt_rb_genes", "nCount_RNA", "nFeature_RNA", "new_row", "nf",
  "pdc.percent", "plot_clusters", "plot_layout", "plots", "ptprc.percent", "preprocess_output",
  "QC", "quality", "rank_mt", "rank_rb", "rb.percent", "rb_plot", "rbcDf", "RBC_number",
  "RBC_percent", "Remaining_cells", "rho_estimate", "seurat", "seurat_check", "seurat_clusters",
  "seurat_damaged", "seurat_mtrb", "spliced", "start_time", "sub_dirs", "subtitle",
  "table_of_counts", "table_of_droplets", "time_output", "title", "tmp_conn", "UMAP_markers",
  "unspliced", "wrap_plots"
))

preprocess_core <- function(
    project_name,
    organism = "Hsap",
    ptprc_threshold = 0,
    filtered_path,
    raw_path = NULL,
    velocyto_path = NULL,
    filter_damaged = TRUE,
    output_dir     = "./"
){

  # Begin timing
  start_time <- Sys.time()
  cat("\nBegin pre-processing for", project_name, "...\n")


  # Check if the input organism is valid
  if (!organism %in% c("Mmul", "Hsap", "Mmus")) {
    stop("Invalid organism. Must be one of 'Mmul', 'Hsap', 'Mmus'.")
  }

  # Read in the corresponding species gene annotation file using the specified organism (Mmul, Hsap, Mmus) for correct gene nomenclature
  if (organism == "Hsap") {annotations <- readRDS("/home/alicen/R/human_annotations.rds")}
  if (organism == "Mmul") {annotations <- readRDS("/home/alicen/R/macaque_annotations.rds")}
  if (organism == "Mmus") {annotations <- readRDS("/home/alicen/R/mouse_annotations.rds")}

  # Verify output directory structure & create if absent
  sub_dirs <- c("RBCQC", "IMCQC", "CellQC", "SampleClustering", "R_objects")

  for (sub_dir in sub_dirs) {
    full_path <- file.path(output_dir, sub_dir)

    if (!dir.exists(full_path)) {
      dir.create(full_path, recursive = FALSE)
    }
  }

  # Initialize the output df (if it's not already empty ~ allows for looping)
  if (is.null(preprocess_output)) {
    preprocess_output <- data.frame(
      Sample_name = character(),
      Initial_cells = numeric(),
      Rho = numeric(),
      RBCs = numeric(),
      Immune_cells = numeric(),
      Potentially_damaged_cells = numeric(),
      Doublets = numeric(),
      Final_cells = numeric(),
      Time = numeric()
    )
  }

  # 1. Ambient RNA correction with SoupX --------------

  if (is.null(raw_path) & is.null(velocyto_path))

  {

    # Using Seurat, read in the matrices from the alignment output
    table_of_counts <- suppressWarnings(Read10X(filtered_path))

    # Create a Seurat object & check the cell number of the sample
    seurat <- suppressWarnings(CreateSeuratObject(
      counts = table_of_counts,
      project = project_name,
      min.cells = 1,
      min.features = 50
    ))

    initial_cells <- length(Cells(seurat))
    rho_estimate <- "NA"

  }

  else # If raw & velocyto counts are available continue with SoupX

  {

    cat("\u2714  Ambient RNA correcting ...\n")

    # Using Seurat read in the matrices from the STARsolo output for filtered (TOC) and raw (TOD) counts  (must be zipped input files)
    table_of_droplets <- suppressWarnings(Read10X(raw_path))
    table_of_counts <- suppressWarnings(Read10X(filtered_path))

    # Create a temp Seurat object to check the cell number of the sample (obvs use the TOC)
    check <- suppressWarnings(CreateSeuratObject(
      counts = table_of_counts,
      project = project_name,
      min.cells = 1))


    seurat <- NULL # initialize
    initial_cells <- length(Cells(check))

    # Only run Soup X if the Seurat object more than 1000 cells (else problems)
    if (initial_cells >= 1000) {

      # Create the soup channel (sc)
      sc <- SoupChannel(table_of_droplets, table_of_counts, calcSoupProfile = FALSE)

      # Estimate the contamination
      sc <- estimateSoup(sc)

      # Use Seurat to cluster the filtered matrix, although not essential it is recommended to get better estimations
      seurat_soup <- suppressWarnings(CreateSeuratObject(table_of_counts, min.cells = 0))
      seurat_soup <- suppressWarnings(SCTransform(seurat_soup, verbose = FALSE) %>%
                                        RunPCA(verbose = FALSE) %>%
                                        RunUMAP(dims = 1:30, verbose = FALSE) %>%
                                        FindNeighbors(dims = 1:30, verbose = FALSE) %>%
                                        FindClusters(verbose = FALSE))

      # Adding the cluster embeddings to the SoupX object
      meta.data <- seurat_soup@meta.data
      umap.embedding <- seurat_soup@reductions$umap@cell.embeddings
      sc <- suppressWarnings(setClusters(sc, setNames(meta.data$seurat_clusters, rownames(meta.data))))
      sc <- suppressWarnings(setDR(sc, umap.embedding, c("UMAP_1", "UMAP_2")))

      # With defined clusters, run the main SoupX function to calculate the contamination fraction rho where rho E (0, 1) and the closer to 1, the more contaminated
      sc <- autoEstCont(sc, verbose = FALSE, doPlot = FALSE)
      rho_estimate <- sc$metaData$rho[1] # record the rho in PreProcess output

      # Silencing the output ...
      # Open a connection to a temporary file for writing
      tmp_conn <- file(tempfile(), open = "wt")

      # Redirect standard output and messages to the temporary file
      sink(tmp_conn)
      sink(tmp_conn, type = "message")

      # Call the (now silenced) verbose function
      # Output integer matrix of soup-corrected reads (unzipped output) where contaminated reads are removed
      adj.matrix <- suppressWarnings(adjustCounts(sc, roundToInt = T))

      # Reset output redirection
      sink(NULL)
      sink(NULL, type = "message")

      # Close the connection
      close(tmp_conn)

      # Output results in Seurat object for continued workflow
      seurat <- suppressWarnings(CreateSeuratObject(counts = adj.matrix, # SoupX corrected count matrix
                                                    min.cells = 1, # At least one cell must express the gene for the gene to be included in the count matrix
                                                    project = project_name))

      cat("\u2714  Correction complete\n")

    }

    else

    {
      # Otherwise, if there are fewer than 1000 cells use the adjusted count matrix
      seurat <- check

    }

    # Finalise Seurat object creation for next steps
    # Add the organism type as meta data (helpful for downstream)
    seurat$organism <- organism
    cat("\u2714  Initial cell number of", initial_cells, "\n")



    # 2. Add nuclear fraction data --------------


    # Calculate the nuclear fraction (nf) from spliced and unspliced counts
    spliced <- ReadMtx(mtx      = paste0(velocyto_path, "spliced.mtx.gz"),
                       cells    = paste0(velocyto_path, "barcodes.tsv.gz"),
                       features = paste0(velocyto_path, "features.tsv.gz"))

    spliced <- suppressWarnings(CreateSeuratObject(
      counts    = spliced,
      project   = "spliced",
      min.cells = 1))

    unspliced <- ReadMtx(mtx      = paste0(velocyto_path, "unspliced.mtx.gz"),
                         cells    = paste0(velocyto_path, "barcodes.tsv.gz"),
                         features = paste0(velocyto_path, "features.tsv.gz"))

    unspliced <- suppressWarnings(CreateSeuratObject(
      counts    = unspliced,
      project   = "unspliced",
      min.cells = 1))

    # Calculate nuclear fraction
    exon_sum   <- Matrix::colSums(spliced[['RNA']]$counts)   # summing over all the genes for each cell (1 reading per cell)
    intron_sum <- Matrix::colSums(unspliced[['RNA']]$counts)
    nuclear_fraction <- intron_sum / (exon_sum + intron_sum)
    nf_csv <- data.frame(barcode = rownames(unspliced@meta.data), nf = nuclear_fraction)

    # Add nf to Seurat object (ensure row order (barcode) correct)
    seurat@meta.data[['nf']] <- nf_csv$nf[match(rownames(seurat@meta.data), nf_csv$barcode)]

    cat("\u2714  Velocyto nuclear fraction output added\n")

  }

  # 3. Remove red blood cells and isolate immune cells --------------

  # Use the specified organism (Mmul, Hsap, Mmus) for correct gene nomenclature
  if (organism == "Hsap") {

    # Define genes of interest
    hemo_gene <- c("HBA1", "HBA2", "HBB")          # haemoglobin subunit genes
    cd45_gene <- c("PTPRC", "PTPRCAP")             # CD45 as proxy for immune cell & CD45 polypeptide-associated protein

  }

  if (organism == "Mmul") {

    # Define genes of interest
    hemo_gene <- c("HBM", "HBA", "HBQ1")           # haemoglobin subunit genes slightly different names
    cd45_gene <- c("PTPRC", "PTPRCAP")

  }

  if (organism == "Mmus") {

    # Define genes of interest
    hemo_gene <- c("Hba-a1", "Hba-a2", "Hbb-bt", "Hbb-bs") # Mmus genes are lower case & different
    cd45_gene <- c("Ptprc", "Ptprcap")

  }

  # Add metadata columns to store the expression values of the genes of interest in each barcode
  seurat[['hbb.percent']] <- PercentageFeatureSet(
    object   = seurat,
    features = intersect(hemo_gene, rownames(seurat@assays$RNA)),
    assay    = "RNA"
  )

  seurat[['ptprc.percent']] <- PercentageFeatureSet(
    object   = seurat,
    features = intersect(cd45_gene, rownames(seurat@assays$RNA)),
    assay    = "RNA"
  )

  # RBC
  # First, label the barcodes as either if they are likely red blood cells or have red blood cell contamination
  seurat[['RBC']] <- ifelse(seurat$hbb.percent !=  0, "RBC", "non-RBC")

  # Calculate the number of red blood cells (account for no hemoglobin expression)
  RBC_number <- NULL

  if (all(seurat$hbb.percent == 0)) {

    # Say there were no RBCs
    RBC_number = 0
    cat("\u2714  No red blood cells found\n")

  } else {

    # Count RBCs if they are present
    RBC <- subset(seurat, hbb.percent > 0)
    RBC_number <- length(Cells(RBC))
    cat("\u2714  Red blood cells removed\n")

  }

  # Plot the hemoglobin percent
  rbcDf <- as.data.frame(seurat@meta.data)
  RBC_percent <- ( RBC_number / initial_cells ) * 100

  # Plot hemo vs CD45
  hemo_feature_QC <- ggplot(rbcDf, aes(x = hbb.percent, y = ptprc.percent, color = RBC)) +
    geom_point(size = 0.6) +
    scale_color_manual(values = c("non-RBC" = "grey","RBC" = "#CF597B")) +
    xlab("Haemoglobin expression") +
    ylab("CD45 expression") +
    labs(title = project_name) +
    annotate("text", x = Inf, y = Inf, label = paste("Percentage of contamination: ", round(RBC_percent, 2), "%"), hjust = 1.1, vjust = 2.4) +
    theme_classic() +
    theme(
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.4), 'cm'),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title.x = element_text(hjust = 0.5, vjust = -1, face = "bold", size = 12),
      axis.title.y = element_text(hjust = 0.5, vjust = 3, face = "bold", size = 12),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = -1.2, size = 12),
      axis.text.y = element_text(angle = 0, hjust = 1, vjust = -1.2, size = 12),
      legend.title = element_text(face = "bold"),
      legend.position = "below",
      legend.box.background = element_rect(colour = "black"),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1)

    )

  # Save output  to directory specified (or default) : here we don't want the individual plots, we want them all together
  saveRDS(hemo_feature_QC, file.path(output_dir, "RBCQC", paste0(project_name, "_RBCQC", ".rds")))

  # Filtering step after visualization, keep those not marked as RBC
  seurat <- subset(seurat, RBC == "non-RBC")

  # IMC
  # Second, label the non-RBC barcodes as either IMC or non-IMC if they do not express CD45 (aka PTPRC)
  if (ptprc_threshold == "FALSE"){

    # If they don't want any cells to be filtered
    seurat[['IMC']] <- "IMC"
    IMC_number <- "NA"
    cat("\u2714  All cells retained\n")

  } else {

    # For any numerical value, such as 0, only barcodes with CD45 percentage above this (not equal to) will be marked as IMC
    seurat[['IMC']] <- ifelse(seurat$ptprc.percent > ptprc_threshold, "IMC", "non-IMC")

    # Calculate the number of immune cells (account for no immune expression)
    IMC <- subset(seurat, ptprc.percent > ptprc_threshold)
    IMC_number <- length(Cells(IMC))
    cat("\u2714  Total immune cells retained ", IMC_number, "\n")

  }


  # Plot the barcodes to visualise how many were removed (non-IMC)
  imcDf <- as.data.frame(seurat@meta.data)

  if (ptprc_threshold == "FALSE"){IMC_percent <- 0}
  else {IMC_percent <- ( IMC_number / initial_cells ) * 100 }


  # Plotting metrics that give good spread of data, UMI & gene counts per barcode
  immune_feature_QC <- ggplot(imcDf, aes(x = nCount_RNA, y = nFeature_RNA, color = IMC)) +
    geom_point(size = 0.6) +
    scale_color_manual(values = c("IMC" = "grey","non-IMC" = "#5372B4")) +
    labs(title = project_name) +
    annotate("text", x = Inf, y = Inf, label = paste("Percentage of immune cells: ", round(IMC_percent, 2), "%"), hjust = 1.1, vjust = 2.4) +
    theme_classic() +
    theme(
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.4), 'cm'),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title.x = element_text(hjust = 0.5, vjust = -1, face = "bold", size = 12),
      axis.title.y = element_text(hjust = 0.5, vjust = 3, face = "bold", size = 12),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = -1.2, size = 12),
      axis.text.y = element_text(angle = 0, hjust = 1, vjust = -1.2, size = 12),
      legend.title = element_text(face = "bold"),
      legend.position = "below",
      legend.box.background = element_rect(colour = "black"),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1)

    )

  # Save output  to directory specified (or default) : here we don't want the individual plots, we want them all together
  saveRDS(immune_feature_QC, file.path(output_dir, "IMCQC", paste0(project_name, "_IMCQC", ".rds")))

  # Filtering step after visualization, keep those that are immune cells
  seurat <- subset(seurat, IMC == "IMC")


  # Only continue if velocyto present
  if (!is.null(velocyto_path)) {

    # 4. Identify droplets to be removed with DropletQC --------------

    # Extract nf meta data & associated cell barcode from Seurat object
    edDf <- data.frame(nf = as.numeric(seurat$nf), umi = seurat$nCount_RNA)

    # Use DropletQC function to identify empty droplets
    edresultsDf <- DropletQC::identify_empty_drops(edDf)

    # Identify damaged cells
    # Create vector of length same as cell number (runs default with no adjustment on groups of cells, need to fill this column requirement while getting results for that sample as a whole)
    n_cells <- length(Cells(seurat))
    cell_type <- rep(1, n_cells)

    dcDf <- data.frame(
      nf = as.numeric(seurat$nf),
      umi = seurat$nCount_RNA,
      cell_status = edresultsDf$cell_status,
      cell_type = cell_type
    )

    # Open a connection to a temporary file for writing
    tmp_conn <- file(tempfile(), open = "wt")

    # Redirect standard output and messages to the temporary file
    sink(tmp_conn)
    sink(tmp_conn, type = "message")

    # Your function call
    dcresultsDf <- DropletQC::identify_damaged_cells(dcDf)

    # Reset output redirection
    sink(NULL)
    sink(NULL, type = "message")

    # Close the connection
    close(tmp_conn)

    # Add DropletQC output to Seurat object metadata
    seurat$DQCStatus <- dcresultsDf$df$cell_status

    # Rescue cells (dataset specific)
    if (organism == "Hsap") {
      seurat$DQCStatus <- ifelse(seurat$DQCStatus == "empty_droplet" & seurat$nf >= 0.02, "cell", seurat$DQCStatus)
      seurat$DQCStatus <- ifelse(seurat$DQCStatus == "empty_droplet" & seurat$nCount_RNA >= 2000, "cell", seurat$DQCStatus)
    }

    if (organism == "Mmul") {
      seurat$DQCStatus <- ifelse(seurat$DQCStatus == "empty_droplet" & seurat$nf >= 0.08, "cell", seurat$DQCStatus)
      seurat$DQCStatus <- ifelse(seurat$DQCStatus == "empty_droplet" & seurat$nCount_RNA >= 10000, "cell", seurat$DQCStatus)
    }


    cat("\u2714  Droplet QC predictions\n")

  }

  else

  {

    # All cells marked as cell
    seurat$DQCStatus <- "cell"

  }

  # 5. Identify damaged cells with MtRb --------------

  if (organism == "Hsap") {
    # Get gene annotations for mitochondrial (MT) & ribosomal (RB) genes
    mt_gene_annotations <- annotations[grep("MT-", annotations$gene_name, perl=TRUE),]
    mt_gene_annotations <- mt_gene_annotations[grepl("protein_coding", mt_gene_annotations$gene_biotype, perl=TRUE),]
    mt_genes <- mt_gene_annotations %>% pull(gene_name)

    # isolate ribosomal genes
    rps_genes <- annotations[grep("RPS", annotations$gene_name, perl=TRUE),]
    rps_genes <- rps_genes[grepl("protein_coding", rps_genes$gene_biotype, perl=TRUE),]
    rps_genes <- rps_genes %>% pull(gene_name)
    rpl_genes <- annotations[grep("RPL", annotations$gene_name, perl=TRUE),]
    rpl_genes <- rpl_genes[grepl("protein_coding", rpl_genes$gene_biotype, perl=TRUE),]
    rpl_genes <- rpl_genes %>% pull(gene_name)
    rb_genes  <- c(rps_genes, rpl_genes)

    # combine mt and rb genes
    mt_rb_genes <- c(mt_genes, rb_genes)
    mt_rb_genes <- unique(mt_rb_genes)
  }

  if (organism == "Mmul") {

    # Get gene annotations for mitochondrial (MT) & ribosomal (RB) genes
    mt_genes <- annotations[grepl("mitochondrial", annotations$description, perl=TRUE),]
    mt_genes <- mt_genes %>% pull(gene_name)

    # isolate ribosomal genes
    rb_genes <- annotations[grepl("ribosomal", annotations$description, perl=TRUE),]
    rb_genes <- rb_genes %>% pull(gene_name)

    # combine mt and rb genes
    mt_rb_genes <- c(mt_genes, rb_genes)
    mt_rb_genes <- unique(mt_rb_genes)

  }

  if (organism == "Mmus") {

    # Get gene annotations for mitochondrial genes
    mt_genes <- annotations[grep("mt-", annotations$gene_name, perl = TRUE), ]
    mt_genes <- mt_genes[grepl("mitochondrially encoded", mt_genes$description, perl = TRUE), ]
    mt_genes <- mt_genes %>% pull(gene_name)


    # isolate ribosomal genes
    rb_genes <- annotations[grepl("ribosomal", annotations$description, perl = TRUE), ]
    rb_genes <- rb_genes[grepl("protein_coding", rb_genes$gene_biotype, perl = TRUE), ]
    rb_genes <- rb_genes %>% pull(gene_name)

    # combine mt and rb genes
    mt_rb_genes <- c(mt_genes, rb_genes)
    mt_rb_genes <- unique(mt_rb_genes)

  }


  # reduce based on mt & rb genes only (this occurs in a separate Seurat object (temp))
  seurat_mtrb <- subset(seurat, features = intersect(mt_rb_genes, rownames(seurat@assays$RNA)))

  seurat_mtrb <- NormalizeData(seurat_mtrb, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = 1:10, verbose = FALSE) %>%
    FindClusters(resolution = 0.7, verbose = FALSE) %>% # high resolution hoping to finely locate population with as little inclusion of non-damaged cells as possible
    RunTSNE(dims = 1:10, verbose = FALSE, check_duplicates = FALSE) # In this case (with far fewer genes) tSNE is prefered

  # Annotations to the reduced object based on the unreduced seurat object (all genes)
  DefaultAssay(seurat_mtrb) <- "RNA"
  DefaultAssay(seurat) <- "RNA"

  # Define mitochondrial expression in the temp Seurat object using the counts from the actual object
  seurat_mtrb$mt.percent <- PercentageFeatureSet(
    object   =  seurat,
    features = intersect(mt_genes, rownames(seurat@assays$RNA)),
    assay    = "RNA"
  )

  # Define ribosomal expression
  seurat_mtrb$rb.percent <- PercentageFeatureSet(
    object   = seurat,
    features = intersect(rb_genes, rownames(seurat@assays$RNA)),
    assay    = "RNA"
  )

  # Calculate complexity
  matrix <- GetAssayData(seurat_mtrb, layer = "data")
  complexity_metric <- colSums(matrix > 0)
  seurat_mtrb$complexity <- complexity_metric
  seurat$complexity <- seurat_mtrb$complexity # Transfer to Seurat

  # Automatically find the damaged cell population
  # Calculate the average mt.percent and rb.percent for each cluster, add ranks, and identify the best cluster
  best_cluster <- seurat_mtrb@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(
      avg_mt_percent = mean(mt.percent, na.rm = TRUE),
      avg_rb_percent = mean(rb.percent, na.rm = TRUE)
    ) %>%
    mutate(
      rank_mt = rank(-avg_mt_percent),  # Negative for descending order (highest mt to be #1)
      rank_rb = rank(avg_rb_percent),
      combined_rank = rank_mt + rank_rb
    ) %>%
    arrange(combined_rank) %>%
    slice(1) %>%
    pull(seurat_clusters) %>%
    as.character()

  # Label all cells belonging to this cluster as "damaged
  seurat_mtrb$seurat_clusters <- ifelse(seurat_mtrb$seurat_clusters == best_cluster,'damaged_cell',seurat_mtrb$seurat_clusters)

  # Add Cell QC meta data to object
  # MtRb results
  seurat_mtrb$MtRb <- ifelse(seurat_mtrb$seurat_clusters == "damaged_cell", "damaged_cell", "cell")

  # Add cluster #s to actual seurat object
  seurat$MtRb <- seurat_mtrb$MtRb

  cat("\u2714  MtRb QC predictions\n")

  # 6. Combine cell QC output --------------

  # "Initialize an empty metadata column
  seurat$QC <- "cell" # all barcodes marked as cell by default at the start

  seurat$QC <- seurat@meta.data %>%
    dplyr::mutate(QC = case_when(
      DQCStatus == "cell" & MtRb == "cell" ~ "cell",
      DQCStatus == "damaged_cell" & MtRb == "damaged_cell" ~ "DQ_MtRb",
      DQCStatus == "damaged_cell" & MtRb == "cell" ~ "DQ",
      DQCStatus == "cell" & MtRb == "damaged_cell" ~ "MtRb",
      DQCStatus == "empty_droplet" ~ "empty"
    ))

  # Transfer labels to reduced object to visualise labels
  seurat_mtrb$QC <- seurat$QC

  # Label damaged cells & empty drops
  Unfiltered_cells <- length(Cells(seurat))
  seurat$quality <- ifelse(seurat$QC == "cell" | seurat$QC == "DQ" | seurat$QC == "mt", "retain", "discard") # only discarding empty drops (empty) & those marked as damaged by both metrics

  # Quantify damaged cells
  Final_cells <-  length(Cells(seurat))
  seurat_damaged <- subset(seurat, quality == "discard")
  Damaged_cells <- length(Cells(seurat_damaged))
  Damaged_percent <- ( Damaged_cells / initial_cells ) * 100

  # View labels in reduced space (MtRb metric)
  clusters <- DimPlot(
    seurat_mtrb, pt.size = 1, group.by = "QC", cols = c("cell" = "grey", "DQ_MtRb" = "#5372B4", "DQ" = "#D27324", "MtRb" = "#A2B2D6", "empty" = "#D98C25")) +
    labs(caption = "MtRb metric") + NoAxes() +
    theme(
      legend.text = element_text(size = 15),
      plot.title = element_blank(),
      plot.caption = element_text(hjust = 0.5, size = 16),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1)
    )

  # tSNE of mt percent
  mt_plot <- FeaturePlot(seurat_mtrb, features = c("mt.percent"), cols = c("#E1E1E1", "#5372B4"), pt.size = 1) +
    NoAxes() + labs(caption = "Mt") +
    theme(
      plot.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
      plot.caption = element_text(hjust = 0.5, size = 16))

  # tSNE of rb percent
  rb_plot <- FeaturePlot(seurat_mtrb, features = c("rb.percent"), cols = c("#E1E1E1", "#5372B4"), pt.size = 1) +
    NoAxes() + labs(caption = "Rb") +
    theme(
      plot.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
      plot.caption = element_text(hjust = 0.5, size = 16))

  # tSNE of complexity
  complexity_plot <- FeaturePlot(seurat_mtrb, features = c("complexity"), cols = c("#E1E1E1", "#5372B4"), pt.size = 1) +
    NoAxes() + labs(caption = "Complexity") +
    theme(
      plot.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
      plot.caption = element_text(hjust = 0.5, size = 16))

  # PLOT COMBINING
  title <- ggdraw() + draw_label(project_name, fontface = 'bold', x = 0.45, y = 0.1, hjust = 0.5, size = 20)
  subtitle <- ggdraw() + draw_label(paste("Estimated", round(Damaged_percent, 2), "% damaged of ", initial_cells, " cells"), x = 0.45, hjust = 0.5, size = 16)
  plots <- ((mt_plot | complexity_plot) / (rb_plot | clusters))

  # Combine everything into the final plot
  cell_QC <- plot_grid(
    title, subtitle, plots,
    ncol = 1,
    rel_heights = c(0.1, 0.1, 1)
  )

  cell_QC <- (mt_plot | complexity_plot) / (rb_plot | clusters)
  ggsave(file.path(output_dir, "CellQC", paste0(project_name, "_CellQC", ".png")), plot = cell_QC, width = 16, height = 12, dpi = 300)


  if (!is.null(velocyto_path)) {
    # View labels in scatter plot nf vs UMI_count (DropletQC metric)
    df <- seurat@meta.data
    df$nf <- as.numeric(df$nf)
    df$nCount_RNA <- as.numeric(df$nCount_RNA)

    scatter <- ggplot(df, aes(x = nf, y = nCount_RNA, color = QC)) +
      geom_point(size = 1) +
      scale_color_manual(values = c("cell" = "grey", "DQ_MtRb" = "#5372B4", "DQ" = "#D27324", "MtRb" = "#A2B2D6", "empty" = "#D98C25")) +
      coord_trans(y = "log10") +
      xlab("nf") + ylab("log10(UMI)") + labs(caption = "DropletQC metric") +
      theme(
        plot.title = element_text(hjust = 0, vjust = 4, size = 20, face = "bold"),
        plot.caption = element_text(hjust = 0.5,  size = 16),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(hjust = 0.5, size = 16),
        axis.text.y = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 16),
        legend.position = "below",
        legend.box.background = element_rect(colour = "black"),
        strip.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth =1)

      )

    # PLOT COMBINING
    title <- ggdraw() + draw_label(project_name, fontface = 'bold', x = 0.45, y = 0.1, hjust = 0.5, size = 20)
    subtitle <- ggdraw() + draw_label(paste("Estimated", round(Damaged_percent, 2), "% damaged of ", initial_cells, " cells"), x = 0.45, hjust = 0.5, size = 16)
    plots <- ((mt_plot | complexity_plot) / (rb_plot | clusters)) / scatter

    # Combine everything into the final plot
    cell_QC <- plot_grid(
      title, subtitle, plots,
      ncol = 1,
      rel_heights = c(0.1, 0.1, 1)
    )

    ggsave(file.path(output_dir, "CellQC", paste0(project_name, "_CellQC", ".png")), plot = cell_QC, width = 12, height = 14, dpi = 300)

  }



  # After visualising, remove damaged cells & empty drops if specified

  if (filter_damaged == TRUE) {

    Unfiltered_cells <- length(Cells(seurat))
    seurat <- subset(seurat, quality == "retain") # As above, only removing empty drops & those marked as damaged by both metrics
    Final_cells <-  length(Cells(seurat))
    Damaged_cells <- Unfiltered_cells - Final_cells

  }


  # 7. DoubletFinder --------------
  cat("\u2714  Doublet Finder running\n")

  # Prepare seurat object (required) for DropletQC
  seurat_DF <- NormalizeData(seurat, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = 1:10, verbose = FALSE) %>%
    FindClusters(verbose = FALSE) %>%
    RunUMAP(dims = 1:10, verbose = FALSE)

  # Open a connection to a temporary file for writing (output of DoubletFinder is RIDICULOUSLY verbose but lacks the ability to quieten it)
  tmp_conn <- file(tempfile(), open = "wt")

  # Redirect standard output and messages to the temporary file
  sink(tmp_conn)
  sink(tmp_conn, type = "message")

  # pK Identification (no ground-truth)
  sweep.res.list <- paramSweep(seurat_DF, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)

  pK <- bcmvn %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  pK <- as.numeric(as.character(pK[[1]]))

  # Homotypic Doublet Proportion Estimate
  annotations <- seurat_DF@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075 * nrow(seurat_DF@meta.data))  # Assuming 7.5% doublet formation rate (can be changed for dataset specificity)
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  nExp_poi.adj <- nExp_poi.adj - (nExp_poi.adj * (1 / 7))

  # Run DF itself
  seurat_DF <- doubletFinder(
    seurat_DF,
    PCs = 1:10,
    pN = 0.25,
    pK = pK,
    nExp = nExp_poi.adj,
    reuse.pANN = FALSE,
  )

  # Reset output redirection
  sink(NULL)
  sink(NULL, type = "message")

  # Close the connection
  close(tmp_conn)



  # Meta data for actual Seurat object
  df_col_name <- grep("^DF\\.classifications", colnames(seurat_DF@meta.data), value = TRUE)
  doublet_classifications <- seurat_DF@meta.data[[df_col_name]]
  seurat$DF <- doublet_classifications

  # Print number of damaged cells
  Cells_before_doublets <- length(Cells(seurat))
  seurat <- subset(seurat, DF == "Singlet")
  Remaining_cells <- length(Cells(seurat))
  Doublet_cells <- Cells_before_doublets - Remaining_cells
  cat("\u2714  Doublet Finder complete with ", Doublet_cells, " doublets detected\n")

  # 8. UMAP check to visualize markers --------------

  # No longer altering the same seurat object, this is just for visualisation of each sample
  seurat_check <- NormalizeData(seurat, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = 1:10, verbose = FALSE) %>%
    FindClusters(verbose = FALSE) %>%
    RunUMAP(dims = 1:10, verbose = FALSE)

  # Use the specified organism (Mmul, Hsap, Mmus) for correct gene nomenclature
  if (organism == "Hsap") {

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

  if (organism == "Mmul") {

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

  if (organism == "Mmus") {

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


  markers_plot <- FeaturePlot(
    seurat_check,
    features = SpeciesFeatures,
    cols = c("#F1F1F1", "#1F78B4")) +
    theme(
      axis.text.x = element_text(angle = 90),
      axis.title.x = element_blank(),
      plot.caption = element_blank()
    )

  for (plot in 1:length(markers_plot)) {
    markers_plot[[plot]] <- markers_plot[[plot]] +
      NoLegend() +
      NoAxes() +
      labs(title = FeaturesLabels[plot]) +
      theme(panel.border = element_rect(colour = "black", fill=NA, linewidth =1))
    #annotate("text", x = Inf, y = Inf, label = FeaturesLabels[plot], hjust = 1.5, vjust = 2.4)
  }

  colours_full <- c("#9FCBE1","#70B1D2", "#A1A9F5","#6C79F0", "#C5E7A7","#ABDC7E", "#ADDBB6", "#72BC7B", "#C9AFDF","#9D71B3",
                    "#9FCBE1","#70B1D2", "#A1A9F5","#6C79F0", "#C5E7A7","#ABDC7E", "#ADDBB6", "#72BC7B", "#C9AFDF","#9D71B3",
                    "#9FCBE1","#70B1D2", "#A1A9F5","#6C79F0", "#C5E7A7","#ABDC7E", "#ADDBB6", "#72BC7B", "#C9AFDF","#9D71B3")


  # Open a connection to a temporary file for writing (output of DoubletFinder is RIDICULOUSLY verbose but lacks the ability to quieten it)
  tmp_conn <- file(tempfile(), open = "wt")

  # Redirect standard output and messages to the temporary file
  sink(tmp_conn)
  sink(tmp_conn, type = "message")

  plot_clusters <- suppressWarnings(DimPlot(seurat_check,
                                            pt.size = 1,
                                            label = TRUE,
                                            label.box = TRUE,
                                            label.size = 5,
                                            label.color = "white",
                                            repel = TRUE) +
                                      scale_color_manual(values = colours_full) +
                                      scale_fill_manual(values = colours_full) +
                                      labs(title = paste(project_name, " clustering and marker visualisation", sep = ""), subtitle = paste("Cells: ", Final_cells)) +
                                      xlab("UMAP 1") +
                                      ylab("UMAP 2") +
                                      theme(
                                        plot.title = element_text(hjust = 0.5, size = rel(1), face = "bold"),
                                        plot.subtitle = element_text(hjust = 0.5, vjust = 1),
                                        panel.border = element_rect(colour = "black", fill=NA, linewidth =1)))

  # Reset output redirection
  sink(NULL)
  sink(NULL, type = "message")

  # Close the connection
  close(tmp_conn)

  UMAP_markers <- plot_clusters + wrap_plots(markers_plot) + plot_layout(ncol = 2)

  ggsave(filename = file.path(output_dir, "SampleClustering", paste0(project_name, "_Clustering",".png")), plot = UMAP_markers, width = 20, height = 9, dpi = 300)

  # 9. Clean & save output --------------

  # Filter columns
  columns <- c("nf", "hbb.percent", "ptprc.percent", "IMC", "RBC", "DQCStatus", "MtRb", "QC", "DF")

  # "ptprc.percent"
  for (column in columns){

    if (column %in% colnames(seurat@meta.data)) {
      seurat[[column]] <- NULL }
  }

  # Only remove the filtered column if the user wants the object to be filtered (then column in redundant: only "retain")
  if (filter_damaged) {

    seurat$quality <- NULL

  }


  # End the timing function
  end_time <- Sys.time()
  time <- (end_time - start_time)
  time_output <- as.character(round(time, digits = 2))
  # cat(paste("TIME:", time_output, "\n"))


  # Append the sample information from the PreProcess run to the output data frame
  new_row <- data.frame(
    Sample_name               = project_name,
    Initial_cells             = initial_cells,
    Rho                       = rho_estimate,
    RBCs                      = RBC_number,
    Immune_cells              = IMC_number,
    Potentially_damaged_cells = Damaged_cells,
    Doublets                  = Doublet_cells,
    Final_cells               = Final_cells,
    Time                      = time_output
  )

  preprocess_output <- rbind(preprocess_output, new_row)

  # Save the processed object itself
  saveRDS(seurat, file.path(output_dir, "R_objects", paste0(project_name, ".rds")))


  cat("\u2714  Complete\n")


  return(list(
    seuratObject = seurat,
    preprocessOutput = preprocess_output
  ))


}
