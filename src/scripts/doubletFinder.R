#!/usr/bin/env Rscript
# Author:		K. Beigel
# Original code: E. Reichenberger
# Date:			6.20.2024
# Purpose: Identify doublets from sample data. The script should be called on each individual sample (on original/raw data).
# Notes: Doublets are identified and the cell barcodes for doublets are written to a .txt file.

sample <- ''
project <- ''
starting_data <- ''
input_path <- ''
output_path <- ''
mito_cutoff <- ''
min_feature_threshold <- ''
components <- ''
organism <- ''
lib_path <- ''
processes <- ''

args <- commandArgs(trailingOnly = TRUE)

# test if there is at least n arguments: if not, return an error
nargs <- 11
if (length(args) < nargs) {
  stop(paste('At least ', nargs, 'arguments must be supplied.'), call. = FALSE)
} else if (length(args) == nargs) {
  sample <- args[1]
  project <- args[2]
  starting_data <- args[3]
  input_path <- args[4]
  output_path <- args[5]
  mito_cutoff <- args[6]
  min_feature_threshold <- args[7]
  components <- args[8]
  organism <- args[9]
  lib_path <- args[10]
  processes <- args[11]
}

# LOAD LIBRARIES
#--------------------------------------------------------------------
suppressMessages(library(future, lib.loc = lib_path))
suppressMessages(library(tidyverse, lib.loc = lib_path))
suppressMessages(library(ggplot2, lib.loc = lib_path))
suppressMessages(library(rcartocolor, lib.loc = lib_path))
suppressMessages(library(cowplot, lib.loc = lib_path))
suppressMessages(library(Seurat, lib.loc = lib_path))
suppressMessages(library(DoubletFinder, lib.loc = lib_path))
#--------------------------------------------------------------------

# PARALLEL w/ FUTURE + SET SEED
# --------------------------------------------------------------------
options(future.globals.maxSize = 20000 * 1024^2)
plan(multisession(workers = as.integer(processes)))

set.seed(42)
# --------------------------------------------------------------------

# ASSIGN GLOBAL VARIABLES
#--------------------------------------------------------------------
pN <- 0.25
#--------------------------------------------------------------------

# OUTPUT PATHS FOR FIGURES AND DATATABLES
#--------------------------------------------------------------------
# Define figure and tables folders under sample folder
fig_dir <- file.path(output_path, 'figures')
tbl_dir <- file.path(output_path, 'tables')
print(paste0('Output path for figures: ', fig_dir))
#--------------------------------------------------------------------

# MAKE OUTPUT DIRS
#--------------------------------------------------------------------
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE)
dir.create(tbl_dir, showWarnings = FALSE)
#--------------------------------------------------------------------

# MAKE THRESHOLDS NUMERIC
#--------------------------------------------------------------------
mito_cutoff = as.numeric(mito_cutoff)
min_feature_threshold = as.numeric(min_feature_threshold)
#--------------------------------------------------------------------

# FUNCTIONS #########################################################
# READ IN FILES & CREATE SEURAT OBJECT
#--------------------------------------------------------------------
make_seuratobj <- function(inpath, starting.data, project, sample)
{
  outs.dir <- ''

  # For IDing doublets, use the STARTING_DATA
  if (starting.data == 'matrix') {
    outs.dir <- paste0(inpath)
  }

  if (starting.data == 'cellranger') {
    outs.dir <- file.path(inpath, 'outs/filtered_feature_bc_matrix')
  }
  print(paste0('Reading in data from: ', outs.dir))
  sample.data <- Read10X(data.dir = outs.dir)
  so <- CreateSeuratObject(counts = sample.data,
    project = paste0(project, '_', sample),
    min.cells = 0,
    min.features = 0
  )

  return(so)
}
#--------------------------------------------------------------------

# GET MULTIPLET RATE FROM NUMBER RECOVERED CELLS
#--------------------------------------------------------------------
determine_multiplet_rate <- function(seurat.object)
{
  print('Determining multiplet rate.')
  cell.count <- length(seurat.object@meta.data[["orig.ident"]])

  doublet.predictor <- NULL # NOTE: will be based on # of recovered cells (cellranger, filtered)

  if (cell.count > 0 && cell.count <= 740)
  { doublet.predictor <- 0.004 }
  if (cell.count > 740 && cell.count <= 1400)
  { doublet.predictor <- 0.008 }
  if (cell.count > 1400 && cell.count <= 2400)
  { doublet.predictor <- 0.016 }
  if (cell.count > 2400 && cell.count <= 3400)
  { doublet.predictor <- 0.023 }
  if (cell.count > 3400 && cell.count <= 4400)
  { doublet.predictor <- 0.031 }
  if (cell.count > 4400 && cell.count <= 5400)
  { doublet.predictor <- 0.039 }
  if (cell.count > 5400 && cell.count <= 6400)
  { doublet.predictor <- 0.046 }
  if (cell.count > 6400 && cell.count <= 7400)
  { doublet.predictor <- 0.054 }
  if (cell.count > 7400 && cell.count <= 8400)
  { doublet.predictor <- 0.061 }
  if (cell.count > 8400 && cell.count <= 9400)
  { doublet.predictor <- 0.069 }
  if (cell.count > 9400) # && cell.count <= 10400
  { doublet.predictor <- 0.076 }

  print(paste0('Mutiplet rate for ', cell.count, ' cells: ', doublet.predictor))
  write.table(
    data.frame(cell_count = cell.count, doublet_rate = doublet.predictor),
    file = file.path(tbl_dir, paste0(project, '_', sample, '_', 'doubletfinder_stats.csv')),
    quote = FALSE, sep = ',', row.names = FALSE
  )

  return(doublet.predictor)
}
#--------------------------------------------------------------------

# FILTER SEURAT OBJECT
#--------------------------------------------------------------------
filter_seuratobj <- function(seurat.object, organism, components)
{
  if (organism == 'human')
  {
    seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, pattern = "^MT-")
  }

  if (organism == 'mouse')
  {
    seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, pattern = "^mt-")
  }

  print('Pre-processing Seurat object.')
  seurat.object <- subset(seurat.object, subset = nFeature_RNA > min_feature_threshold & percent.mt < mito_cutoff)
  seurat.object <- NormalizeData(seurat.object)
  seurat.object <- FindVariableFeatures(seurat.object, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(seurat.object)
  seurat.object <- ScaleData(seurat.object, features = all.genes)
  seurat.object <- RunPCA(seurat.object, features = VariableFeatures(object = seurat.object))

  seurat.object <- RunUMAP(seurat.object, dims = 1:components)
  seurat.object <- FindNeighbors(seurat.object, dims = 1:components)
  seurat.object <- FindClusters(seurat.object, resolution = 0.5)

  return(seurat.object)
}
#--------------------------------------------------------------------

# GET PK VALUE FROM SEURAT OBJECT
#--------------------------------------------------------------------
get_pk <- function(seurat.object, sample, project)
{
  print('Calculating pk.')
  # pK Identification (no ground-truth) --------------------------
  sweep.res.list <- paramSweep(seurat.object, PCs = 1:components, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)

  # Need to plot the default pK plot from find.pK() somewhere
  # otherwise it just makes a default img in the working dir
  png(file = paste0(fig_dir, '/', project, '_', sample, '_pk_sweep_plot_default.png'),
      width = 11, height = 8.5, res = 300, units = 'in')
  bcmvn <- find.pK(sweep.stats)
  dev.off()

  # Plot a nicer pK plot with the actual axis labels
  png(file = paste0(fig_dir, '/', project, '_', sample, '_pk_sweep_plot.png'),
      width = 11, height = 8.5, res = 300, units = 'in')
  print(ggplot(bcmvn, aes(x = pK, y = BCmetric)) + geom_point())
  dev.off()

  write.table(bcmvn,
              file = paste0(tbl_dir, '/', project, '_', sample, '_pk_values.txt'),
              append = FALSE, sep = '\t', dec = '.',
              quote = FALSE, row.names = FALSE, col.names = TRUE)

  # Method to get pK at maximum BCmetric:
  bcmvn.max.row <- bcmvn[which.max(bcmvn$BCmetric), ]  # may have mutliple max, potentially need bcmvn[which(bcmvn$BCmetric == max(bcmvn$BCmetric)),]
  optimal.pk <- as.numeric(as.character(bcmvn.max.row$pK)) # character to numeric to remove factoring/levels

  write.table(optimal.pk,
              file = paste0(tbl_dir, '/', project, '_', sample, '_optimal_pk_value.txt'),
              append = FALSE, sep = '\t', dec = '.',
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  return(optimal.pk)
}
#--------------------------------------------------------------------

# RUN DATA THROUGH DOUBLETFINDER
#--------------------------------------------------------------------
run_doubletfinder <- function(seurat.object, doublet.rate, pk.value)
{
  print('Running doubletFinder.')
  ## Homotypic Doublet Proportion Estimate ------------------------
  annotations <- seurat.object@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp.poi <- round(doublet.rate * nrow(seurat.object@meta.data))
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

  seurat.object <- doubletFinder(seurat.object, PCs = 1:components, pN = pN, pK = pk.value, nExp = nExp.poi, reuse.pANN = NULL, sct = FALSE)
  seurat.object <- doubletFinder(seurat.object, PCs = 1:components, pN = pN, pK = pk.value, nExp = nExp.poi.adj, reuse.pANN = NULL, sct = FALSE)

  return(seurat.object)
}
#--------------------------------------------------------------------

# WRITE DOUBLET IDs TO FILE
#--------------------------------------------------------------------
write_doublet_ids <- function(seurat.object, sample, project)
{
  # Get the names of the doublet columns
  doublet.col.list <- colnames(seurat.object@meta.data)[grepl("DF.classification", colnames(seurat.object@meta.data))] #there could be more than 1

  print('Plotting  doublets.')
  f.name <- paste0(fig_dir, '/', project, '_', sample, '_doublet.png')
  d.name <- paste0(tbl_dir, '/', project, '_', sample, '_doublet_ids.txt')

  if (file.exists(d.name)) {
    file.remove(d.name)
  }

  # If any of the DoubletFinder columsn classify a cell as a doublet,
  # define it as a Doublet in the final Doublet_Classification
  doublet.class = seurat.object@meta.data[, doublet.col.list] %>%
    mutate(DF_Classification = case_when(
		if_any(all_of(doublet.col.list), ~ .x == 'Doublet') ~ 'Doublet',  TRUE ~ 'Singlet')  # Check across columns # added 3.4.25
      #if_any(c(doublet.col.list), ~ .x == 'Doublet') ~ 'Doublet', .default = 'Singlet')
    )
  
  # Add final Doublet classification to the Seurat object
  seurat.object = AddMetaData(seurat.object, metadata = doublet.class$DF_Classification, col.name = 'DF_Classification')

  # ------- Plotting -------
  # Colors for doublet plotting
  # colors from: https://github.com/Nowosad/rcartocolor
  dim.pal = carto_pal(7, "ag_GrnYl")[c(1, 4)]
  
  # Set the vln.pal differently based on if there are less than 2 doublets
  # When there's not enough cells to make a violin, the colors shift
  vln.pal = ''
  if (nrow(seurat.object@meta.data %>% filter(DF_Classification == "Doublet")) < 2) {
    vln.pal = rev(dim.pal)
  } else {
    vln.pal = dim.pal
  }

  png(file = f.name, width = 11, height = 8.5, res = 300, units = 'in')
  umap = DimPlot(seurat.object, group.by = 'DF_Classification', cols = dim.pal) + NoAxes()
  vln = VlnPlot(seurat.object, features = 'nFeature_RNA', group.by = 'DF_Classification', cols = vln.pal, pt.size = 0.1)
  print(plot_grid(umap, vln, rel_widths = c(2.5, 1)))
  dev.off()
  # ------------------------

  # Get cell barcodes from cells ID'd as doublets (from seurat object metadata)
  doublet.ids <- rownames(seurat.object@meta.data[seurat.object@meta.data[, 'DF_Classification'] == "Doublet", ])

  # Single file of doublet IDs, write IDs to file
  print('Writing doublet IDs (cell barcodes) to file.')
  write.table(doublet.ids, file = d.name, quote = FALSE, sep = '', col.names = 'doublet_ids', row.names = FALSE)
  print('Doublet IDs written.')
}
#--------------------------------------------------------------------
#####################################################################


# FUNCTION CALLS
#--------------------------------------------------------------------
# MAKE SEURAT OBJECT, GET MULTIPLET RATE, FILTER SEU OBJ, GET PK VAL
S <- make_seuratobj(input_path, starting_data, project, sample)
doublet_predictor <- determine_multiplet_rate(S)
S <- filter_seuratobj(S, organism, components)
pk <- get_pk(S, sample, project)

# RUNNING DOUBLETFINDER TO IDENTIFY DOUBLETS, WRITING IDS TO FILE
S_df <- run_doubletfinder(S, doublet_predictor, pk)
write_doublet_ids(S_df, sample, project)
#--------------------------------------------------------------------
