#!/usr/bin/env Rscript
# Author:    K. Beigel
# Date:      9.30.2024
# Modified by: M. Brown
# Mod note: M. Brown add preprocessing directive, converted positional args to flags with input checks and a default lib path for compatibility with argparse and nextflow integration
# Modified by: E. Reichenberger
# Mod note: E. Reichenberger updated to start with unfiltered seurat object + to use per-sample thresholds from passed file for filtering
# Date:      8.14.2026

# set a default library path for optparse to give help
lib_path <- '/usr/local/lib/R/site-library'

args <- commandArgs(trailingOnly = TRUE)
tryCatch(
  expr = {
    suppressMessages(library("optparse", lib.loc = lib_path))
  },
  error = function(e) {
    if (dir.exists(tail(args, n=1))) lib_path <- tail(args, n=1) else stop("Valid library path specified at the end as a positional arg required.")
    # if lib path provided, override default
    suppressMessages(library("optparse", lib.loc = lib_path))
  }
)

option_list <- list(
  make_option(c("--seurat_obj"), type="character",
              help="Path to the unfiltered Seurat object (.qs2) from the QC metrics step"),
  make_option(c("--thresholds_file"), type="character",
              help="TSV with per-sample filtering thresholds (Sample, MIN_FEATURE_THRESHOLD, MAX_FEATURE_THRESHOLD, MITO_THRESHOLD, RIBO_THRESHOLD)"),
  make_option(c("--project"), type="character",
              help="Name of the project"),
  make_option(c("--seurat_file_name"), type="character",
              help="Name of the filtered Seurat object file to save")
)

library(tidyverse, lib.loc = lib_path)
library(rcartocolor, lib.loc = lib_path)
library(Seurat, lib.loc = lib_path)
library(qs2, lib.loc = lib_path)
library(patchwork, lib.loc = lib_path)
library(dplyr, lib.loc = lib_path)

# Since the last args is positional, object sticks the options in a separate key
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments=TRUE, args=args)
seurat_obj_path <- if (is.null(opt$options$seurat_obj) | !file.exists(opt$options$seurat_obj)) stop("Valid --seurat_obj is required. See --help for all opts") else opt$options$seurat_obj
thresholds_file_path <- if (is.null(opt$options$thresholds_file) | !file.exists(opt$options$thresholds_file)) stop("Valid --thresholds_file is required. See --help for all opts") else opt$options$thresholds_file
project <- if (is.null(opt$options$project)) stop("--project is required. See --help for all opts") else opt$options$project
seurat_file_name <- if (is.null(opt$options$seurat_file_name)) stop("--seurat_file_name is required. See --help for all opts") else opt$options$seurat_file_name

# CREATE DIRECTORIES
#--------------------------------------------------------------------
# Output
base_directory = file.path('data/endpoints', project, '/analysis/')
figure_dir = file.path(base_directory, '/figures')
rds_dir = file.path(base_directory, '/RDS')
table_dir = file.path(base_directory, '/tables')
report_dir = file.path(base_directory, '/report/qc_report')

for (dir in c(figure_dir, rds_dir, table_dir, report_dir)) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}
#--------------------------------------------------------------------

# QC PLOTS: PRE- AND POST-FILTERING
#--------------------------------------------------------------------
seu_qc_plots = function(seu.obj, plot.set, title, caption)
{  
  # colors from: https://github.com/Nowosad/rcartocolor
  palette = carto_pal(7, "ag_Sunset")
  n_samples = length(unique(seu.obj@meta.data$Sample))
  if (n_samples <= 2) {
    palette = palette[c(2, 4)]
  } else {
    palette = NULL
  }
  
  n_cells <- ncol(seu.obj)
  f <- paste0(figure_dir, '/', project, '_qc_', plot.set)
  png(file = paste0(f, '_vln.png'), width = 11, height = 8.5, res = 300, units = 'in')
  print(
    VlnPlot(
      seu.obj,
      group.by = 'Sample',
      features = c('nFeature_RNA', 'nCount_RNA', 'percent.mito', 'percent.ribo'),
      ncol = 4,
      cols = palette
    ) +
    plot_annotation(
      title = title,
      subtitle = paste0('Number of cells: ', n_cells),
      caption = caption,
      theme = theme(
        plot.title = element_text(size = 26, hjust = 0.5),
        plot.subtitle = element_text(size = 18, hjust = 0.5),
        plot.caption = element_text(size = 15, hjust = 0.5)
      )
    )
  )
  dev.off()

  plot1 <- FeatureScatter(
    seu.obj,
    group.by = 'Sample',
    feature1 = 'nCount_RNA',
    feature2 = 'percent.mito',
    cols = palette,
    shuffle = TRUE,
    pt.size = 0.3
  )
  plot2 <- FeatureScatter(
    seu.obj,
    group.by = 'Sample',
    feature1 = 'nCount_RNA',
    feature2 = 'nFeature_RNA',
    cols = palette,
    shuffle = TRUE,
    pt.size = 0.3
  )
  patchwork_1 = plot1 + plot2
  png(file = paste0(f, '_scatter1.png'), width = 11, height = 8.5, res = 300, units = 'in')
  print(
    patchwork_1 +
    plot_annotation(
      title = title,
      subtitle = paste0('Number of cells: ', n_cells),
      caption = caption,
      theme = theme(
        plot.title = element_text(size = 26, hjust = 0.5),
        plot.subtitle = element_text(size = 18, hjust = 0.5),
        plot.caption = element_text(size = 15, hjust = 0.5)
      )
    )
  )
  dev.off()

  plot3 <- FeatureScatter(
    seu.obj,
    group.by = 'Sample',
    feature1 = 'nCount_RNA',
    feature2 = 'percent.ribo', 
    cols = palette,
    shuffle = TRUE,
    pt.size = 0.3
  )
  plot4 <- FeatureScatter(
    seu.obj,
    group.by = 'Sample',
    feature1 = 'nFeature_RNA',
    feature2 = 'percent.ribo',
    cols = palette,
    shuffle = TRUE,
    pt.size = 0.3
  )
  patchwork_2 = plot3 + plot4
  png(file = paste0(f, '_scatter2.png'), width = 11, height = 8.5, res = 300, units = 'in')
  print(
    patchwork_2 +
    plot_annotation(
      title = title,
      subtitle = paste0('Number of cells: ', n_cells),
      caption = caption,
      theme = theme(
        plot.title = element_text(size = 26, hjust = 0.5),
        plot.subtitle = element_text(size = 18, hjust = 0.5),
        plot.caption = element_text(size = 15, hjust = 0.5)
      )
    )
  )
  dev.off()
}
#--------------------------------------------------------------------

# FILTER SEURAT OBJECT — per-sample thresholds, fail loudly on missing samples
#--------------------------------------------------------------------
filter_save_seuobj = function(seu.obj, thresholds)
{
  # pre-filtered plots
  print('Plotting QC figures for data (not filtered).')
  seu_qc_plots(seu.obj, '1', 'Unfiltered', 'Data before filtering')

  # every sample in the object must have a matching row in the thresholds table
  obj_samples <- unique(seu.obj@meta.data$Sample)
  missing_samples <- setdiff(obj_samples, thresholds$Sample)
  if (length(missing_samples) > 0) {
    stop(paste0(
      'No thresholds found for sample(s): ', paste(missing_samples, collapse = ', '),
      '. Check the thresholds file: ', thresholds_file_path
    ))
  }

  # join each cell's Sample to its per-sample thresholds
  meta <- seu.obj@meta.data
  meta$cell_id <- rownames(meta)
  meta_joined <- left_join(meta, thresholds, by = 'Sample')

  keep_cells <- meta_joined$cell_id[
    meta_joined$nFeature_RNA > meta_joined$MIN_FEATURE_THRESHOLD &
    meta_joined$nFeature_RNA < meta_joined$MAX_FEATURE_THRESHOLD &
    meta_joined$percent.mito < meta_joined$MITO_THRESHOLD &
    meta_joined$percent.ribo < meta_joined$RIBO_THRESHOLD
  ]

  # per-sample before/after cell counts, for the log
  before_counts <- table(seu.obj@meta.data$Sample)
  seu.obj.filt <- subset(seu.obj, cells = keep_cells)
  after_counts <- table(seu.obj.filt@meta.data$Sample)

  count_summary <- data.frame(
    Sample = names(before_counts),
    before = as.integer(before_counts),
    after = as.integer(after_counts[names(before_counts)])
  )
  print('Cells retained per sample after filtering:')
  print(count_summary)

  # post-filtered plots
  print('Plotting QC figures for data (filtered).')
  seu_qc_plots(seu.obj.filt, '2', 'Filtered', 'Data after per-sample filtering')

  # save filtered object
  print('Saving Seurat object.')
  qs2::qs_save(seu.obj.filt, file = seurat_file_name)
}
#--------------------------------------------------------------------

# FUNCTION CALLS
#--------------------------------------------------------------------
seu.obj <- qs2::qs_read(seurat_obj_path)
thresholds <- read.table(thresholds_file_path, sep = '\t', header = TRUE)

filter_save_seuobj(seu.obj, thresholds)
#--------------------------------------------------------------------