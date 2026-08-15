#!/usr/bin/env Rscript
# Author:    R. Reichenberger

# Purpose:   Generate QC metrics + plots, suggest filtering thresholds 
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
  make_option(c("--sample_file"), type="character",
              help="TSV with header, sample name in first column, condition in second column, path to starting data in third column"),
  make_option(c("--project"), type="character",
              help="Name of the project"),
  make_option(c("--organism"), type="character",
              help="Organism (human or mouse)"),
  make_option(c("--seurat_creation_source"), type="character",
              help="Use count matrices from cellranger or soupX"),
  make_option(c("--run_doubletfinder"), type="character",
              help="'y' to incorporate doubletFinder data'"),
  make_option(c("--mito_cutoff"), type="integer",
              help="Mitochondrial percentage cutoff for filtering as an integer"),
  make_option(c("--ribo_cutoff"), type="integer",
              help="Ribosomal cutoff for filtering as an integer"),
  make_option(c("--min_feature_threshold"), type="integer",
              help="Minimum number of features per cell for filtering as an integer"),
  make_option(c("--seurat_file_name"), type="character",
              help="Name of the Seurat object file to save"),
  make_option(c("--base_directory"), type="character",
              help="Name of output directory"),
  make_option(c("--threshold_file"), type="character",
              help="Path to write the suggested thresholds TSV")
)

library(tidyverse, lib.loc = lib_path)
library(rcartocolor, lib.loc = lib_path)
library(Seurat, lib.loc = lib_path)
library(qs2, lib.loc = lib_path)
library(patchwork, lib.loc = lib_path)
library(dplyr, lib.loc = lib_path)
library(ggplot2, lib.loc = lib_path)

# Since the last args is positional, object sticks the options in a separate key
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments=TRUE, args=args)
sample_file <- if (is.null(opt$options$sample_file) | !file.exists(opt$options$sample_file)) stop("Valid --sample_file is required. See --help for all opts") else opt$options$sample_file
project <- if (is.null(opt$options$project)) stop("--project is required. See --help for all opts") else opt$options$project
organism <- if (is.null(opt$options$organism)) stop("--organism is required. See --help for all opts") else tolower(opt$options$organism)
seurat_creation_source <- if (is.null(opt$options$seurat_creation_source)) stop("--seurat_creation_source is required. See --help for all opts") else opt$options$seurat_creation_source
run_doubletfinder <- if (is.null(opt$options$run_doubletfinder)) stop("--run_doubletfinder is required. See --help for all opts") else tolower(opt$options$run_doubletfinder)
mito_cutoff <- if (is.null(opt$options$mito_cutoff)) stop("--mito_cutoff is required. See --help for all opts") else opt$options$mito_cutoff
ribo_cutoff <- if (is.null(opt$options$ribo_cutoff)) stop("--ribo_cutoff is required. See --help for all opts") else opt$options$ribo_cutoff
min_feature_threshold <- if (is.null(opt$options$min_feature_threshold)) stop("--min_feature_threshold is required. See --help for all opts") else opt$options$min_feature_threshold
seurat_file_name <- if (is.null(opt$options$seurat_file_name)) stop("--seurat_file_name is required. See --help for all opts") else opt$options$seurat_file_name
base_directory <- if (is.null(opt$options$base_directory)) stop("--base_directory is required. See --help for all opts") else opt$options$base_directory
threshold_file <- if (is.null(opt$options$threshold_file)) stop("--threshold_file is required. See --help for all opts") else opt$options$threshold_file

# CREATE DIRECTORIES
#--------------------------------------------------------------------
# Output
#base_directory = file.path('data/endpoints', project, '/analysis/qc_metrics')
figure_dir = file.path(base_directory, '/figures')
table_dir = file.path(base_directory, '/tables')

for (dir in c(figure_dir, table_dir)) 
{
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}
#--------------------------------------------------------------------

# MITO + RIBOSOMAL PREFIX (based on organism)
#--------------------------------------------------------------------
mito = '^MT-'
ribo = '^RP[LS]'

if (tolower(organism) == 'mouse') {
  mito = '^mt-'
  ribo = '^Rp[ls]'
}
#--------------------------------------------------------------------

# MAKE THRESHOLDS NUMERIC
#--------------------------------------------------------------------
mito_cutoff = as.numeric(mito_cutoff)
ribo_cutoff = as.numeric(ribo_cutoff)
min_feature_threshold = as.numeric(min_feature_threshold)
#max_feature_threshold = 0 #as.numeric(max_feature_threshold)

# Standard threshold line styling — change once here, applies everywhere below
THRESH_COL_1  <- '#3B9AB2'  # nFeature threshold (99th pct)
THRESH_COL_2  <- '#E1AF00'  # mito threshold (MAD)
THRESH_LTY    <- 'dotted'
THRESH_LWD    <- 1.0
#--------------------------------------------------------------------

# IMPORT DATA, CREATE SEURAT OBJECT, ADD METADATA, CALCULATE MITO%
#--------------------------------------------------------------------
get_sample_list = function(sample.file)
{ 
  sample.list = read.table(sample.file, sep = '\t', header = TRUE)
}

# CREATE LIST FOR MERGING
#--------------------------------------------------------------------
make_list_seu.objs = function(sample.list, merge.data.source)
{
  seu.obj.list = vector(mode = 'list')

  for (row in 1:nrow(sample.list)) {

    sample = sample.list[row, 1]
    experiment = sample.list[row, 2]

    data.path <- ''
    data.tenx <- ''

    # determine location of files
    if (merge.data.source == 'matrix' || merge.data.source == 'soupX') {
      data.path = paste0('data/endpoints/', project, '/', sample, '/', merge.data.source, '/')
    }

    if (merge.data.source == 'cellranger') {
      data.path = paste0('data/endpoints/', project, '/', sample, '/', merge.data.source, '/outs/filtered_feature_bc_matrix/')
    }

    print(paste('Loading 10X data for', sample, 'from', data.path))
    data.tenx = Read10X(data.path)

    print(paste('Creating Seurat object for', sample))
    seu.obj = CreateSeuratObject(counts = data.tenx, project = project, min.cells = 0, min.features = 0)
    seu.obj = AddMetaData(seu.obj, metadata = experiment, col.name = 'Experiment')
    seu.obj = AddMetaData(seu.obj, metadata = sample, col.name = 'Sample')
    seu.obj[['percent.mito']] = PercentageFeatureSet(seu.obj, pattern = mito)
    seu.obj[['percent.ribo']] = PercentageFeatureSet(seu.obj, pattern = ribo)
    print(paste(project, sample, 'doublet_ids', sep = '_'))
    if (run_doubletfinder == 'y') {
      doublet.file = paste0('data/endpoints/', project, '/', sample, '/doubletFinder/tables/', project, '_', sample, '_', 'doublet_ids.txt')
      doublet.ids = read.table(doublet.file, header = TRUE)
      doublet.ids = doublet.ids[, 'doublet_ids']
      print(paste0('Number of doublets: ', length(doublet.ids)))
      if (length(doublet.ids) >= 1) {
        print(paste0('Removing ', length(doublet.ids), ' doublets.'))
        seu.obj <- seu.obj[, !colnames(seu.obj) %in% doublet.ids]
      }
    }

    seu.obj.list[[sample]] = seu.obj
  }

  return(seu.obj.list)
}
#--------------------------------------------------------------------

# MAKE A MERGED SEURAT OBJECT
#--------------------------------------------------------------------
make_merged_seu.obj = function(seu.obj.list)
{
  # create merged seurat object
  print(paste('Merging samples:', paste(names(seu.obj.list), collapse = ', ')))
  S.merged <- merge(seu.obj.list[[1]], y = c(seu.obj.list[2:length(seu.obj.list)]), add.cell.ids = names(seu.obj.list), project = project)

  return(S.merged)
}
#--------------------------------------------------------------------

# FOR ONE SAMPLE, PULL SINGLE OBJ FROM LIST
#--------------------------------------------------------------------
make_single_seu.obj = function(seu.obj.list)
{
  # create merged seurat object
  S.single = seu.obj.list[[1]]
  S.single = RenameCells(object = S.single, add.cell.id = names(seu.obj.list)[1])

  return(S.single)
}
#--------------------------------------------------------------------

# PER-SAMPLE DISTRIBUTION PLOTS (unfiltered)
#--------------------------------------------------------------------
plot_qc_distributions_by_sample <- function(seu.obj, plot.set, title, floor_nfeature = min_feature_threshold, percentile = 0.99, mito_mad_multiplier = 3, ribo_mad_multiplier = 3, samples_per_page = 12)
{
  meta <- seu.obj@meta.data
  all_samples <- sort(unique(meta$Sample))
  n_samples <- length(all_samples)

  f <- paste0(figure_dir, '/', project, '_qc_', plot.set)

  # ---- compute per-sample nFeature threshold (excluding likely-empty droplets) ----
  threshold_df <- meta %>%
    filter(nFeature_RNA > floor_nfeature) %>%
    group_by(Sample) %>%
    summarise(
      pct_upper = quantile(nFeature_RNA, percentile),
      .groups = 'drop'
    )

  # ---- compute per-sample mito MAD threshold (upper/high-mito direction only) ----
  mito_threshold_df <- meta %>% group_by(Sample) %>%
    summarise(
      mito_med = median(percent.mito),
      mito_mad = mad(percent.mito),
      mito_mad_upper = mito_med + mito_mad_multiplier * mito_mad,
      .groups = 'drop'
    )

  # ---- compute per-sample ribo MAD threshold (lower direction) ----
  ribo_threshold_df <- meta %>% group_by(Sample) %>%
  summarise(
    ribo_med = median(percent.ribo),
    ribo_mad = mad(percent.ribo),
    ribo_mad_upper = ribo_med + ribo_mad_multiplier * ribo_mad,
    .groups = 'drop'
  )

  # ---- split samples into pages ----
  n_pages <- ceiling(n_samples / samples_per_page)
  sample_pages <- split(all_samples, ceiling(seq_along(all_samples) / samples_per_page))

  get_ncol <- function(n) 
  {
    if (n <= 2) return(1)
    if (n <= 6) return(2)
    if (n <= 12) return(3)
    return(4)
  }

  metrics <- c('nFeature_RNA', 'nCount_RNA', 'percent.mito', 'percent.ribo')

  for (page_i in seq_len(n_pages))
  {
    page_samples <- sample_pages[[page_i]]
    meta_page <- meta %>% filter(Sample %in% page_samples)
    threshold_page <- threshold_df %>% filter(Sample %in% page_samples)
    mito_threshold_page <- mito_threshold_df %>% filter(Sample %in% page_samples)
    ribo_threshold_page <- ribo_threshold_df %>% filter(Sample %in% page_samples)

    palette <- carto_pal(7, 'ag_Sunset')
    if (length(page_samples) <= 2) {
      palette <- palette[c(2, 4)]
    } else {
      palette <- NULL
    }

    ncol_page <- get_ncol(length(page_samples))
    page_suffix <- if (n_pages > 1) paste0('_page', page_i) else ''
    n_rows <- ceiling(length(page_samples) / ncol_page)
    plot_height <- max(6, n_rows * 2.5)

    for (metric in metrics)
    {
      p <- ggplot(meta_page, aes(x = .data[[metric]], fill = Sample)) +
        geom_density(alpha = 0.5) +
        facet_wrap(~ Sample, scales = 'free_y', ncol = ncol_page) +
        labs(title = paste0(title, ': ', metric), x = metric, y = 'Density') +
        theme_minimal(base_size = 12) +
        theme(legend.position = 'none')

      if (!is.null(palette)) p <- p + scale_fill_manual(values = palette)

      if (metric == 'nFeature_RNA') 
      {
        p <- p +
          geom_vline(
            data = threshold_page,
            aes(xintercept = pct_upper),
            color = THRESH_COL_1, linetype = THRESH_LTY, linewidth = THRESH_LWD
          ) +
          geom_text(
            data = threshold_page,
            aes(x = pct_upper, y = Inf, label = paste0('99th pct: ', round(pct_upper, 0))),
            color = THRESH_COL_1, angle = 90, vjust = -0.5, hjust = 1.1, size = 2.5, inherit.aes = FALSE
          )
      }

      if (metric == 'percent.mito') 
      {
        p <- p +
        geom_vline(
          data = mito_threshold_page,
          aes(xintercept = mito_mad_upper),
          color = THRESH_COL_2, linetype = THRESH_LTY, linewidth = THRESH_LWD
        ) +
        geom_text(
          data = mito_threshold_page,
          aes(x = mito_mad_upper, y = Inf, label = paste0(mito_mad_multiplier, 'xMAD: ', round(mito_mad_upper, 1))),
          color = THRESH_COL_2, angle = 90, vjust = -0.5, hjust = 1.1, size = 2.5, inherit.aes = FALSE
        )
      }

      if (metric == 'percent.ribo') 
      {
        p <- p +
        geom_vline(
          data = ribo_threshold_page,
          aes(xintercept = ribo_mad_upper),
          color = THRESH_COL_2, linetype = THRESH_LTY, linewidth = THRESH_LWD
        ) +
        geom_text(
          data = ribo_threshold_page,
          aes(x = ribo_mad_upper, y = Inf, label = paste0(ribo_mad_multiplier, 'xMAD: ', round(ribo_mad_upper, 1))),
          color = THRESH_COL_2, angle = 90, vjust = -0.5, hjust = 1.1, size = 2.5, inherit.aes = FALSE
        )
      }

      png(file = paste0(f, '_density_', metric, page_suffix, '.png'), width = 10, height = plot_height, res = 300, units = 'in')
      print(p)
      dev.off()
    }

    # ---- scatter: nCount_RNA vs nFeature_RNA ----
    plot2 <- ggplot(meta_page, aes(x = nCount_RNA, y = nFeature_RNA, color = Sample)) +
      geom_point(size = 0.3, alpha = 0.5) +
      facet_wrap(~ Sample, ncol = ncol_page) +
      theme_minimal(base_size = 12) +
      theme(legend.position = 'none') +
      geom_hline(
        data = threshold_page,
        aes(yintercept = pct_upper),
        color = THRESH_COL_1, linetype = THRESH_LTY, linewidth = THRESH_LWD
      ) +
      geom_text(
        data = threshold_page,
        aes(x = -Inf, y = pct_upper, label = paste0('99th pct: ', round(pct_upper, 0))),
        color = THRESH_COL_1, hjust = -0.1, vjust = -0.5, size = 2.5, inherit.aes = FALSE
      )
  
    if (!is.null(palette)) plot2 <- plot2 + scale_color_manual(values = palette)

    png(file = paste0(f, '_scatter_bysample', page_suffix, '.png'), width = 10, height = plot_height, res = 300, units = 'in')
    print(plot2)
    dev.off()

    # ---- scatter: percent.mito vs nFeature_RNA ----
    plot3 <- ggplot(meta_page, aes(x = nFeature_RNA, y = percent.mito, color = Sample)) +
      geom_point(size = 0.3, alpha = 0.5) +
      facet_wrap(~ Sample, ncol = ncol_page) +
      theme_minimal(base_size = 12) +
      theme(legend.position = 'none') +
      geom_vline(
        data = threshold_page,
        aes(xintercept = pct_upper),
        color = THRESH_COL_1, linetype = THRESH_LTY, linewidth = THRESH_LWD
      ) +
      geom_text(
        data = threshold_page,
        aes(x = pct_upper, y = Inf, label = paste0('99th pct: ', round(pct_upper, 0))),
        color = THRESH_COL_1, angle = 90, vjust = -0.5, hjust = 1.1, size = 2.5, inherit.aes = FALSE
      ) +
      geom_hline(
        data = mito_threshold_page,
        aes(yintercept = mito_mad_upper),
        color = THRESH_COL_2, linetype = THRESH_LTY, linewidth = THRESH_LWD
      ) +
      geom_text(
        data = mito_threshold_page,
        aes(x = -Inf, y = mito_mad_upper, label = paste0(mito_mad_multiplier, 'xMAD: ', round(mito_mad_upper, 1))),
        color = THRESH_COL_2, hjust = -0.1, vjust = -0.5, size = 2.5, inherit.aes = FALSE
      )

    if (!is.null(palette)) plot3 <- plot3 + scale_color_manual(values = palette)

    png(file = paste0(f, '_scatter_mito_vs_nfeature', page_suffix, '.png'), width = 10, height = plot_height, res = 300, units = 'in')
    print(plot3)
    dev.off()
  }

  # ---- combine both threshold tables into one summary ----
  full_threshold_df <- threshold_df %>%
  left_join(mito_threshold_df, by = 'Sample') %>%
  left_join(ribo_threshold_df, by = 'Sample')

  return(full_threshold_df)
}
#--------------------------------------------------------------------

# CORRELATION SCATTER PLOTS (mito/ribo vs features/counts, whole dataset)
#--------------------------------------------------------------------
plot_qc_correlations <- function(seu.obj, plot.set, title)
{
  meta <- seu.obj@meta.data
  f <- paste0(figure_dir, '/', project, '_qc_', plot.set)

  palette <- carto_pal(7, 'ag_Sunset')
  n_samples <- length(unique(meta$Sample))
  if (n_samples <= 2) palette <- palette[c(2, 4)] else palette <- NULL

  make_scatter <- function(x, y)
  {
    p <- ggplot(meta, aes(x = .data[[x]], y = .data[[y]], color = Sample)) +
      geom_point(size = 0.3, alpha = 0.5) +
      theme_minimal(base_size = 12)
    if (!is.null(palette)) p <- p + scale_color_manual(values = palette)
    p
  }

  # ---- mito: nFeature vs mito, nCount vs mito, nCount vs nFeature ----
  mito_combined <- make_scatter('nFeature_RNA', 'percent.mito') +
    make_scatter('nCount_RNA', 'percent.mito') +
    make_scatter('nCount_RNA', 'nFeature_RNA') +
    plot_annotation(
      title = paste0(title, ': percent.mito correlations'),
      subtitle = paste0('Number of cells: ', ncol(seu.obj))
    )

  png(file = paste0(f, '_mito_correlations.png'), width = 15, height = 5, res = 300, units = 'in')
  print(mito_combined)
  dev.off()

  # ---- ribo: nFeature vs ribo, nCount vs ribo, nCount vs nFeature ----
  ribo_combined <- make_scatter('nFeature_RNA', 'percent.ribo') +
    make_scatter('nCount_RNA', 'percent.ribo') +
    make_scatter('nCount_RNA', 'nFeature_RNA') +
    plot_annotation(
      title = paste0(title, ': percent.ribo correlations'),
      subtitle = paste0('Number of cells: ', ncol(seu.obj))
    )

  png(file = paste0(f, '_ribo_correlations.png'), width = 15, height = 5, res = 300, units = 'in')
  print(ribo_combined)
  dev.off()
}
#--------------------------------------------------------------------

# BIMODALITY / VALLEY DETECTION HELPERS
#--------------------------------------------------------------------
bimodality_coefficient <- function(x)
{
  x <- x[!is.na(x)]
  n <- length(x)
  m <- mean(x)
  s <- sd(x)

  skewness <- (sum((x - m)^3) / n) / s^3
  kurtosis <- (sum((x - m)^4) / n) / s^4

  bc <- (skewness^2 + 1) / (kurtosis + (3 * (n - 1)^2) / ((n - 2) * (n - 3)))
  return(bc)
}

find_valleys <- function(x, ...)
{
  x <- x[!is.na(x)]
  dens <- density(x, ...)
  d1 <- diff(dens$y)
  minima_idx <- which(diff(sign(d1)) == 2) + 1
  data.frame(x = dens$x[minima_idx], y = dens$y[minima_idx])
}

nfeature_bimodality_check <- function(seu.obj, metric)
{
  samples <- unique(seu.obj@meta.data$Sample)

  do.call(rbind, lapply(samples, function(s) {
    x <- seu.obj@meta.data[seu.obj@meta.data$Sample == s, metric]
    bc <- bimodality_coefficient(x)
    valleys <- find_valleys(x)

    data.frame(
      Sample = s,
      metric = metric,
      n = length(x),
      bimodality_coef = bc,
      likely_bimodal = bc > 0.555,
      n_valleys = nrow(valleys),
      valley_x = if (nrow(valleys) > 0) paste(round(valleys$x, 0), collapse = ', ') else NA
    )
  }))
}
#--------------------------------------------------------------------

# qc summary table output
#--------------------------------------------------------------------
write_qc_summaries <- function(bimodality_summary, output_thresholds, seu.obj, table_dir, project, threshold_file)
{
  write.table(
    bimodality_summary,
    file = file.path(table_dir, paste0(project, '_nFeature_bimodality_summary.tsv')),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
  )

  # Save per-sample threshold suggestions for use in filtering step
  #--------------------------------------------------------------------
  threshold_table <- output_thresholds %>%
    mutate(
      MIN_FEATURE_THRESHOLD = 200,
      MAX_FEATURE_THRESHOLD = round(pct_upper, 0),
      MITO_THRESHOLD = round(mito_mad_upper, 1),
      RIBO_THRESHOLD = round(ribo_cutoff, 1),
      RIBO_MAD_VALUE = round(ribo_mad_upper, 1)
    ) %>%
    select(Sample, MIN_FEATURE_THRESHOLD, MAX_FEATURE_THRESHOLD, MITO_THRESHOLD, RIBO_THRESHOLD, RIBO_MAD_VALUE)

  write.table(threshold_table, file = threshold_file, sep = '\t', row.names = FALSE, quote = FALSE)

  print(paste0('Suggested thresholds written to: ', threshold_file))

  # Compare sequencing depth (nCount_RNA) between samples
  #--------------------------------------------------------------------
  depth_summary <- seu.obj@meta.data %>%
    group_by(Sample) %>%
    summarise(
      n_cells        = n(),
      median_nCount  = median(nCount_RNA),
      mean_nCount    = mean(nCount_RNA),
      median_nFeature = median(nFeature_RNA),
      .groups = 'drop'
    )

  print(depth_summary)
}
#--------------------------------------------------------------------

# FUNCTION CALLS
#--------------------------------------------------------------------
sample_list = get_sample_list(sample_file)
seu.objs_list = make_list_seu.objs(sample_list, seurat_creation_source)

seu.obj = ''

if (nrow(sample_list) > 1) 
{
  seu.obj = make_merged_seu.obj(seu.objs_list)
}

if (nrow(sample_list) == 1) 
{
  seu.obj = make_single_seu.obj(seu.objs_list)
}

qs2::qs_save(seu.obj, file = seurat_file_name)

output_thresholds <- plot_qc_distributions_by_sample(seu.obj, '1', 'Unfiltered')
plot_qc_correlations(seu.obj, '1', 'Unfiltered')
bimodality_summary <- nfeature_bimodality_check(seu.obj, 'nFeature_RNA')
write_qc_summaries(bimodality_summary, output_thresholds, seu.obj, table_dir, project, threshold_file)
#--------------------------------------------------------------------
