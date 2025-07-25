#!/usr/bin/env Rscript
# Author: K. Beigel
# Date: 7.22.2024
# Purpose: Processing of Seurat object. This script will handle either an individual sample Seurat object or merged Seurat.
# Details: Runs normalization, integration (if more than once sample), PCA, UMAP.

initial_seurat_object = ''
project = ''
organism = ''
n_components = ''
mito_regression = ''
ribo_regression = ''
cc_regression = ''
cc_method = ''
num_var_features = ''
scale_data_features = ''
split_layers_by = ''
normalization_config = ''
integration_config = ''
ref_based_integration = ''
ref_samples = ''
run_azimuth = ''
azimuth_ref = ''
run_transferdata = ''
transferdata_ref_file = ''
transferdata_reduction = ''
transferdata_annocol = ''
resolution_config = ''
include_tsne = ''
analyzed_seurat_object = ''
report_path_figures = ''
lib_path = ''
processes = ''
memory_file = ''
max_memory_config = ''
regression_file = ''

args = commandArgs(trailingOnly = TRUE)

nargs = 30
if (length(args) != nargs) {
  stop(paste(nargs, 'arguments must be supplied.'), call. = FALSE)
} else if (length(args) == nargs) {
  initial_seurat_object = args[1]
  project = args[2]
  organism = args[3]
  n_components = args[4]
  mito_regression = args[5]
  ribo_regression = args[6]
  cc_regression = args[7]
  cc_method = args[8]
  num_var_features = args[9]
  scale_data_features = args[10]
  split_layers_by = args[11]
  normalization_config = args[12]
  integration_config = args[13]
  ref_based_integration = args[14]
  ref_samples = args[15]
  run_azimuth = args[16]
  azimuth_ref = args[17]
  run_transferdata = args[18]
  transferdata_ref_file = args[19]
  transferdata_reduction = args[20]
  transferdata_annocol = args[21]
  resolution_config = args[22]
  include_tsne = args[23]
  analyzed_seurat_object = args[24]
  report_path_figures = args[25]
  lib_path = args[26]
  processes = args[27]
  memory_file = args[28]
  max_memory_config = args[29]
  regression_file = args[30]
}

# LOAD PACKAGES
# --------------------------------------------------------------------
suppressMessages(library(future, lib.loc = lib_path))
suppressMessages(library(tools, lib.loc = lib_path))
suppressMessages(library(qs, lib.loc = lib_path))
suppressMessages(library(Seurat, lib.loc = lib_path))
suppressMessages(library(clustree, lib.loc = lib_path))
if (run_azimuth == 'y')
{
  suppressMessages(library(SeuratData, lib.loc = lib_path))
  suppressMessages(library(Azimuth, lib.loc = lib_path))
}
# --------------------------------------------------------------------

# GET MEMORY DEMAND
# --------------------------------------------------------------------
memory_f <- read.table(memory_file, header=FALSE)
memory <- as.numeric(memory_f$V1)

#memory <- memory_f
#
#if (memory_f > as.numeric(max_memory_config))
#{
#	memory <- as.numeric(max_memory_config)
#}
# --------------------------------------------------------------------

# PARALLEL w/ FUTURE + SET SEED
# --------------------------------------------------------------------
print(paste0('Assigning ', memory, ' bytes of memory...'))
options(future.globals.maxSize = memory)
#options(future.globals.onReference = 'error')
plan(multisession(workers = as.integer(processes)))
#plan(multisession(workers = 1)) # if higher number of samples and sct.rpca, runs out of memory (future issue)

set.seed(42)
# --------------------------------------------------------------------

# FUNCTION: PARSE COMMA-SEP STRINGS TO LISTS
# --------------------------------------------------------------------
single_or_list = function(variable, comma = ',')
{
  if (grepl(comma, variable, fixed = TRUE) == TRUE)
  {
    variable = strsplit(variable, comma)[[1]]
    return(variable)
  } else {
    return(variable)
  }
}
# --------------------------------------------------------------------

# PARSING CONFIG STRINGS
# --------------------------------------------------------------------
# PARSE THE LIST ARGUMENTS
normalization_config_list = single_or_list(normalization_config)
integration_config_list = single_or_list(integration_config)
resolution_config_list = as.numeric(single_or_list(resolution_config))
ref_samples = single_or_list(ref_samples)
# --------------------------------------------------------------------

# ASSIGN VARIABLES
# --------------------------------------------------------------------
# Normalization: list of which assay each config term for normalization will need
normalization_assay_dict = list('standard' = 'RNA', 'sct' = 'SCT')

# Subset to a ist of the norm methods from config
normalization_method_list = normalization_assay_dict[normalization_config_list]

# INTEGRATION METHODS: list of what each config term for integration will be
integration_dict = list('cca' = 'CCAIntegration', 'rpca' = 'RPCAIntegration', 'harmony' = 'HarmonyIntegration',
                        'pca' = 'PCA') # need this here so this passes properly to fns that use the integration_method_list

# Subset to a ist of the input methods from config with Seurat int names
integration_method_list = integration_dict[integration_config_list]
# --------------------------------------------------------------------

# OUTPUT PATHS
#--------------------------------------------------------------------
# Output that will be unchanged by component #s go here
base_directory = paste0('data/endpoints/', project, '/analysis')
figure_dir = paste0(base_directory, '/figures')
rds_dir = paste0(base_directory, '/RDS')
table_dir = paste0(base_directory, '/tables')

# Create directories
dir.create(rds_dir, showWarnings = FALSE)
dir.create(table_dir, showWarnings = FALSE)
dir.create(figure_dir, showWarnings = FALSE)
dir.create(report_path_figures, recursive = TRUE, showWarnings = FALSE)

# If there is more than one resolution, make dir for clustree fig(s)
if (length(resolution_config_list) > 1)
{
  dir.create(paste0(report_path_figures, '/clustree'), showWarnings = FALSE)
}
#--------------------------------------------------------------------

# FUNCTIONS
######################################################################
# IMPORT DATA
# --------------------------------------------------------------------
import_data = function(filename, filetype)
{
  print('Loading Seurat object.')
  seurat.object = qread(file = filename)

  return(seurat.object)
}
# --------------------------------------------------------------------

# GET A VECTOR OF THE VARIABLES (IN SEURAT METADATA) TO BE USED FOR REGRESSION
# --------------------------------------------------------------------
get_regression_variables = function(cc.regression, cc.method, mito.regression, ribo.regression, gene.module.file)
{
  regression.vars = NULL
  # If mito regression is specified, add to list of regression varaibles
  if (mito.regression == 'y')
  {
    regression.vars = append(regression.vars, 'percent.mito')
  }

  # If ribo regression is specified, add to list of regression varaibles
  if (ribo.regression == 'y')
  {
    regression.vars = append(regression.vars, 'percent.ribo')
  }

  # If cell cycle regression is specified, add to list of regression varaibles
  if (cc.regression == 'y')
  {
    if (cc.method == 'standard')
    {
      regression.vars = append(regression.vars, c('S.Score', 'G2M.Score'))
    }
    if (cc.method == 'alternative')
    {
      regression.vars = append(regression.vars, 'CC.Difference')
    }
  }

  if (gene.module.file != 'does_not_exist' & file.exists(gene.module.file))
  {
    regression.vars = append(regression.vars, 'gene.module.score1')
  }

  return(regression.vars)
}
# --------------------------------------------------------------------

# CONVERT LIST OF GENES FROM HUMAN TO MOUSE MGI SYMBOLS
# --------------------------------------------------------------------
# human.hgnc.list should be a list of genes (human HGNC symbols)
# https://github.com/satijalab/seurat/issues/2493#issuecomment-1334132532
convert_genes_hs2mm = function(human.hgnc.list)
{
  suppressMessages(library('gprofiler2', lib.loc = lib_path))
  mouse.mgi.list = gorth(human.hgnc.list, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

  return(mouse.mgi.list)
}
# --------------------------------------------------------------------

# DETERMINE CELL CYCLE SCORING (& OPTIONAL GENE LIST MODULE SCORING)
# --------------------------------------------------------------------
get_scoring = function(seurat.object, s.genes, g2m.genes, regression.vars, gene.module.file)
{
  set.seed(42)
  print('Joining layers of Seurat object.')
  # Join layers for CellCycleScoring to work in Seurat v5
  seurat.object = JoinLayers(seurat.object)

  print('Normalizing data for cell cycle scoring.')
  seurat.object = NormalizeData(seurat.object)

  print('Calculating cell cycle scoring.')
  seurat.object = CellCycleScoring(seurat.object, assay = 'RNA', s.features = s.genes, g2m.features = g2m.genes)
  seurat.object@meta.data[['CC.Difference']] = seurat.object@meta.data[['S.Score']] - seurat.object@meta.data[['G2M.Score']]

  pdf(file = paste0(figure_dir, '/', project, '_cell_cycle_scoring.pdf'), width = 11, height = 8.5)
  print(VlnPlot(seurat.object, features = c('S.Score', 'G2M.Score', 'CC.Difference'), group.by = 'Phase', ncol = 3, pt.size = 0.1) + ggtitle('Cell Cycle Phase'))
  dev.off()

  # Gene module scoring, if we added 'gene.module.score1' to the regression.vars
  if ('gene.module.score1' %in% regression.vars)
  {
    print('Loading list of genes to use for module scoring.')
    gene.module = as.vector(unique(read.table(gene.module.file, header = FALSE)$V1))

    print('Calculating gene list module scoring.')

    mod.length = length(gene.module)
    if (mod.length == 1)
    {
      mod.length = 2
    }

    seurat.object = AddModuleScore(object = seurat.object, features = list(gene.module), ctrl = mod.length, name = 'gene.module.score')
  }

  return(seurat.object)
}
# --------------------------------------------------------------------

# NORMALIZE DATA
# --------------------------------------------------------------------
seurat_processing = function(seurat.object, normalization.method.list, split.layers.by, regression.vars, num.var.features)
{
  set.seed(42)
  print('---- Processing Seurat object. ----')
  if (n_samples > 1) {
    # Split the object into layers according to whatever is specificed by split.layers.by
    print(paste0('Splitting Seurat object by ', split.layers.by, '.'))
    seurat.object[["RNA"]] = split(seurat.object[["RNA"]], f = get(split.layers.by, seurat.object@meta.data))
  }

  # Normalization (Standard Approach)
  if ('standard' %in% names(normalization.method.list)) {

    print('Running standard approach normalization.')

    # Note that NormalizeData() was run on the object with layers joined if get.cc.scoring() was run
    # Running it again here with layers split
    #seurat.object = NormalizeData(seurat.object)
    seurat.object = FindVariableFeatures(seurat.object, selection.method = 'vst', nfeatures = as.numeric(num.var.features))

    # get variable genes
    variable_genes <- VariableFeatures(seurat.object)

    variable_file <- paste0(table_dir, '/', project, '_standard_variable_genes.txt')
    write.table(variable_genes, file=variable_file, sep='\t', quote=F, col.names=TRUE, row.names=TRUE, append=FALSE)

    # Determine which features to use for scaling data based on config
    if (scale_data_features == 'all')
    {
      scale.features = rownames(seurat.object)
    }
    if (scale_data_features == 'variable')
    {
      scale.features = VariableFeatures(seurat.object)
    }

    # If the list of regression variables is null, no regression
    if (is.null(regression.vars))
    {
      seurat.object = ScaleData(seurat.object, features = scale.features)
    }
    # If regression variable list is not null, use those vars as vars.to.regress
    if (!is.null(regression.vars))
    {
      seurat.object = ScaleData(seurat.object, features = scale.features, vars.to.regress = regression.vars)
    }
  }

  # Normalization (SCTransform)
  if ('sct' %in% names(normalization.method.list))
  {
    print('Running SCT normalization.')

    # Determine which features to use for scaling data based on config
    if (scale_data_features == 'all')
    {
      num.scale.features = as.numeric(nrow(seurat.object))
    }
    if (scale_data_features == 'variable')
    {
      num.scale.features = as.numeric(num.var.features)
    }

    if (is.null(regression.vars))
    {
      seurat.object = SCTransform(seurat.object, method = "glmGamPoi", variable.features.n = num.scale.features)
    }
    if (!is.null(regression.vars))
    {
      seurat.object = SCTransform(seurat.object, method = "glmGamPoi", variable.features.n = num.scale.features, vars.to.regress = regression.vars)
    }

    sct.vargenes = VariableFeatures(seurat.object, assay = 'SCT')

    variable_sct_file <- paste0(table_dir, '/', project, '_sct_variable_genes.txt')
    write.table(sct.vargenes, file=variable_sct_file, sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE, append=FALSE)
  }

  # Running PCA and plotting ElbowPlot for each method
  for (normalization.method in names(normalization.method.list))
  {
    print(paste('PCA for', normalization.method, 'normalization.'))

    norm.assay = normalization.method.list[[normalization.method]]
    pca.name = paste0(normalization.method, '.pca')
    pc.key = paste0(normalization.method, 'PC_')

    # Run PCA
    seurat.object = RunPCA(seurat.object, assay = norm.assay, reduction.name = pca.name, npcs = as.numeric(n_components), reduction.key = pc.key)

    pdf(file = paste0(figure_dir, '/', project, '_', normalization.method, '_elbow_plot.pdf'), width = 11, height = 8.5)
    print(ElbowPlot(seurat.object, ndims = n_components, reduction = pca.name))
    dev.off()
  }

  return(seurat.object)
}
# --------------------------------------------------------------------

# INTEGRATE DATA
# --------------------------------------------------------------------
seurat_integration = function(seurat.object, integration.method.list, normalization.method.list)
{
  set.seed(42)
  print('---- Integration of Seurat object. ----')
  # This will run each integration method for each specified normalization method
  # Each normalization x integration will be stored as a separate reduction result
  for (normalization.method in names(normalization.method.list))
  {
    # Set normalization method to use for IntegrateLayers
    if (normalization.method == 'standard')
    {
      norm.method = 'LogNormalize'
    }

    if (normalization.method == 'sct')
    {
      norm.method = 'SCT'
    }

    norm.assay = normalization.method.list[[normalization.method]]
    pca.name = paste0(normalization.method, '.pca')

    for (integration.method in names(integration.method.list)) {

      print(paste(integration.method, 'integration for', normalization.method, 'normalization.'))

      # Get the actual name of the integration method
      integration.method.name = integration.method.list[[integration.method]]

      # Make the name that will be used for this reduction
      int.reduction.name = paste0(normalization.method, '.', integration.method)

      # If ref_based_integration = 'y' and normalization.method = 'SCT'
      if (ref_based_integration == 'y')
      {
        # For whichever normalization method is used, determine the index position for the sample(s) to be used as references
        if (normalization.method == 'standard')
        {
          # Need to get the position of the reference samples in the RNA assay
          assay.samples = gsub('counts\\.', '', names(seurat.object@assays$RNA@layers)[grep('counts', names(seurat.object@assays$RNA@layers))])
          ref.sample.index = grep(ref_samples, assay.samples)
        }

        if (normalization.method == 'sct')
        {
          assay.samples = names(S@assays$SCT@SCTModel.list)
          ref.sample.index = grep(ref_samples, assay.samples)
        }

        print('Reference-based integration')
        # Run integration (reference-based)
        seurat.object = IntegrateLayers(object = seurat.object, method = integration.method.name,
                                        assay = norm.assay, normalization.method = norm.method,
                                        orig.reduction = pca.name, new.reduction = int.reduction.name,
                                        reference = ref.sample.index)
      }

      if (ref_based_integration == 'n')
      {
        print('Integration')
        # Run integration (normal)
        seurat.object = IntegrateLayers(object = seurat.object, method = integration.method.name,
                                        assay = norm.assay, normalization.method = norm.method,
                                        orig.reduction = pca.name, new.reduction = int.reduction.name)
      }
    }
  }

  # re-join RNA assay layers after integration
  print('Rejoining RNA assay layers...')
  seurat.object[["RNA"]] <- JoinLayers(seurat.object[["RNA"]])

  return(seurat.object)
}
# --------------------------------------------------------------------

# RUN ANNOTATION WITH AZIMUTH
# --------------------------------------------------------------------
run_azimuth_anno = function(seurat.object, azimuth.ref)
{
  print('---- Running Azimuth for Seurat object. ----')
  seurat.object = RunAzimuth(seurat.object, reference = azimuth.ref)
  
  return(seurat.object)
}
# --------------------------------------------------------------------

# LOAD AND VERIFY CUSTOM REFERENCE FOR TRANSFERDATA
# --------------------------------------------------------------------
load_transferdata_ref = function(transferdata.ref.file, transferdata.reduction)
{
  print(paste0('---- Loading custom reference ', transferdata.ref.file, '. ----'))

  transferdata.ref = ''

  # Call load fxn based on extension
  transferdata.ref = switch(
    file_ext(transferdata.ref.file),
    qs = qread(transferdata.ref.file),
    rds = readRDS(transferdata.ref.file),
    stop('The reference file provided for TransferData is not a qs or rds file.')
  )

  # Check that object is Seurat class
  if (class(transferdata.ref)[1] != 'Seurat') {
    stop('Please provide a Seurat object saved as a qs or rds file.')
  }

  if (!transferdata.reduction %in% names(transferdata.ref@reductions)) {
    stop(
      paste0(transferdata.reduction, ' is not a reduction in the provided TransferData reference.\n',
            'Available reductions: ', paste0(names(transferdata.ref@reductions), collapse = ', ')
            )
    )
  }

  return(transferdata.ref)
}
# --------------------------------------------------------------------

# RUN ANNOTATION WITH TRANSFERDATA
# --------------------------------------------------------------------
run_seurat_transferdata = function(seurat.object, transferdata.ref.file, transferdata.reduction, transferdata.annocol)
{
  transferdata.ref = load_transferdata_ref(transferdata.ref.file, transferdata_reduction)

  print('---- Running FindTransferAnchors. ----')
  # Get transfer anchors between query and reference
  transfer.anchors = FindTransferAnchors(
    reference = transferdata.ref,
    query = seurat.object,
    dims = 1:length(transferdata.ref@reductions[[transferdata.reduction]]),
    reference.reduction = transferdata.reduction
  )

  print('---- Running TransferData. ----')
  predictions = TransferData(
    anchorset = transfer.anchors,
    refdata = transferdata.ref@meta.data[, transferdata.annocol])
  
  # colnames(predictions) = paste0('transferdata.', colnames(predictions))

  seurat.object = AddMetaData(
    seurat.object,
    metadata = predictions
    )

  return(seurat.object)
}
# --------------------------------------------------------------------

# GET NUMBER OF PCs QUANTITATIVELY
# --------------------------------------------------------------------
# PCs determined in a quantative fashion, method from:
# https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/elbow_plot_metric.md
get_sigPCs <- function(seurat.object, pca.name)
{
  print('Calculating number of PCs.')

  # Determine percent of variation associated with each PC
  pct <- seurat.object[[pca.name]]@stdev / sum(seurat.object[[pca.name]]@stdev) * 100

  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  df <- cbind(pct, cumu)
  
  # write to file
  print('Writing PC values.')
  write.csv(df, file = paste0(table_dir, '/', project, '_', gsub('\\.', '_', pca.name), '_sigPCs.txt'), quote = FALSE, row.names = FALSE)

  # Determine which PC exhibits cumulative percent greater than 90% and % variation
  # associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]

  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

  pcs <- min(co1, co2)

  print('Writing out final PC value.')
  write.csv(pcs, file = paste0(table_dir, '/', project, '_', gsub('\\.', '_', pca.name), '_sigPC.txt'), quote = FALSE, row.names = FALSE)

  return(pcs)
}
# --------------------------------------------------------------------

# RUN FINDNEIGHBORS, FINDCLUSTERS, UMAP/TSNE
# --------------------------------------------------------------------
find_neighbors_clusters_reductions = function(seurat.object, sig.pcs, reduction.name, resolution.config.list, include.tsne)
{
  # FIND NEIGHBORS
  set.seed(42)
  seurat.object = FindNeighbors(seurat.object, reduction = reduction.name, dims = 1:sig.pcs,
                                graph.name = paste0(reduction.name, c('_nn', '_snn')))

  umap.name = paste0(reduction.name, '.umap')
  umap.key = paste0(gsub('\\.', '', reduction.name), 'UMAP_') # RunUMAP doesn't like '.'

  # RUN UMAP
  print(paste('Running UMAP for', reduction.name))
  seurat.object = RunUMAP(seurat.object, reduction = reduction.name, dims = 1:sig.pcs,
                            reduction.key = umap.key, reduction.name = umap.name)

  if (include.tsne == 'y')
  {
    # RUN TSNE
    print(paste('Running tSNE for', reduction.name))
    seurat.object = RunTSNE(seurat.object, reduction = reduction.name, dims = 1:sig.pcs,
                            reduction.key = paste0(gsub('\\.', '', reduction.name), 'tSNE_'),
                            reduction.name = paste0(reduction.name, '.tsne'))
  }

  # Define the graph.name that is needed for FindClusters
  graph.name = paste0(reduction.name, '_snn')

  # FOR EACH RESOLUTION
  for (res in resolution.config.list)
  {
    print(paste('---- Clustering', graph.name, 'at resolution', res, '. ----'))

    # FIND CLUSTERS BASED ON THE GRAPH AND RES (graph.name_snn_resN.N)
    seurat.object = FindClusters(seurat.object, resolution = res, graph.name = graph.name, verbose = FALSE)

    clust.res = paste0(graph.name, '_res.', res) # if periods need to be removed: gsub('\\.', '_', graph.name)

    print(paste('Plotting', clust.res, 'ON', umap.name))

    pdf(file = paste0(report_path_figures, '/', project, '_', clust.res, '_clusters_umap.pdf'), width = 10, height = 10)
    print(DimPlot(seurat.object, reduction = umap.name, group.by = clust.res, label = TRUE, label.size = 6, repel = TRUE) + NoLegend())
    dev.off()
  }

  # CLUSTREE
  # If there are multiple resolutions in the list, run clustree
  if (length(resolution.config.list) > 1)
  {
    print(paste('---- Running Clustree for all', reduction.name, 'resolution', '(', paste(resolution.config.list, collapse = ', '), '). ----'))

	 title_name = paste0('Normalization + Integration: ', reduction.name)
    clustree.result = clustree(seurat.object, prefix = paste0(graph.name, '_res.'))

    # PLOT: Clustree results
    pdf(file = paste0(report_path_figures, '/clustree/', project, '_', reduction.name, '_clustree_results.pdf'), width = 11, height = 8.5)
    print(clustree.result + ggtitle(title_name) + scale_colour_discrete(name = 'resolution') + scale_size_continuous(name = "cells/cluster", range = c(3, 15)) +
	 guides(colour = guide_legend(order = 1), size = guide_legend(order = 2), alpha = guide_legend(order = 3), fill = guide_legend(order = 4)) +
	 theme(legend.title = element_text(size = 14), legend.text  = element_text(size = 13), legend.key.size = unit(1.2, "lines")))
    dev.off()
  }

  return(seurat.object)
}
# --------------------------------------------------------------------

# RUN FINDNEIGHBORS, FINDCLUSTERS, UMAP/TSNE
# --------------------------------------------------------------------
run_clustering = function(seurat.object, integration.method.list = NULL, normalization.method.list, resolution.config.list, include.tsne)
{
  set.seed(42)
  for (normalization.method in names(normalization.method.list)) {

    if (normalization.method == 'standard')
    {
      pca.name = paste0(normalization.method, '.pca')
    }

    if (normalization.method == 'sct')
    {
      pca.name = paste0(normalization.method, '.pca')
    }

    # Get the sig_pcs here because this will be dependent on the pca reduction being used
    print('---- Calculating number of PCs to use for downstream operations. ----')
    sig.pcs = get_sigPCs(seurat.object, pca.name)

    if (normalization.method == 'sct')
    {
      # need to set sig.pcs to total number of components if +10 to sig.pcs makes it > n_components
      sig.pcs = min(n_components, sig.pcs + 10)
    }

    # Run each integration method
    for (integration.method in names(integration.method.list))
    {
      print(paste('---- Finding Neighbors, Clusters, and Visualization Reductions for',
                  normalization.method, 'normalization  and', integration.method, '. ----'))
      # Make a list of all of the reduction names from the integration (everything but standard.pca and sct.pca)
      reduction.name = paste0(normalization.method, '.', integration.method)
      seurat.object = find_neighbors_clusters_reductions(seurat.object, sig.pcs, reduction.name, resolution.config.list, include.tsne)
    }
  }

  return(seurat.object)
}
# --------------------------------------------------------------------
######################################################################


# SEURAT ANALYSIS
# --------------------------------------------------------------------
# LOAD THE MERGED OBJECT (rds or qs)
S = import_data(initial_seurat_object, 'qs')

# GET NUMBER OF SAMPLES
n_samples = length(unique(S@meta.data[['Sample']]))

# GET LISTS OF CELL CYCLE GENES (convert to mouse IDs if needed)
s_genes <- cc.genes.updated.2019$s.genes
g2m_genes <- cc.genes.updated.2019$g2m.genes
if (tolower(organism) == 'mouse') {
  s_genes <- convert_genes_hs2mm(cc.genes.updated.2019$s.genes)
  g2m_genes <- convert_genes_hs2mm(cc.genes.updated.2019$g2m.genes)
}

# GET A VECTOR OF THE REGRESSION VARIABLES TO USE IN DATA SCALING + ADD USER-SUPPLIED REGRESSION VARS
regression_vars = get_regression_variables(cc_regression, cc_method, mito_regression, ribo_regression, regression_file)

print(regression_vars)

# GET SCORING
S = get_scoring(S, s_genes, g2m_genes, regression_vars, regression_file)

# RUN SEURAT PROCESSING (NORMALIZATION, SCALING, PCA)
S = seurat_processing(S, normalization_method_list, split_layers_by, regression_vars, num_var_features)

# RUN INTEGRATION (IF THERE IS MORE THAN ONE SAMPLE)
if (n_samples > 1) {
  S = seurat_integration(S, integration_method_list, normalization_method_list)
}

# RUN AZIMUTH IF SPECIFIED
if (run_azimuth == 'y') {
  S = run_azimuth_anno(S, azimuth.ref = azimuth_ref)
}

# RUN TRANSFERDATA IF SPECIFIED
if (run_transferdata == 'y') {
  S = run_seurat_transferdata(S, transferdata_ref_file, transferdata_reduction, transferdata_annocol)
}

# FIND NEIGHBORS FOR ALL COMBOS OF INETEGRATION METHODS AND NORMALIZATION METHODS
# FIND CLUSTERS FOR EACH SNN GRAPH AT EACH RESOLUTION
S = run_clustering(S, integration_method_list, normalization_method_list, resolution_config_list, include_tsne)

# SAVE OBJECT
qsave(S, analyzed_seurat_object)
# --------------------------------------------------------------------
