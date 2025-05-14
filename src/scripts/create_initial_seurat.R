# Author:		K. Beigel
# Date:			9.30.2024

sample_file <- ''
project <- ''
organism <- ''
seurat_creation_source <- ''
run_doubletfinder <- ''
mito_cutoff <- ''
ribo_cutoff <- ''
min_feature_threshold <- ''
max_feature_threshold <- ''
seurat_file_name <- ''
lib_path <- ''

args <- commandArgs(trailingOnly = TRUE)

# test if there is at least 10 arguments: if not, return an error
nargs <- 11
if (length(args) < nargs) {
  stop(paste('At least ', nargs, 'arguments must be supplied.'), call. = FALSE)
} else if (length(args) == nargs) {
	sample_file = args[1]
	project = args[2]
	organism = args[3]
	seurat_creation_source = args[4]
	run_doubletfinder = args[5]
	mito_cutoff = args[6] #mito cut-off for scRNA
	ribo_cutoff = args[7] #mito cut-off for scRNA
	min_feature_threshold = args[8]
	max_feature_threshold = args[9]
	seurat_file_name = args[10]
	lib_path = args[11]
}

library(tidyverse, lib.loc = lib_path)
library(rcartocolor, lib.loc = lib_path)
library(Seurat, lib.loc = lib_path)
library(qs, lib.loc = lib_path)
library(patchwork, lib.loc = lib_path)
library(dplyr, lib.loc = lib_path)

# CREATE DIRECTORIES
#--------------------------------------------------------------------
# Output
base_directory = paste0('data/endpoints/', project, '/analysis/')
figure_dir = paste0(base_directory, '/figures')
rds_dir = paste0(base_directory, '/RDS')
table_dir = paste0(base_directory, '/tables')

for (dir in c(figure_dir, rds_dir, table_dir)) {
	dir.create(dir, showWarnings = FALSE)
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
max_feature_threshold = as.numeric(max_feature_threshold)
#--------------------------------------------------------------------

# IMPORT DATA, CREATE SEURAT OBJECT, ADD METADATA, CALCULATE MITO%
#--------------------------------------------------------------------
get_sample_list = function(sample.file)
{ 
	sample.list = read.table(sample.file, sep = '\t', header = TRUE)
}

# CREATE LIST FOR MERGING
#--------------------------------------------------------------------
make_list_seuobjs = function(sample.list, merge.data.source)
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
make_merged_seuobj = function(seu.obj.list)
{
	# create merged seurat object
	print(paste('Merging samples:', paste(names(seu.obj.list), collapse = ', ')))
	S.merged <- merge(seu.obj.list[[1]], y = c(seu.obj.list[2:length(seu.obj.list)]), add.cell.ids = names(seu.obj.list), project = project)

	return(S.merged)
}
#--------------------------------------------------------------------

# FOR ONE SAMPLE, PULL SINGLE OBJ FROM LIST
#--------------------------------------------------------------------
make_single_seuobj = function(seu.obj.list)
{
	# create merged seurat object
	S.single = seu.obj.list[[1]]
	S.single = RenameCells(object = S.single, add.cell.id = names(seu.obj.list)[1])

	return(S.single)
}
#--------------------------------------------------------------------


# QC PLOTS: PRE- AND POST-FILTERING
#--------------------------------------------------------------------
seu_qc_plots = function(seu.obj, plot.set, title, caption)
{	
	# commented out cols for plots; scheme below is hard to distinguish samples
	# colors from: https://github.com/Nowosad/rcartocolor
	palette = carto_pal(7, "ag_Sunset")
	if (nrow(sample_list) <= 2) {
		palette = palette[c(2, 4)]
	} else if (nrow(sample_list) > 2) {
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

# FILTER SEURAT OBJECT
#--------------------------------------------------------------------
filter_save_seuobj = function(seu.obj)
{
	# pre-filtered plots
	print('Plotting QC figures for data (not filtered).')
	seu_qc_plots(seu.obj, '1', 'Unfiltered', 'Data before filtering')

	# filter
	filter_info = paste0('nFeature_RNA > ', min_feature_threshold, ', nFeature_RNA < ', max_feature_threshold, ', percent.mito < ', mito_cutoff, ', percent.ribo < ', ribo_cutoff)
	print(paste0('Filtering Seurat object: ', filter_info))
	seu.obj.filt <- subset(seu.obj, subset = nFeature_RNA > min_feature_threshold & nFeature_RNA < max_feature_threshold & percent.mito < mito_cutoff & percent.ribo < ribo_cutoff)

	# post-filtered plots
	print('Plotting QC figures for data (filtered).')
	seu_qc_plots(seu.obj.filt, '2',  'Filtered', paste0('Filtered: ', filter_info))

	# save filtered object
	print('Saving Seurat object.')
	qsave(seu.obj.filt, file = seurat_file_name)
}
#--------------------------------------------------------------------

# FUNCTION CALLS
#--------------------------------------------------------------------
sample_list = get_sample_list(sample_file)
seuobjs_list = make_list_seuobjs(sample_list, seurat_creation_source)

seuobj = ''

if (nrow(sample_list) > 1) {
	seuobj = make_merged_seuobj(seuobjs_list)
}

if (nrow(sample_list) == 1) {
	seuobj = make_single_seuobj(seuobjs_list)
}

filter_save_seuobj(seuobj)
#--------------------------------------------------------------------
