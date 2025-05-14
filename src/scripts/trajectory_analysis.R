# Author: K. Beigel
# Date: 8.18.2024
# Purpose: Trajectory analysis using monocle3.

# DEFINE ARGUMENTS
#--------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)

project <- ''
seurat_object <- ''
assignment_file <- ''
normalization_method <- '' 
integration_method <- '' 
resolution <- ''
partition_trajectory <- ''
final_storage <- ''
lib_path <- ''
provide <- ''
annotate <- ''

# test if there is at least 17 arguments: if not, return an error
nargs = 11
if (length(args) != nargs) 
{
  stop(paste(nargs, 'arguments must be supplied.'), call. = FALSE)
}

if (length(args) == nargs) 
{
	project <- args[1]
	seurat_object <- args[2]
	assignment_file <- args[3]
	normalization_method <- args[4]
	integration_method <- args[5]
	resolution <- args[6]
	partition_trajectory <- args[7]
	final_storage <- args[8]
	lib_path <- args[9]
	provide <- args[10]
	annotate <- args[11]
}

# determine if argument is single value or list
# --------------------------------------------------------------------
single_or_list <- function(variable, comma = ',')
{
   if (grepl(comma, variable, fixed = TRUE))
   {
      variable = as.list(strsplit(variable, comma)[[1]])

      return(variable)
   }

   else
   {
      return(variable)
   }
}
# --------------------------------------------------------------------

# [DR]EFINE VARIABLES
# --------------------------------------------------------------------
final_storage_method = single_or_list(final_storage)
# --------------------------------------------------------------------

# SET CONFIG VAR AS LIST 
# --------------------------------------------------------------------
normalization_method = single_or_list(normalization_method)
integration_method = single_or_list(integration_method)
resolution = as.numeric(single_or_list(resolution))
# --------------------------------------------------------------------

# IMPORT LIBS
# --------------------------------------------------------------------
suppressPackageStartupMessages(library(tidyverse, lib.loc=lib_path))
suppressPackageStartupMessages(library(stringr, lib.loc=lib_path))
suppressPackageStartupMessages(library(xfun, lib.loc=lib_path))
suppressPackageStartupMessages(library(future, lib.loc=lib_path))
suppressPackageStartupMessages(library(Seurat, lib.loc=lib_path))
suppressPackageStartupMessages(library(SeuratWrappers, lib.loc=lib_path))
suppressPackageStartupMessages(library(monocle3, lib.loc=lib_path))
suppressPackageStartupMessages(library(ggplot2, lib.loc=lib_path))
suppressPackageStartupMessages(library(qs, lib.loc=lib_path))
# --------------------------------------------------------------------

# DIRECTORIES
# --------------------------------------------------------------------
base_directory = paste0('data/endpoints/', project, '/analysis/final_analysis/trajectory_analysis')
cds_dir <- paste0(base_directory, '/cds/')
fig_dir <- paste0(base_directory, '/figures/')
dir.create(cds_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
# --------------------------------------------------------------------

# FUNCTIONS
######################################################################
# DEFINE PARTITIONS BASED ON CLUSTER ANNOTATION FILE
# --------------------------------------------------------------------
define_partitions <- function(seurat.object, partition.trajectory, cluster.scheme, assignment.file)
{
	print('Defining partitions.')

	# Need to define the paritions
	partitions <- ''
	
	# If there is a user-supplied partition file, use that info to define partitions
	if (partition.trajectory == 'y')
	{
		# read the table (cluster_annotation_file)
		partitions <- read.table(assignment.file, header = TRUE, sep = '\t') %>%
			mutate_if(is.character, as.factor) %>%
			mutate_if(is.integer, as.factor)

		# Check that there are no clusters ('celltypes') assigned to multiple partitions
		partitions <- partitions %>%
			select(!!sym(cluster.scheme), partition) %>%
			distinct()
		
		if (anyDuplicated(partitions[, cluster.scheme]) != 0) {
			stop(
				paste0("Each 'celltype' must only belong to ONE partition. Please fix '", assignment.file, "' file and re-run.")
			)
		} else {
			print(partitions)
		}
	}

	# If there is not a user-supplied partition file, assign everything to the same partition
	if (partition.trajectory == 'n')
	{	
		# Using the levels of clusters in the chosen cluster scheme, make a table where all clusters are listed and have partition = 1
		partitions <- tibble(
			!!cluster.scheme := levels(seurat.object@meta.data[[cluster.scheme]]), partition = 1) %>%
				mutate_if(is.character, as.factor) %>% mutate_if(is.integer, as.factor)
	}

	# use the partition info to add a new metadata column of the partition designation
	meta.new <- seurat.object@meta.data %>%
		left_join(
			partitions %>% select(celltypes, partition),
			by = join_by(celltypes == celltypes)
		)
	
	seurat.object <- AddMetaData(seurat.object, metadata = meta.new$partition, col.name = 'partition')

	return(seurat.object)
}
# --------------------------------------------------------------------

# CONSTRUCT A cell_data_set OBJECT FOR MONOCLE3 ANALYSIS
# --------------------------------------------------------------------
make_cds <- function(seurat.object, cluster.scheme, norm.method, int.method, redux.umap)
{
	# Make a new reduction named 'pca' from the specified reduction (monocle3 doesn't like custom reduction names)
	seurat.object@reductions$pca <- seurat.object@reductions[[paste0(norm.method, '.pca')]]

	# Make a new reduction named 'umap' from the specified reduction (monocle3 doesn't like custom reduction names)
	seurat.object@reductions$umap <- seurat.object@reductions[[redux.umap]]

	# Make cell_data_set
	cds <- as.cell_data_set(seurat.object)

	# Get partition info from the seurat object
	recreate.partitions <- seurat.object@meta.data$partition

	# Get cluster info from the seurat object
	list.cluster <- seurat.object@meta.data[[cluster.scheme]]

	# Set names and make factors
	#print(cds@colData@rownames)
	names(recreate.partitions) <- cds@colData@rownames
	recreate.partitions <- as.factor(recreate.partitions)

	# Assign the partititons and clusters to the UMAP
	cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
	cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

	# Add the cell embeddings from the seurat object to the cell_data_set
	cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- seurat.object@reductions$umap@cell.embeddings

	return(cds)
}
# --------------------------------------------------------------------

# RUN MONOCLE TRAJECTORY ANALYSIS
# --------------------------------------------------------------------
trajectory_analysis <- function(cds, cluster.scheme, partition.trajectory)
{
	print('Running trajectory analysis.')
	
	if (partition.trajectory == 'y')
	{
		cds <- learn_graph(cds, use_partition = TRUE)
	}

	if (partition.trajectory == 'n')
	{
		cds <- learn_graph(cds, use_partition = FALSE)
	}

	png(filename = paste0(fig_dir, project, '_monocle3_clusters.png'), height = 2000, width = 2700, res = 300)
	# Generate plot
	print(
		plot_cells(cds,
			color_cells_by = cluster.scheme,
			label_groups_by_cluster = FALSE,
			label_branch_points = TRUE,
			label_roots = TRUE,
			label_leaves = TRUE,
			group_label_size = 6,
			graph_label_size = 4
		)
	)
	dev.off()

	png(filename = paste0(fig_dir, project, '_monocle3_partitions.png'), height = 2000, width = 2700, res = 300)
	# Generate plot
	print(
		plot_cells(cds,
			color_cells_by = 'partition',
			label_groups_by_cluster = FALSE,
			label_branch_points = TRUE,
			label_roots = TRUE,
			label_leaves = TRUE,
			group_label_size = 6,
			graph_label_size = 4
		)
	)
	dev.off()

	return(cds)
}
# --------------------------------------------------------------------

# SAVE THE cell_Data_set OBJECT AS qs AND/OR rds
# --------------------------------------------------------------------
save_cds <- function(cds, final.storage)
{
	core_name <- paste0(cds_dir, project, '_monocle3_cds')
	
	for (s in final.storage)
	{
		print(s)

		if (s == 'qs')
		{
			print('You have requested a qs file...OK')
			f_name = paste0(core_name, '.qs')
			qsave(cds, f_name)
		}

		if (s == 'rds')
		{
			print('You have requested a rds file...OK')
			f_name = paste0(core_name, '.RDS')
			saveRDS(cds, file = f_name)
		}
	}
}
# --------------------------------------------------------------------
######################################################################

# GET SEURAT OBJECT EXTENSION AND LOAD SEURAT OBJECT
# --------------------------------------------------------------------
# NEED to FIX THIS -- THIS WILL OPEN BOTH QS AND RDS FILES ERER 12.9.24
print('Importing seurat object.')
S <- ''

extension = file_ext(seurat_object)

if (extension == 'qs')
{
	print('reading in qs file...')
	S = qread(seurat_object)
}

if (extension == 'rds')
{
	print('reading in rds file...')
	S = readRDS(seurat_object)
}
# --------------------------------------------------------------------

# DEFINE ARGUMENTS + VARIABLES
# --------------------------------------------------------------------
# Use the 'celltypes' ident from final_analysis.R
ident <- 'celltypes'

# use final_configs info to define the reduction
redux_umap <- paste0(normalization_method, '.', integration_method, '.umap')

# Determine the assay to be used
assay = ''
if (normalization_method == 'sct') {
	assay = 'SCT'
}

if (normalization_method == 'standard') {
	assay = 'RNA'
}

# Set assay and idents
DefaultAssay(S) <- assay
Idents(object = S) <- ident
# --------------------------------------------------------------------


# MAIN FUNCTION CALLS
# --------------------------------------------------------------------
# Trajectory analysis - define partitions
S = define_partitions(S, partition_trajectory, ident, assignment_file)

# Run Monocle3
cds_obj = make_cds(S, ident, normalization_method, integration_method, redux_umap)
cds_obj = trajectory_analysis(cds_obj, ident, partition_trajectory)
save_cds(cds_obj, final_storage_method)
# --------------------------------------------------------------------
