#!/usr/bin/env Rscript
# Author:	E. Reichenberger
# Date:		8.08.2024

# Purpose: 	Create cluster images and look at conserved and DE markers

# Updated dec24: added optional visualization outputs (user supplied genes), added totals to 

#(r message=FALSE)

# DEFINE ARGUMENTS
#--------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
print('arguments --------------------')
print(args)
print(length(args))
print('arguments --------------------')

project <- ''
lib_path <- ''
analyzed_seurat_object <- ''
normalization_method <- '' 
integration_method <- '' 
resolution <- ''
celltype_assignment_file <- ''
tsne <- ''
user_gene_file <- ''
processes <- ''
filtering_threshold <- ''
avg2fc <- ''
min.pct <- ''
conserved_genes <- ''
organism <- ''
final_storage <- ''
provide <- ''
annotate <- ''
visualization <- ''
meta_sample <- ''
meta_experiment <- ''
meta_annotation <- ''
umap_reduction <- ''
tnse_reduction <- ''
meomry <- ''

arg_length <- 25

# test for correct argument length
if (length(args) < arg_length) 
{
  statement = paste('At least', arg_length, 'arguments must be supplied.', sep = ' ')
  stop(statement, call.=FALSE)
} 

if (length(args) == arg_length) 
{
	project = args[1]
	lib_path = args[2]
	analyzed_seurat_object = args[3] 
	normalization_method = args[4]
	integration_method = args[5]
	resolution = args[6]
	celltype_assignment_file = args[7]
	tsne = args[8]
	user_gene_file = args[9]
	processes = args[10]
	filtering_threshold = args[11]
	avg2fc = args[12]
	min.pct = args[13]
	conserved_genes = args[14]
	organism = args[15]
	final_storage = args[16]
	provide = args[17]
	annotate = args[18]
	visualization = args[19]
	meta_sample = args[20]
	meta_experiment = args[21]
	meta_annotation = args[22]
	umap_reduction = args[23]
	tnse_reduction = args[24]
	memory = args[25]
}
# --------------------------------------------------------------------

# DEFINE VARIABLES
filtering_threshold <- as.numeric(filtering_threshold)
avg2fc <- as.numeric(avg2fc)
processes <- as.integer(processes)
min.pct <- as.numeric(min.pct)
memory <- as.numeric(memory)
# --------------------------------------------------------------------

# determine if argument is single value or list
# --------------------------------------------------------------------
single_or_list <- function(variable, comma = ',')
{
   if ( grepl( comma, variable, fixed = TRUE) )
   {
      variable = as.list(strsplit(variable, comma)[[1]])

      return(variable)
   }

   else
   {
      return(variable)
   }/loup
}
# --------------------------------------------------------------------

# SET CONFIG VAR AS LIST 
# --------------------------------------------------------------------
resolution = as.numeric(single_or_list(resolution))
normalization_method = single_or_list(normalization_method)
integration_method = single_or_list(integration_method)
# visualization method defined towards end of script

final_storage_method = ''
if (is.null(final_storage) == TRUE)
{
	final_storage_method = as.list('qs')
}

if (is.null(final_storage) == FALSE)
{
	final_storage_method = single_or_list(final_storage)
}
# --------------------------------------------------------------------

# AUTOACCEPT THE LOUPER EULA
# --------------------------------------------------------------------
if ('cloupe' %in% unlist(final_storage_method)) 
{
	options(
		repos = c(CRAN = "internalrepo"),
		download.file.method = "curl",
		download.file.extra = "-k -L")

	# Load library	
	library(loupeR, lib.loc=lib_path))

	# Set the EULA auto accept env var to "y"
	Sys.setenv(AUTO_ACCEPT_EULA = "y")
	# AUTO_ACCEPT_ENV_VAR = 'y'

	# Assign env var to another env var to comply with:
	# https://github.com/10XGenomics/loupeR/blob/d7994c5c5bd1fc7a6202e23f84ecf9df43831e6e/R/eula.R
	AUTO_ACCEPT_ENV_VAR = "AUTO_ACCEPT_EULA"

	# Run the loupeR:::eula() function which will check AUTO_ACCEPT_ENV_VAR with auto_accepted_eula()
	loupeR:::eula()
	loupeR::setup()
}
# --------------------------------------------------------------------

# IMPORT LIBS
#--------------------------------------------------------------------
suppressPackageStartupMessages(library(Seurat, lib.loc=lib_path))
suppressPackageStartupMessages(library(RColorBrewer, lib.loc=lib_path))
suppressPackageStartupMessages(library(dplyr, lib.loc=lib_path))
suppressPackageStartupMessages(library(ggrepel, lib.loc=lib_path))
suppressPackageStartupMessages(library(stringr, lib.loc=lib_path))
suppressPackageStartupMessages(library(sctransform, lib.loc=lib_path))
suppressPackageStartupMessages(library(ggplot2, lib.loc=lib_path))
suppressPackageStartupMessages(library(reshape2, lib.loc=lib_path))
suppressPackageStartupMessages(library(data.table, lib.loc=lib_path))
suppressPackageStartupMessages(library(qs, lib.loc=lib_path))
suppressPackageStartupMessages(library(future, lib.loc=lib_path))
suppressPackageStartupMessages(library(progressr, lib.loc=lib_path))
suppressPackageStartupMessages(library(presto, lib.loc=lib_path))
suppressPackageStartupMessages(library(fgsea, lib.loc=lib_path))
suppressPackageStartupMessages(library(msigdbr, lib.loc=lib_path))
suppressPackageStartupMessages(library(viridis, lib.loc=lib_path))
suppressPackageStartupMessages(library(purrr, lib.loc=lib_path))
suppressPackageStartupMessages(library(readr, lib.loc=lib_path))
suppressPackageStartupMessages(library(tidyverse, lib.loc=lib_path))
suppressPackageStartupMessages(library(tools, lib.loc=lib_path))
#--------------------------------------------------------------------

# PARALLEL w/ FUTURE + SET SEED
#--------------------------------------------------------------------
print(paste0('Assigning ', memory, ' bytes of memory...'))
options(future.globals.maxSize = memory)
#options(future.globals.maxSize = 10001 * 1024^2) 
plan(multisession(workers = as.integer(processes)))
set.seed(42)
#--------------------------------------------------------------------

# COLOR SCHEME
#--------------------------------------------------------------------
col=colorRampPalette(c("#FDE725FF","#29AF7FFF",'#39558CFF', "#440154FF"))(100)
#--------------------------------------------------------------------

# DIRECTORIES
#--------------------------------------------------------------------
base_directory = paste('data/endpoints/', project, '/analysis/final_analysis', sep='')

# -----------------------------------------------------
# want to delete any output that was previously created.
if (dir.exists(base_directory)) {
	  unlink(base_directory, recursive = TRUE)  
}

# Create a new, empty directory
dir.create(base_directory, showWarnings=FALSE)
# -----------------------------------------------------

rds_dir <- paste('data/endpoints/', project, '/analysis/RDS/', sep='')
fig_dir <- paste0(base_directory, '/figures/')
tbl_dir <- paste0(base_directory, '/tables/')
pathway_dir <- paste0(base_directory, '/pathway_analysis/')
pathway_dir_table <- paste0(base_directory, '/pathway_analysis/tables/')
pathway_dir_figure <- paste0(base_directory, '/pathway_analysis/figures/')

dir.create(fig_dir, showWarnings=FALSE)
dir.create(tbl_dir, showWarnings=FALSE)
dir.create(pathway_dir, showWarnings=FALSE)
dir.create(pathway_dir_table, showWarnings=FALSE)
dir.create(pathway_dir_figure, showWarnings=FALSE)
#--------------------------------------------------------------------

# FUNCTIONS ---------------------------------------------------------

# CHANGE META COLUMN NAMES IN USER-SUPPLIED S.O.
#--------------------------------------------------------------------
providing_object <- function(provide, annotate, seurat_object, meta_sample, meta_experiment, meta_annotation)
{
	if (provide == 'y')
	{
		print('A user-supplied seurat object has been indicated ')
		print('Converting meta names in seurat object...')
		meta_length = length(colnames(seurat_object@meta.data))

		for (m in 1:length(colnames(seurat_object@meta.data)))
		{
			if (colnames(seurat_object@meta.data)[m] == meta_sample)
			{
				colnames(seurat_object@meta.data)[m] <- 'Sample'
			}
			
			if (colnames(seurat_object@meta.data)[m] == meta_experiment)
			{
				colnames(seurat_object@meta.data)[m] <- 'Experiment'
			}
			
			if (annotate == 'n')
			{
				if (colnames(seurat_object@meta.data)[m] == meta_annotation)
				{
					colnames(seurat_object@meta.data)[m] <- 'celltypes'
				}
			}
		}
	}

	return(seurat_object)
}
# ------------------------------------------------

# RENAME CLUSTERS + PROPORTIONS + DIMPLOT
#--------------------------------------------------------------------
rename_and_visualize <- function(seurat_object, celltype_file, ident, genes, markers, umap, tsneR, tsne_yn, f_dir, t_dir, rename_flag, visi, project, normalization_method)
{
	if (rename_flag == 0)
	{
		print('renaming clusters...')

		annotation <- read.table(celltype_file, header=TRUE, sep='\t')

		old.cluster.ids <- annotation$cluster
		new.cluster.ids <- annotation$celltypes
		print(old.cluster.ids)
		print(new.cluster.ids)

		#print(levels(seurat_object))
		 #[1] "0"  "1"  "10" "11" "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9" 
		# levels must human-readable ordered
		levels(seurat_object) <- seq(0, length(annotation$celltypes)-1)

		names(new.cluster.ids) <- levels(seurat_object)
		seurat_object <- RenameIdents(seurat_object, new.cluster.ids)

		# create new meta.data column and set equal to active ident
		seurat_object@meta.data$celltypes <- Idents(seurat_object) # original
	}

	print('Making Heatmap of Variable Genes...')
	variable_genes_file <- paste0('data/endpoints/', project, '/analysis/tables/', project, '_', normalization_method, '_variable_genes.txt')
	variable_genes <- read.table(variable_genes_file, header=TRUE)
	variable_genes <- variable_genes$x

	heat_file <- paste0(f_dir, project, '_', normalization_method, '_variable_genes_heatmap.png')
	png(filename=heat_file, width=2700,height=2200,res=300)
	print(DoHeatmap(subset(seurat_object, downsample = 100), features = variable_genes[1:75], size = 3) + 
		theme(axis.text.y = element_text(size = 7)) + 
		theme(legend.position = 'none'))
	dev.off()

	print('Calculating Z-scores for each cluster...')
	#-----------------------------------------------------------
	ae <- AggregateExpression(object = seurat_object, group.by = 'celltypes')$RNA
	z_ae <- as.data.frame(scale(ae))
	z_ae <- tibble::rownames_to_column(z_ae, "gene")
	z_scores_report = paste(t_dir, project, '_z_scores.', ident, '_celltypes.txt', sep='')
	write.table(z_ae, file=z_scores_report, sep='\t', quote=F, col.names=TRUE, row.names=FALSE)

	print('Calculating Z-scores for each cluster AND experiment...')
	#-----------------------------------------------------------
	ae2 <- AggregateExpression(object = seurat_object, group.by = c('celltypes', 'Experiment'))$RNA
	z_ae2 <- as.data.frame(scale(ae2))
	z_ae2 <- tibble::rownames_to_column(z_ae2, "gene")
	z_scores_report2 = paste(t_dir, project, '_z_scores.', ident, '_celltypes_experiment.txt', sep='')
	write.table(z_ae2, file=z_scores_report2, sep='\t', quote=F, col.names=TRUE, row.names=FALSE)
	#-----------------------------------------------------------

	print('Calculating Z-scores for each cluster AND experiment AND sample...')
	#-----------------------------------------------------------
	ae3 <- AggregateExpression(object = seurat_object, group.by = c('celltypes', 'Experiment', 'Sample'))$RNA
	z_ae3 <- as.data.frame(scale(ae3))
	z_ae3 <- tibble::rownames_to_column(z_ae3, "gene")
	z_scores_report3 = paste(t_dir, project, '_z_scores.', ident, '_celltypes_experiment_sample.txt', sep='')
	write.table(z_ae3, file=z_scores_report3, sep='\t', quote=F, col.names=TRUE, row.names=FALSE)
	#-----------------------------------------------------------

	print('Calculating cell counts/proportions...')
	print('...by experiment...')
	#-----------------------------------------------------------

	# add in totals (w/ + without %)  #added 12.4.2024 ERR
	number_perCluster_experiment <- table(seurat_object@meta.data$Experiment, seurat_object@meta.data[['celltypes']])
	number_perCluster_experiment_prop <- round(proportions(as.matrix(number_perCluster_experiment), 1), 3)

	# calculate column sums (a), calculate sum  of column sum (b), append b to a
	exp_column_sums <- colSums(number_perCluster_experiment)
	exp_column_sums <- append(exp_column_sums, sum(exp_column_sums))

	# calculate row sums for number_perCluster_experiment 
	exp_row_sums <- rowSums(number_perCluster_experiment)

	# merge number_perCluster_experiment + number_perCluster_experiment_prop
	Xe = as.matrix(number_perCluster_experiment)
	Ye = as.matrix(number_perCluster_experiment_prop)
	Ze <- as.data.frame(matrix(paste0(Xe, " (", Ye, ")"), nrow = nrow(Xe), dimnames = dimnames(Xe)))

	# add row sums (exp_row_sums) to Ze
	Ze$totals <- exp_row_sums

	# add column sums (exp_column_sums) to Ze
	Ze_bound <- rbind(Ze, exp_column_sums)

	# add `totals` as final rowname
	rownames(Ze_bound)[rownames(Ze_bound) == nrow(Ze_bound) ] <- 'totals'
	
	# write counts/props/totals to file
	cluster_counts_proportions_experiment <- paste0(t_dir, project, '_', ident, '_final_cluster_counts_proportions_experiment.txt')
	write.table(Ze_bound, file=cluster_counts_proportions_experiment, sep='\t', quote=FALSE)

	print('...making barplot...')
	df_e <- reshape2::melt(number_perCluster_experiment)
	names(df_e)[1] <- 'Experiment'
	names(df_e)[2] <- 'Cluster'

	df_e$Cluster <- as.factor(df_e$Cluster)

	bplot_experiment_fname <- paste0(f_dir, project, '_barplot_proportions_experiment_', ident, '.png')
	png(filename=bplot_experiment_fname, width=1600, height=2000, res=300)
	print(
		ggplot(df_e, aes(x=Experiment, y=value, fill=Cluster)) +
			geom_bar(stat='identity', position='fill', width=0.8)  +
			ylab('Cluster Proportion') +
			xlab('Experiment') +
			guides(fill=guide_legend(title='Cluster')) +
			theme_bw() +
			theme(
				axis.text = element_text(size=12),
				axis.title = element_text(size=14),
				axis.text.x = element_text(angle=45, hjust=1, vjust=1)
			)
	)
	dev.off()

	print('...by sample...')
	number_perCluster_sample <- table(seurat_object@meta.data$Sample, seurat_object@meta.data[['celltypes']])
	number_perCluster_sample_prop <- round(proportions(as.matrix(number_perCluster_sample), 1), 3) # nolint: line_length_linter.

	# calculate column sums (a), calculate sum  of column sum (b), append b to a
	sample_column_sums <- colSums(number_perCluster_sample)
	sample_column_sums <- append(sample_column_sums, sum(sample_column_sums))

	# calculate row sums for number_perCluster_sample
	sample_row_sums <- rowSums(number_perCluster_sample)

	# merge number_perCluster_sample + number_perCluster_sample_prop
	Xs = as.matrix(number_perCluster_sample)
	Ys = as.matrix(number_perCluster_sample_prop)
	Zs <- as.data.frame(matrix(paste0(Xs, " (", Ys, ")"), nrow = nrow(Xs), dimnames = dimnames(Xs)))

	# add row sums (sample_row_sums) to Ze
	Zs$totals <- sample_row_sums

	# add column sums (sample_column_sums) to Ze
	Zs_bound <- rbind(Zs, sample_column_sums)

	# add `totals` as final rowname
	rownames(Zs_bound)[rownames(Zs_bound) == nrow(Zs_bound) ] <- 'totals'

	# write counts/props/totals to file
	cluster_counts_proportions_sample <- paste0(t_dir, project, '_', ident, '_final_cluster_counts_proportions_sample.txt')
	write.table(Zs_bound, file=cluster_counts_proportions_sample, sep='\t', quote=FALSE)

	print('...making barplot...')
	df_s <- reshape2::melt(number_perCluster_sample)
	names(df_s)[1] <- 'Sample'
	names(df_s)[2] <- 'Cluster'

	df_s$Cluster <- as.factor(df_s$Cluster)

	bplot_sample_fname <- paste0(f_dir, project, '_barplot_proportions_sample_', ident, '.png')
	png(filename=bplot_sample_fname, width=1600, height=2000, res=300)
	print(
		ggplot(df_s, aes(x=Sample, y=value, fill=Cluster)) +
			geom_bar(stat='identity', position='fill', width=0.8)  +
			ylab('Cluster Proportion') +
			xlab('Sample') +
			guides(fill=guide_legend(title='Cluster')) +
			theme_bw() +
			theme(
				axis.text = element_text(size=12),
				axis.title = element_text(size=14),
				axis.text.x = element_text(angle=45, hjust=1, vjust=1)
			)
	)
	dev.off()
	#-----------------------------------------------------------

	# phase information (added sept 2025) ----------------------

	# Experiment
	print('making table of phases across all clusters (by experiment)...')
	cc_phase_cluster <- as.data.frame(table(seurat_object@meta.data$Experiment, seurat_object@meta.data[['celltypes']], seurat_object@meta.data$Phase))
	colnames(cc_phase_cluster) <- c('Experiment', 'Cluster', 'Phase', 'Frequency')

	cc_phase_cluster_prop <- as.data.frame(prop.table(table(seurat_object@meta.data$Experiment, seurat_object@meta.data[['celltypes']], seurat_object@meta.data$Phase), margin = 1) * 100)
	colnames(cc_phase_cluster_prop) <- c('Experiment', 'Cluster', 'Phase', 'Frequency')

	cc_phase_cluster$Combined_Frequency <- paste0(cc_phase_cluster$Frequency, ' (', round(cc_phase_cluster_prop$Frequency, 1), '%)')

	# write table
	cc_phase_cluster <- cc_phase_cluster %>% arrange(Cluster, Phase)
	phase_table = paste(t_dir, project, '_', ident, '_phase_experiment.txt', sep='')
	write.table(cc_phase_cluster, file=phase_table, sep='\t', quote=F, col.names=TRUE, row.names=FALSE)

	# Sample
	print('making table of phases across all clusters (by sample)...')
	cc_phase_cluster_sample <- as.data.frame(table(seurat_object@meta.data$Sample, seurat_object@meta.data[['celltypes']], seurat_object@meta.data$Phase))
	colnames(cc_phase_cluster_sample) <- c('Sample', 'Cluster', 'Phase', 'Frequency')

	cc_phase_cluster_sample_prop <- as.data.frame(prop.table(table(seurat_object@meta.data$Sample, seurat_object@meta.data[['celltypes']], seurat_object@meta.data$Phase), margin=1) * 100)
	colnames(cc_phase_cluster_sample_prop) <- c('Sample', 'Cluster', 'Phase', 'Frequency')

	cc_phase_cluster_sample$Combined_Frequency <- paste0(cc_phase_cluster_sample$Frequency, ' (', round(cc_phase_cluster_sample_prop$Frequency, 1), '%)')

	# write
	phase_table_sample = paste(t_dir, project, '_', ident, '_phase_sample.txt', sep='')
	write.table(cc_phase_cluster_sample, file=phase_table_sample, sep='\t', quote=F, col.names=TRUE, row.names=FALSE)
	# ---------------------------------------------------------


	# dimplots -------------------------------------------------
	print('creating dimplot images...')
	#-----------------------------------------------------------
	fname1 <- paste0(f_dir, project, '_', ident, '_umap.png')
	png(filename=fname1, width=2700,height=2000,res=300)
	print(DimPlot(seurat_object, reduction = umap, label=TRUE, group.by="celltypes") + NoLegend()) #showing clusters with new names
	dev.off()

	fname2 <- paste0(f_dir, project, '_', ident, '_umap_experiment.png')
	png(filename=fname2, width=2400, height=2000, res=300)
	print(
		DimPlot(seurat_object, reduction = umap, label=TRUE, split.by='Experiment') + 
		NoLegend() + 
		theme(
			panel.spacing=unit(2, "lines"), # increase space between samples
			strip.text.x.top=element_text(size=20) # increase the labels for the split.by var
		)
	)
	dev.off()

	fname3 <- paste0(f_dir, project, '_', ident, '_umap_sample.png')
	png(filename=fname3, width=2400, height=2000, res=300)
	print(
		DimPlot(
			seurat_object,
			reduction = umap,
			label=FALSE,
			split.by='Sample'
		) + 
		theme(
			strip.text.x.top=element_text(size=16, angle=90), # increase the labels for the split.by var
			axis.text.x = element_text(size=8)
		)
	)
	dev.off()

	fname4 <- paste0(f_dir, project, '_', ident, '_umap_phase.png')
	png(filename=fname4, width=2700,height=2000,res=300)
	print(DimPlot(seurat_object, reduction = umap, label=FALSE, group.by="Phase") + labs(title = 'Cell Cycle Phase')) #showing clusters with new names
	dev.off()

	# added additional phase plots sept 2025
	fname4e <- paste0(f_dir, project, '_', ident, '_umap_phase_experiment.png')
	png(filename=fname4e, width=2700,height=2000,res=300)
	print(DimPlot(seurat_object, reduction = umap, label=FALSE, group.by="Phase", split.by='Experiment') + labs(title = 'Cell Cycle Phase (split by experiment)')) #showing clusters with new names
	dev.off()

	fname4s <- paste0(f_dir, project, '_', ident, '_umap_phase_sample.png')
	png(filename=fname4s, width=2700,height=2000,res=300)
	print(DimPlot(seurat_object, reduction = umap, label=FALSE, group.by="Phase", split.by='Sample') + 
		labs(title = 'Cell Cycle Phase (split by sample)') +  #showing clusters with new names
		theme(
			strip.text.x.top=element_text(size=16, angle=90), # increase the labels for the split.by var
			axis.text.x = element_text(size=8)
		))
	dev.off()

	fname4s_b <- paste0(f_dir, project, '_', ident, '_geom_point_phase_sample.png')
	png(filename=fname4s_b, width=2700,height=2000,res=300)
	print(ggplot(cc_phase_cluster_sample, aes(x=Cluster, y=Frequency, shape=Phase, color=Sample)) + 
		geom_point(size=3) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		labs(x='Cluster', y='Frequency', shape='Cell Cycle Phase', color='Sample'))
	dev.off()

	fname4e_b <- paste0(f_dir, project, '_', ident, '_geom_point_phase_experiment.png')
	png(filename=fname4e_b, width=2700,height=2000,res=300)
	print(ggplot(cc_phase_cluster, aes(x=Cluster, y=Frequency, shape=Phase, color=Experiment)) + 
		geom_point(size=3) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		labs(x='Cluster', y='Frequency', shape='Cell Cycle Phase', color='Experiment'))
	dev.off()

	f1 = paste(f_dir, project, '_', ident, '_final_cluster_plots.pdf', sep='')
	pdf(file=f1, onefile=TRUE, width=11, height=8.5)

	print(DimPlot(seurat_object, reduction = umap, label=TRUE) + NoLegend() ) #showing clusters with new names

	print(ggplot(df_e, aes(x=Experiment, y=value, fill=Cluster)) +
	geom_bar(stat='identity', position='fill')  +
	ylab('Cluster Proportion') + xlab('Experiment') + guides(fill=guide_legend(title='Cluster')))

	print(ggplot(df_s, aes(x=Sample, y=value, fill=Cluster)) +
	geom_bar(stat='identity', position='fill')  +
	ylab('Cluster Proportion') + xlab('Sample') + guides(fill=guide_legend(title='Cluster')))

	if (tsne_yn == 'y')
	{
		print(DimPlot(seurat_object, reduction = tsneR, label=TRUE)) #showing clusters with new names
		print(DimPlot(seurat_object, reduction = tsneR, label=TRUE, split.by='Experiment')) #showing clusters with new names by experiment
		print(DimPlot(seurat_object, reduction = tsneR, label=TRUE, split.by='Sample')) #showing clusters with new names by sample
	}
	#-----------------------------------------------------------

	# Visualize user-defined genes ----------------------------
	if (is.null(genes) == FALSE)
	{
		# make new meta column (to have expression legend in dotplots)
		seurat_object@meta.data[["experiment_split"]] <- paste0(Idents(seurat_object), "_", seurat_object@meta.data[["Experiment"]]) 
		seurat_object@meta.data[["sample_split"]] <- paste0(Idents(seurat_object), "_", seurat_object@meta.data[["Sample"]], "_", seurat_object@meta.data[["Experiment"]]) 

		temp_vec <- ''

		# works well for 12 or fewer genes.
		# split big list into smaller bites, visualize bites
		if (length(genes) > 12)
		{
			print('Spliting gene list into smaller, manageable bites....')

			i <- 1
			j <- 1
			bin_count <- ceiling(length(genes)/12)
			temp_vec <- vector("list", bin_count)

			while( i <= length(markers$V1) )
			{
				if( length(temp_vec[[j]]) < 12 )
				{
					temp_vec[[j]] <- c(temp_vec[[j]], markers$V1[i])
					i <- i+1
				}

				else
				{
					j <- j+1
				}
			}

			for (t in temp_vec)
			{
				print('Making plots highlighting user-provided genes...')
				print(t)
				#try(print(FeaturePlot(seurat_object, label=TRUE, repel=TRUE, label.size=2.5, reduction=umap, features=c(t))))

				if (!is.null(visi) == TRUE)
				{
					print('You have opted to make some additional visualalization plots.')
					for (v in visi)
					{
						if (v == 'feature')
						{
							print('making feature plots...')
							try(print(FeaturePlot(seurat_object, label=TRUE, repel=TRUE, label.size=2.5, reduction=umap, features=c(t))))
						}

						if (v == 'dot')
						{
							print('making dot plots...')
	                      	try(print(DotPlot(seurat_object, features=c(t), group.by='experiment_split') + xlab('Genes') + ylab('Clusters') + RotatedAxis()))
 	                      	try(print(DotPlot(seurat_object, features=c(t), group.by='sample_split') + xlab('Genes') + ylab('Clusters') + RotatedAxis()))
						}

						if (v == 'violin')
						{
							print('making violin plots...')
							try(plot(VlnPlot(seurat_object, features=c(t), split.by='Experiment') + xlab(label = 'Clusters') + ylab('LogNormalized Expression')))
						}

						if (v == 'ridge')
						{
							print('making ridge plots...')
							try(print(RidgePlot(seurat_object, features=t, stack=TRUE) + ylab(label='Clusters')))
						}
					}
				}
			}
		}

		else
		{
			print('Making plots highlighting user-provided genes...')
			print(genes)

			if (!is.null(visi) == TRUE)
			{
				for (v in visi)
				{
					if (v == 'feature')
					{
						print('making feature plots...')
						try(print(FeaturePlot(seurat_object, label=TRUE, repel=TRUE, label.size=2.5, reduction=umap, features=c(genes))))
					}

					if (v == 'dot')
					{
						print('making dot plots...')
						try(print(DotPlot(seurat_object, features=c(genes), group.by='experiment_split') + xlab('Genes') + ylab('Clusters') + RotatedAxis()))
						try(print(DotPlot(seurat_object, features=c(genes), group.by='sample_split') + xlab('Genes') + ylab('Clusters') + RotatedAxis()))
					}

					if (v == 'violin')
					{
						print('making violin plots...')
						try(plot(VlnPlot(seurat_object, features=c(genes), split.by='Experiment') + xlab(label = 'Clusters') + ylab('LogNormalized Expression')))
					}

					if (v == 'ridge')
					{
						print('making ridge plots...')
						try(print(RidgePlot(seurat_object, features=genes, stack=TRUE) + ylab(label='Clusters')))
					}
				}
			}
		}
	}

	dev.off() # dimplot + proportions + visualize user-provided genes
	# ---------------------------------------------------------

	return(seurat_object)
}
#--------------------------------------------------------------------

# CREATE PATHWAY IMAGES --------------------------------------------------------
bioimages <- function(lollipop_data, f_header, ident, pathway_dir, cluster, project, rank, comparison)
{
	print('making pathway images...')
	fname = paste0(pathway_dir_figure, project, '_', ident, '_bio_processess_cluster_', f_header, '_', comparison, '_', rank, '.png')
	fname_pdf = paste0(pathway_dir_figure, project, '_', ident, '_bio_processess_cluster_', f_header, '_', comparison, '_', rank, '.png')
	print(fname)
	p_size = 8

	if (length(lollipop_data$NES) <= 15)
	{
		print('less than 15')
		print(c('number of bio processess:', length(lollipop_data$NES)))
		png(filename=fname, width=2700,height=2000,res=300)
		print(fname)
	}

	if (length(lollipop_data$NES) > 15 & length(lollipop_data$NES) < 40)
	{
		print('more than 15, less than 40')
		print(c('number of bio processess:', length(lollipop_data$NES)))
		png(filename=fname, width=5700,height=3750,res=300)
		p_size = 14
	}

	if (length(lollipop_data$NES) > 40 & length(lollipop_data$NES) < 75)
	{
		print('more than 40, less than 75')
		print(c('number of bio processess:', length(lollipop_data$NES)))
		png(filename=fname, width=8000,height=5000,res=300)
		p_size = 20
	}

	if (length(lollipop_data$NES) > 75)
	{
		print('more than 75')
		print(c('number of bio processess:', length(lollipop_data$NES)))
		png(filename=fname, width=4000,height=7500,res=300)
	}

	if (length(lollipop_data$NES) > 0)
	{
		print('plotting...')
		print(ggplot(lollipop_data, aes(x = Description, y = NES)) +
		geom_segment( aes(x = Description, xend = Description, y = 0, yend = NES, color = p.adjust)) +
		geom_point(aes(color = p.adjust, size = abs(NES))) +
		scale_color_viridis(limits = c(0, 0.1)) +
		coord_flip() +
		theme_light() +
		scale_x_discrete(limits = rev, labels = function(x) str_wrap(x, width = 35)) +
		geom_hline(yintercept=0, color = "black") +
		theme(
			panel.border = element_rect(colour = "black"),
			axis.ticks = element_line(color = "black"),
			axis.title.y = element_blank(),
			legend.key.height = unit(1, 'cm')
			) + 
		theme(strip.background =element_rect(fill="grey90", color = "black"),
		strip.text = element_text(colour = 'black', face = "bold")) +
		theme(text = element_text(size=p_size)) +
		facet_grid(database ~ ., scales = "free_y", space = "free_y") + 
		labs(color = "Adjusted p-value", size = "Absolute Value\nNormalized\nEnrichment\nScore", y = "Normalized Enrichment Score"))
	}

	dev.off()
	print('finished plotting image')
	print('')
}
#--------------------------------------------------------------------

# https://github.com/ctlab/fgsea/issues/151#issuecomment-2088857387
# this will need to be corrected with this issue is addressed by fgsea authors
#fgsea(pathways=dbr_df, stats = rankData+rnorm(length(rankData), sd=0.001), minSize = 10, maxSize = 500) 

f_gsea <- function(df, gmt_name, rank, filtering_threshold, pathway_dir, ident, curated_gene_sets, hallmarks, ontology_gene_sets, project, comp)
{
	print('looking for pathways....')
	set.seed(42)
	print(length(df$gene))
	print(rank)

	rankData = ''

	if (rank == 'avg_log2fc')
	{
		print('ranking by avg logFC')
		df <- df[order(df$avg_log2FC, decreasing = TRUE), ] # want highest to lowest
		rankData <- df$avg_log2FC
	}

	if (rank == 'p_val_adj')
	{
		print('ranking by adj p.value')
		df <- df[order(df$p_val_adj, decreasing = FALSE), ] # want from lowest to highest
		rankData <- df$p_val_adj
	}

	names(rankData) <- df$gene
	cluster = unique(df$cluster)
	print(cluster)

	#curated
	msigdbr_list_curated <- split(x=curated_gene_sets$gene_symbol, f=curated_gene_sets$gs_name)
	#gsea_output_curated <- fgsea(pathways=msigdbr_list_curated, rankData)
	gsea_output_curated <- fgsea(pathways=msigdbr_list_curated, stats = rankData+rnorm(length(rankData), sd=0.001), minSize = 10, maxSize = 500) 
	gsea_output_curated_filtered <- subset(gsea_output_curated, padj < filtering_threshold)

	colnames(gsea_output_curated_filtered)[1] <- 'Description'
	colnames(gsea_output_curated_filtered)[3] <- 'p.adjust'
	gsea_output_curated_filtered <- gsea_output_curated_filtered[, c('Description', 'p.adjust', 'ES', 'NES', 'size', 'leadingEdge')]
	gsea_output_curated_filtered$p.adjust <- round(gsea_output_curated_filtered$p.adjust, 3)
	gsea_output_curated_filtered$ES <- round(gsea_output_curated_filtered$ES, 3)
	gsea_output_curated_filtered$NES <- round(gsea_output_curated_filtered$NES, 3)
	gmt_curated <- paste0(gmt_name, '_curated')
	gsea_output_curated_filtered$cluster <- cluster
	gsea_output_curated_filtered$database <- 'curated'
	gsea_output_curated_filtered$comparison <- comp

	if (length(gsea_output_curated_filtered$NES) > 0)
	{
		bioimages(gsea_output_curated_filtered, gmt_curated, ident, pathway_dir, cluster, project, rank, comp)
		fname1=paste0(pathway_dir_table, project, '_', ident, '_', 'cluster_', gmt_curated, '_filtered_results_', rank, '.csv', sep='')
		print(fwrite(gsea_output_curated_filtered, file=fname1))
		print('finished curated')
	}

	#hallmark
	msigdbr_list_hallmark <- split(x=hallmarks$gene_symbol, f=hallmarks$gs_name)
	set.seed(42)
	#gsea_output_hallmark <- fgsea(pathways=msigdbr_list_hallmark, rankData)
	gsea_output_hallmark <- fgsea(pathways=msigdbr_list_hallmark, stats = rankData+rnorm(length(rankData), sd=0.001), minSize = 10, maxSize = 500) 
	gsea_output_hallmark_filtered <- subset(gsea_output_hallmark, padj < filtering_threshold)

	colnames(gsea_output_hallmark_filtered)[1] <- 'Description'
	colnames(gsea_output_hallmark_filtered)[3] <- 'p.adjust'
	gsea_output_hallmark_filtered <- gsea_output_hallmark_filtered[, c('Description', 'p.adjust', 'ES', 'NES', 'size', 'leadingEdge')]
	gsea_output_hallmark_filtered$p.adjust <- round(gsea_output_hallmark_filtered$p.adjust, 3)
	gsea_output_hallmark_filtered$ES <- round(gsea_output_hallmark_filtered$ES, 3)
	gsea_output_hallmark_filtered$NES <- round(gsea_output_hallmark_filtered$NES, 3)
	gmt_hallmark <- paste0(gmt_name, '_hallmark')
	gsea_output_hallmark_filtered$cluster <- cluster
	gsea_output_hallmark_filtered$database <- 'hallmark'
	gsea_output_hallmark_filtered$comparison <- comp

	if (length(gsea_output_hallmark_filtered$NES) > 0)
	{
		bioimages(gsea_output_hallmark_filtered, gmt_hallmark, ident, pathway_dir_figure, cluster, project, rank, comp)
		fname2=paste0(pathway_dir_table, project, '_', ident, '_', 'cluster_', gmt_hallmark, '_filtered_results_', rank, '.csv', sep='')
		print(fwrite(gsea_output_hallmark_filtered, file=fname2))
		print('finished hallmark')
	}

	#ontology
	msigdbr_list_ontology <- split(x=ontology_gene_sets$gene_symbol, f=ontology_gene_sets$gs_name)
	set.seed(42)
	#gsea_output_ontology <- fgsea(pathways=msigdbr_list_ontology, rankData) 
	gsea_output_ontology <- fgsea(pathways=msigdbr_list_ontology, stats = rankData+rnorm(length(rankData), sd=0.001), minSize = 10, maxSize = 500) 
	gsea_output_ontology_filtered <- subset(gsea_output_ontology, padj < filtering_threshold)

	colnames(gsea_output_ontology_filtered)[1] <- 'Description'
	colnames(gsea_output_ontology_filtered)[3] <- 'p.adjust'
	gsea_output_ontology_filtered <- gsea_output_ontology_filtered[, c('Description', 'p.adjust', 'ES', 'NES', 'size', 'leadingEdge')]
	gsea_output_ontology_filtered$p.adjust <- round(gsea_output_ontology_filtered$p.adjust, 3)
	gsea_output_ontology_filtered$ES <- round(gsea_output_ontology_filtered$ES, 3)
	gsea_output_ontology_filtered$NES <- round(gsea_output_ontology_filtered$NES, 3)
	gmt_ontology <- paste0(gmt_name, '_ontology')
	gsea_output_ontology_filtered$cluster <- cluster
	gsea_output_ontology_filtered$database <- 'ontology'
	gsea_output_ontology_filtered$comparison <- comp

	if (length(gsea_output_ontology_filtered$NES) > 0)
	{
		bioimages(gsea_output_ontology_filtered, gmt_ontology, ident, pathway_dir, cluster, project, rank, comp)
		fname3=paste0(pathway_dir_table, project, '_', ident, '_', 'cluster_', gmt_ontology, '_filtered_results_', rank, '.csv', sep='')
		print(fwrite(gsea_output_ontology_filtered, file=fname3))
		print('finished ontology')
	}
}
#--------------------------------------------------------------------

# Find DE Genes in each cluster
#--------------------------------------------------------------------
dge <- function(seurat_object, tbl_dir, normalization_method, ident, f_threshold, avg2fc, min_pct, pathway_dir, curated_gene_sets, hallmarks, ontology_gene_sets, project, conserved_genes)
{
	print(normalization_method)

	Idents(object=seurat_object) <- 'celltypes'

	samples <- length(unique(seurat_object@meta.data[['Sample']]))

	if (samples > 1)
	{
		print('calculating difference in gene expression by cluster...')

		clusters <- unique(seurat_object@active.ident) #this is what was here originally
		#clusters <- unique(seurat_object@meta.data[['celltypes']])

		for (c in clusters)
		{
			print(c)
			cluster <- subset(seurat_object, idents = c)
			
			# Conserved Genes -----------------------------------------
			if (conserved_genes == 'y' || conserved_genes == TRUE)
			{
	 			# define conserved output dir
				conserved_dir=paste0(tbl_dir, 'conserved_genes/')
				dir.create(conserved_dir, showWarnings=FALSE)

				# need to check that at least one condition in cluster has >= 3 cells
				condition_flag = 0

				cells_per_condition <- table(cluster@meta.data$Experiment)
				max_cells_per_condition <- max(cells_per_condition)

				if (max_cells_per_condition > 3)
				{
					print('Finding conserved DGEs...')

					set.seed(42)
					project.conserved.markers <- FindConservedMarkers(seurat_object, ident.1=c, min.cells.group = 3, grouping.var='Experiment', verbose=TRUE)
					# check to see if project.conserved.markers is an empty dataframe
					if (nrow(project.conserved.markers) != 0 && ncol(project.conserved.markers) != 0) 
					{
						project.conserved.markers$cluster <- c
						project.conserved.markers <- rownames_to_column(project.conserved.markers, var = 'gene')
						conserved_markers_filename = paste0(conserved_dir, '/', project, '_', ident, '_conservedMarkers_cluster_', c, '.txt')
						write.table(project.conserved.markers, file=conserved_markers_filename, sep='\t', quote=F, col.names=TRUE, row.names=TRUE, append=FALSE)
					}
				}
			}
			# ---------------------------------------------------------

			if (normalization_method == 'sct')
			{
				set.seed(42)
				print('Preping seurat object (sct assay)...')
				cluster <- PrepSCTFindMarkers(object = cluster)
			}

			# commented lines below were here originally -- moved up to obviate repetition
			cluster$celltype.exp <- paste(Idents(cluster), cluster$Experiment, sep = '_')
			cluster$celltype <- Idents(cluster)
			Idents(cluster) <- 'celltype.exp'

			conditions <- unique(cluster$celltype.exp)

			#in the event there is more than 2 conditions
			for (i in conditions)
			{
				ci <- subset(x = cluster, subset = celltype.exp == i)
				cix <- nrow(ci@meta.data)

				for (j in conditions)
				{
				  cj <- subset(x = cluster, subset = celltype.exp == j)
				  cjx <- nrow(cj@meta.data)

				  if (cix > 3 && cjx > 3) #need at least 3 cells per cluster
				  {
						if (i != j)
						{
							comparison <- paste0(i, '_vs_', j)
							print('finding markers for conditions:')
							print(c(i, j))

							table_name = paste(tbl_dir, project, '_', ident, '_final_markers_cluster_', i, '_vs_', j, '.txt', sep='')
							inverted_table_name = paste(tbl_dir, project, '_', ident, '_final_markers_cluster_', j, '_vs_', i, '.txt', sep='')
							print(table_name)

							project.markers = '' # initialize empty var

							if ( file.exists(inverted_table_name) == FALSE) # i+j = j+i
							{
								if (normalization_method == 'sct')
								{
									set.seed(42)
									project.markers <- FindMarkers(cluster, ident.1 = i, ident.2 = j, verbose = TRUE, only.pos = FALSE, slot='data', recorrect_umi=FALSE, logfc.threshold=avg2fc, min.pct=min_pct)
								}

								if (normalization_method == 'standard')
								{
									set.seed(42)
									project.markers <- FindMarkers(cluster, ident.1 = i, ident.2 = j, verbose = TRUE, only.pos = FALSE, slot='data', logfc.threshold=avg2fc, min.pct=min_pct)
								}

								if (length(project.markers[1]) != 0)
								{
									setDT(project.markers, keep.rownames = TRUE)[]
									names(project.markers)[1] <- paste('gene')
									project.markers$cluster <- c

									print('this is the comparison....')
									print(comparison)

									project.markers$comparison <- comparison
									write.table(project.markers, file=table_name, sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

									# get pathways for cluster
									f_gsea(project.markers, c, 'p_val_adj', filtering_threshold, pathway_dir, ident, curated_gene_sets, hallmarks, ontology_gene_sets, project, comparison)
									f_gsea(project.markers, c, 'avg_log2fc', filtering_threshold, pathway_dir, ident, curated_gene_sets, hallmarks, ontology_gene_sets, project, comparison)
								}
							}
						}
					}
				}
			}
		}
	}
}
#--------------------------------------------------------------------

# COMBINE AND CLEAN -------------------------------------------------
combine_clean <- function(pathway_dir, project, ident, tbl_dir)
{
	# combine all marker text files into one
	print('combine all marker text files into one')
	marker_file <- paste0(tbl_dir, project, '_', ident, '_final_markers_total.txt')
	
	# marker_df <- list.files(path = tbl_dir, pattern = '*final_markers_cluster_', full.names = TRUE) %>%
	# 	lapply(read_table) %>% bind_rows
	marker_df <- list.files(path = tbl_dir, pattern = '*final_markers_cluster_', full.names = TRUE) %>%
	  lapply(function(file) {
		df <- read_table(file)
			 # Ensure the 'cluster' column is of type 'character' (you can change to 'factor' if necessary)
				  df$cluster <- as.character(df$cluster)
						return(df)
						  }) %>%
							 bind_rows()
	if (dim(marker_df)[1] > 1) {
		write.table(marker_df, marker_file, quote=FALSE, sep='\t', row.names = FALSE)
	}

	# combine pathway tbls into single csv
	# print('combine pathway tbls into single csv')
	# get list of csv files to delete before combined csv has been created

	# delete old copy of csv file (if exists) -- prevents concat file
	print('delete old copy of csv file (if exists)...')
	csv_file_adj <- paste0(pathway_dir_table, project, '_', ident, '_final_filtered_results_combined_pathways_p_val_adj.csv')
	csv_file_log <- paste0(pathway_dir_table, project, '_', ident, '_final_filtered_results_combined_pathways_avg_log2fc.csv')

	if (file.exists(csv_file_adj))
	{
			file.remove(csv_file_adj)
	}

	if (file.exists(csv_file_log))
	{
			file.remove(csv_file_log)
	}

	# must define cvs list here!!!
	#delete_csv  <- list.files(path = pathway_dir_table, pattern = '\\.csv$', full.names =TRUE)

	# combine all csv files into one
	print('combine all csv files into one (adj.p.value)')
	#csv_df <- list.files(path = pathway_dir_table, pattern = '\\.csv$', full.names =TRUE) %>% 
	# csv_df_adj <- list.files(path = pathway_dir_table, pattern = '*p_val_adj.csv$', full.names =TRUE) %>% 
	# 	lapply(read_csv) %>% bind_rows
	csv_df_adj <- list.files(path = pathway_dir_table, pattern = '*p_val_adj.csv$', full.names = TRUE) %>% 
	  lapply(function(file) {
		      df <- read_csv(file)
				    # Ensure 'cluster' column is of type 'character' (change to 'factor' if you prefer)
					     df$cluster <- as.character(df$cluster)
						      return(df)
								  }) %>% 
								    bind_rows()
	if (dim(csv_df_adj)[1] >= 1) {
		write.csv(csv_df_adj, csv_file_adj)
	}

	print('combine all csv files into one (avglog2FC)')
	#csv_df_log <- list.files(path = pathway_dir_table, pattern = '*avg_log2fc.csv$', full.names =TRUE) %>% 
	#	lapply(read_csv) %>% bind_rows 
	csv_df_log <- list.files(path = pathway_dir_table, pattern = '*avg_log2fc.csv$', full.names = TRUE) %>% 
	  lapply(function(file) {
		      df <- read_csv(file)
				    # Ensure 'cluster' column is of type 'character' (change to 'factor' if you prefer)
					     df$cluster <- as.character(df$cluster)
						      return(df)
								  }) %>% 
								    bind_rows()
	if (dim(csv_df_log)[1] >= 1) {
		write.csv(csv_df_log, csv_file_log)
	}

	# Conserved Genes -----------------------------------------
	if (conserved_genes == 'y')
	{
		print('consolidating conserved DGEs...')

		print('delete old copy of txt file (if exists)...')
		conserved_dir <- paste0(tbl_dir, 'conserved_genes/')
		conserved_combined_file <- paste0(conserved_dir, project, '_', ident, '_combined_conserved_markers.txt')
		print(conserved_combined_file)

		if (file.exists(conserved_combined_file))
		{
			file.remove(conserved_combined_file)
		}

		conserved_files <- list.files(path = conserved_dir, pattern = '*conservedMarkers.*\\.txt', full.names = TRUE) %>% 
			lapply(function(file) 
			{
				# Read file into data frame
				df <- read.table(file, header=TRUE, sep='\t') 
				
				# Check if number of columns is >= 10
				if (ncol(df) > 10) 
				{
					return(df)
				} 
				
				else 
				{
					return(NULL)  # If the file has fewer than 10 columns, return NULL
				}
			}) %>% compact() %>% bind_rows() # remove NULL elements (files w/ < 10 columns) + bind all remaining d.f.s together

		write.table(conserved_files, conserved_combined_file, quote=FALSE, row.names=TRUE, sep='\t')
	}
	# ----------------------------------------------

	# may want to keep the files....
	# delete individual csv files
	#for (d in delete_csv)
	#{
#		if (d != csv_file)
#		{
#			df <- ''
#			#file.remove(d)
#		}
#	}
}
#--------------------------------------------------------------------

# REMOVE UNUSED METADATA --------------------------------------------
seurat_object_reduction <- function(seurat_object, ident, redux_umap, redux_tsne, normalization_method, integration_method)
{
	temp_seurat_object <- seurat_object
	meta_names <- names(seurat_object@meta.data)

	meta_tags <- c('orig.ident', 'nCount_RNA', 'nFeature_RNA','Experiment', 'Sample', 
		'percent.mito', 'percent.ribo', 'S.Score', 'G2M.Score', 'Phase', 'CC.Difference', 
		'nCount_SCT', 'nFeature_SCT', 'celltypes', ident)

	print('removing unused meta data...')

	for (m in names(seurat_object@meta.data))
	{
		if (m %in% meta_tags == FALSE)
		{
			temp_seurat_object[[m]] <- NULL
		}
	}

	pca_reduction <- paste0(normalization_method, '.pca')
	reduction_tags = c(redux_umap, redux_tsne, pca_reduction)

	print('removing unused reductions...')
	for (m in names(seurat_object@reductions))
	{
		if (m %in% reduction_tags == FALSE)
		{
			temp_seurat_object[[m]] <- NULL
			temp_seurat_object@reductions[[m]] <- NULL
		}
	}

	nn_graph <- paste0(normalization_method, '.', integration_method, '_nn')
	snn_graph <- paste0(normalization_method, '.', integration_method, '_snn')
	graph_tags = c(nn_graph, snn_graph)

	print('removing unused graphs...')
	for (m in names(seurat_object@graphs))
	{
		if (m %in% graph_tags == FALSE)
		{
			temp_seurat_object@graphs[[m]] <- NULL
		}
	}

	print(names(temp_seurat_object@meta.data))
	print(names(temp_seurat_object@reductions))
	print(names(temp_seurat_object@graphs))

	return(temp_seurat_object)
}
#--------------------------------------------------------------------

# SAVE SEURAT OBJECT ------------------------------------------------
save_seurat_object <- function(seurat_object, rds_dir, project, final_storage_method, ident)
{
	core_name <- paste0(rds_dir, project, '_final_analyzed_seurat_object')
	SCEO <- '' 

	print('these are the final storage options')
	print(final_storage_method)

	for (s in final_storage_method)
	{
		if (s == 'sceo')
		{
			print('You have requested a Single Cell Experiment object...OK')
			SCEO <- as.SingleCellExperiment(seurat_object)

			if ('rds' %in% final_storage_method) 
			{
				print('You have requested a rds file...OK')
				f_name = paste0(core_name, '_sceo.rds')
				saveRDS(SCEO, file = f_name)
			}
		
			if ('qs' %in% final_storage_method)
			{
				print('saving sceo as qs file')
				f_name = paste0(core_name, '_sceo.qs')
				qsave(SCEO, file = f_name)
			}
		}
	}

	for (s in final_storage_method)
	{
		print(s)

		# saving s.o. in qs format 
		#print('saving object as qs file...')
		#f_name_qs = paste0(core_name, '.qs')
		#qsave(seurat_object, f_name_qs)

		if (s == 'cellphone')
		{
			print('You have requested CellPhoneDB input files...OK')
			normalized_counts <- GetAssayData(seurat_object, assay = 'RNA', slot = 'data')

			cellphoneDB_dir <- paste0(base_directory, '/cellphoneDB/')
			dir.create(cellphoneDB_dir, showWarnings=FALSE)

			normalized_counts_file <- paste0(cellphoneDB_dir, project, '_final_analyzed_cellphoneDB.txt')
			write.table(normalized_counts, file=normalized_counts_file, sep='\t', quote=F, col.names=TRUE, row.names=FALSE)
		}

		if (s == 'cellchat')
		{
			print('You have requested a CellChat object...OK')
			suppressPackageStartupMessages(library(CellChat, lib.loc=lib_path))
			
			cellChat <- createCellChat(object = seurat_object, group.by = 'celltypes', assay = 'RNA')

			if ('rds' %in% final_storage_method) 
			#if (grepl( 'rds', final_storage_method, fixed = TRUE) == TRUE)
			{
				f_name = paste0(core_name, '_cellchat.rds')
				saveRDS(cellChat, file = f_name)
			}

			if ('qs' %in% final_storage_method)
			#if (grepl( 'rds', final_storage_method, fixed = TRUE) == FALSE)
			{
				f_name = paste0(core_name, '_cellchat.qs')
				qsave(cellChat, f_name)
			}
		}

		# Always save as .qs !!!
		print('You have requested a qs file...OK')
		f_name = paste0(core_name, '.qs')
		qsave(seurat_object, f_name)

		if (s == 'rds')
		{
			print('You have requested a rds file...OK')
			f_name = paste0(core_name, '.rds')
			saveRDS(seurat_object, file = f_name)
		}

		if (s == 'cloupe')
		{
			print('You have requested a cloupe object...OK')
			Idents(object=seurat_object) <- 'celltypes'
			suppressPackageStartupMessages(library(loupeR, lib.loc=lib_path))

			print('joining layers...')
			seurat_object <- JoinLayers(seurat_object, assay = 'RNA')
			
			print('getting sparse matrix....')
			sparse_counts <- GetAssayData(seurat_object, assay = 'RNA', layer = 'counts')

			print('getting clusters...')
			clusters <- select_clusters(seurat_object)

			print('getting projections...')
			projections <- select_projections(seurat_object)

			cloupe_name <- paste0(project, '_final_analyzed_seurat_object')
			#create_bugreport_from_seurat(seurat_object)

			create_loupe(count_mat=sparse_counts, projections=projections, clusters=clusters, output_dir = rds_dir, output_name = cloupe_name)
			#create_loupe_from_seurat(seurat_object, output_dir = rds_dir, output_name = cloupe_name)
		}
	}
}
#--------------------------------------------------------------------

# IMPORT DATA
#--------------------------------------------------------------------
print('importing seurat object')
S <- ''

# GET SEURAT OBJECT EXTENSION
# User may supply s.o., may be in RDS format
#--------------------------------------------------------------------
extension = file_ext(analyzed_seurat_object)

if (extension == 'qs')
{
	S = qread(analyzed_seurat_object)
}

if (extension == 'rds')
{
	S = readRDS(analyzed_seurat_object)
}

if (extension != 'rds' & extension != 'qs')
{
	print('The provided seurat file is not a qs or rds file.')
	print('This script will fail. Please save the seurat file as a qs or rds file.')
	stop()
}
#--------------------------------------------------------------------

# CHECK SEURAT VERSION -----------------------------------------------
seurat_version <- Version(S) #e.g., '5.0.2'
seurat_version_u <- as.numeric(unlist(seurat_version)) 

if (seurat_version_u[1] < 5.0)
{
	print('The version of the seurat object is less than 5.0')
	print('The pipeline only works with seurat objects >= 5.*')
	print('Supply a correctly versioned object and try this script again')
	stop()
}

if (seurat_version_u[1] >= 5.0)
{
	print('This is the Seurat Object Version:')
	print(seurat_version)
}
#--------------------------------------------------------------------

# DEFINE ARGUMENTS + VARIABLES
#--------------------------------------------------------------------
assay = ''
filtering_threshold <- as.numeric(filtering_threshold)
ident <- paste0(normalization_method, '.', integration_method, '_snn_res.', resolution)
redux_umap <- paste0(normalization_method, '.', integration_method, '.umap')
redux_tsne <- paste0(normalization_method, '.', integration_method, '.tsne')

# if user is supplying seurat object
if	(provide == 'y')
{
	redux_umap <- umap_reduction
	redux_tsne <- tnse_reduction
}

if (normalization_method == 'sct')
{
	assay = 'SCT'
}

if (normalization_method == 'standard')
{
	assay = 'RNA'
}

# Set Assay + Idents
DefaultAssay(S) <- assay
Idents(object=S) <- ident
#--------------------------------------------------------------------

# DETERMINE IF USER HAS SUPPLIED LIST OF GENES
# ------------------------------------------------
genes <- NULL
markers <- ''
visualization_method <- ''

if (user_gene_file != 'does_not_exist' & file.exists(user_gene_file))
{
	print('user has supplied a marker file')
	markers <- read.table(user_gene_file, header=FALSE) #this is a list of cell type file locations
	markers$V1 <- unique(markers$V1)
	genes = as.vector(markers$V1)

	if (is.null(visualization) == FALSE)
	{
		print('we are assigning visualization methods...')
		visualization_method = single_or_list(visualization)
	}
}
# ------------------------------------------------

# CREATE MSIG LIST -------------------------------
curated_gene_sets <- msigdbr(species=organism, category='C2', subcategory='CP')
hallmarks=msigdbr(species=organism, category='H')
ontology_gene_sets <- msigdbr(species=organism, category='C5')
# ------------------------------------------------

# RENAME CLUSTERS + PROP + DIMPLOTS (still need to update for single sample)?????
# ------------------------------------------------
# if providing final seurat object and no annotation is required
flag = 0
if (provide == 'y' & annotate == 'n')
{
	flag = 1
}

# FUNCTION CALLS ----------------------------------------------------
S <- providing_object(provide, annotate, S, meta_sample, meta_experiment, meta_annotation)

print('Entering rename_and_visualize function...')
S <- rename_and_visualize(S, celltype_assignment_file, ident, genes, markers, redux_umap, redux_tsne, tsne, fig_dir, tbl_dir, flag, visualization_method, project, normalization_method)

# DGEA + FGSEA -----------------------------------
print('Entering dge function...')
dge(S, tbl_dir, normalization_method, ident, filtering_threshold, avg2fc, min.pct, pathway_dir, curated_gene_sets, hallmarks, ontology_gene_sets, project, conserved_genes)

# COMBINE AND CLEAN ------------------------------
print('consolidating output...')
combine_clean(pathway_dir, project, ident, tbl_dir)

# REMOVE UNUSED METADATA
print('removing unused metatdata...')
S <- seurat_object_reduction(S, ident, redux_umap, redux_tsne, normalization_method, integration_method)

#, SAVE, UPDATE, SO
print('savings...')
save_seurat_object(S, rds_dir, project, final_storage_method, ident)
