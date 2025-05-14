# Author:	E. Reichenberger
# Date:		7.31.2024

# Purpose:	Create cluster UMAP/TsNE images, create cells/cluster proportions plots+table, and find conserved and DE markers (upregulated only)
# for all clusters across all possible resolutions, normalization, and integration methods.

#(r message=FALSE)

args = commandArgs(trailingOnly=TRUE)
project = '' 
lib_path = ''
storage = ''
normalization_method = '' 
integration_method = '' 
resolution = ''
conserved_genes = ''
analyzed_seurat_object = ''
processes = ''
tsne_plot = ''
report_table_path = ''
user_gene_file <- '' #may not exist
visualization <- '' #may not exist

nargs = 13
# test for correct number of arguments
if (length(args) < nargs) 
{
	  stop('At least 13 arguments must be supplied.', call.=FALSE)
} 

if (length(args)==nargs) 
{
	project = args[1] 
	lib_path = args[2]
	storage = args[3]
	normalization_method = args[4] 
	integration_method = args[5] 
	resolution = args[6]
	conserved_genes = args[7]
	analyzed_seurat_object = args[8]
	processes = args[9]
	tsne_plot = args[10]
	report_table_path = args[11]
	user_gene_file = args[12]
	visualization = args[13]
}

suppressPackageStartupMessages(library(RColorBrewer, lib.loc=lib_path))
suppressPackageStartupMessages(library(dplyr, lib.loc=lib_path))
suppressPackageStartupMessages(library(ggrepel, lib.loc=lib_path))
suppressPackageStartupMessages(library(Seurat, lib.loc=lib_path))
suppressPackageStartupMessages(library(ggplot2, lib.loc=lib_path))
suppressPackageStartupMessages(library(reshape2, lib.loc=lib_path))
suppressPackageStartupMessages(library(data.table, lib.loc=lib_path))
suppressPackageStartupMessages(library(qs, lib.loc=lib_path))
suppressPackageStartupMessages(library(future, lib.loc=lib_path))
suppressPackageStartupMessages(library(progressr, lib.loc=lib_path))
suppressPackageStartupMessages(library(presto, lib.loc=lib_path))
suppressPackageStartupMessages(library(tidyverse, lib.loc=lib_path)) #added 12.3.24 ERR

# PARALLEL w/ FUTURE + SET SEED
#--------------------------------------------------------------------
options(future.globals.maxSize = 210000 * 1024^2) #may way to make that a variable that user can increase if there is a failure or base it on the dataset size???
plan(multisession(workers = as.integer(processes)))

set.seed(42)
#--------------------------------------------------------------------

# FUNCTION: determine if argument is single value or list
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
	}
}
# --------------------------------------------------------------------

# SET CONFIG VAR AS LIST 
# --------------------------------------------------------------------
resolution = single_or_list(resolution)
#resolution = as.numeric(single_or_list(resolution))
normalization_method = single_or_list(normalization_method)
integration_method = single_or_list(integration_method)
print(resolution)
print(normalization_method)
print(integration_method)
print(visualization)
visualization_method <- ''
if (!is.null(visualization))
{
	visualization_method = single_or_list(visualization)
}
print(visualization_method)

# --------------------------------------------------------------------

# CREATE DIRECTORIES
#--------------------------------------------------------------------
base_directory=paste('data/endpoints/', project, '/analysis/normalization/', sep='')

dir.create(base_directory, showWarnings=FALSE)
dir.create(report_table_path, recursive = TRUE, showWarnings = FALSE)

for (n in normalization_method)
{
	normalization_base = paste0(base_directory, n)
	dir.create(normalization_base, showWarnings=FALSE)

	dir_figures=paste(normalization_base, '/figures/', sep='')
	dir_tables=paste(normalization_base, '/tables/', sep='')

	dir.create(dir_figures, showWarnings=FALSE)
	dir.create(dir_tables, showWarnings=FALSE)
}
#--------------------------------------------------------------------

# FUNCTION: calculate proportions, plot reduction dimplots, find upregulated markers 
# --------------------------------------------------------------------

proportions_UMAP_DGE <- function(seurat_object, num_samples, visi, genes=genes, markers=markers, resolut=resolution, normal=normalization_method, integration=integration_method, tsne=tsne_plot, start_directory=base_directory, store=storage)
{
	count_normal = 1
	for (n in normal)
	{
		print(n)
		# set assay, if sct <- SCT, standard <- RNA
		assay <- ''

		dir_fig=paste0(start_directory, n, '/figures/')
		dir_table=paste0(start_directory, n, '/tables/')

		print('Setting assay....')

		if (n == 'sct')
		{
			assay = 'SCT'
		}

		if (n == 'standard')
		{
			assay = 'RNA'
		}
		# --------------------------------------

		DefaultAssay(seurat_object) <- assay

		count_integration = 1
		for (i in integration)
		{
			print(i)
			
			# set reductions
			print('Setting reduction(s)...')
			redux_umap <- paste0(n, '.', i, '.umap')
			redux_tsne <- paste0(n, '.', i, '.tsne')
			print(redux_umap)
			print(redux_tsne)

			if (num_samples == 1) #assuming integration was set to 'pca' 
			{
				redux_umap <- paste0(n, '.umap')
				redux_tsne <- paste0(n, '.tsne')
			}

			for (r in resolut)
			{
				# was getting odd error were i = 45 & 47????
				if ( n != 'standard' && i != 'sct')
				{
					n = normal[count_normal]
				}

				if ( i != 'cca' && i != 'harmony' && i != 'rpca')
				{
					i = integration[count_integration]
				}

				r <- as.numeric(r)
				name = paste0(n, '.', i, '_snn_res.', r)
				print(name)
				Idents(object=seurat_object) <- name
				
				# z-scores by cluster  ------------------------------------
				print('calculating z-scores across all clusters...')
				ae <- AggregateExpression(object = seurat_object, group.by = name)$RNA
				z_ae <- as.data.frame(scale(ae))
				z_ae <- tibble::rownames_to_column(z_ae, "gene")
				z_ae <- reshape2::melt(z_ae)
				z_ae$variable <- gsub('g', '', z_ae$variable)
				colnames(z_ae) <- c('gene', 'cluster', 'z.score') # changed to z.score 2.25.25
				z_ae$z.score <- round(z_ae$z.score, 3) # added 2.25.25

				z_scores_report = paste(report_table_path, '/', project, '_z_scores.', name, '.txt', sep='')
				write.table(z_ae, file=z_scores_report, sep='\t', quote=F, col.names=TRUE, row.names=FALSE)
				# ---------------------------------------------------------

				# DimPlots -----------------------------------------------
				print('Making dimplots...')
				dimplot_fig_name = paste(dir_fig, '/', project, '_DimPlot_Proportions_', name, '.pdf', sep='')
				pdf(file=dimplot_fig_name, onefile=TRUE, width=11, height=8.5)

				print(DimPlot(seurat_object, group.by = name, reduction = redux_umap, label=TRUE, repel=TRUE))
				print(DimPlot(seurat_object, group.by = name, reduction = redux_umap, split.by = 'Experiment', repel=TRUE, label=TRUE)) + NoLegend()
				print(DimPlot(seurat_object, group.by = name, reduction = redux_umap, split.by = 'Sample', repel=TRUE, label=TRUE)) + NoLegend()
				print(DimPlot(seurat_object, reduction = redux_umap, group.by = 'Phase'))

				if (tsne == 'y')
				{
					print('Making dimplots, tsne edition...')
					print(DimPlot(seurat_object, group.by = name, reduction = redux_tsne, label=TRUE, repel=TRUE))
					print(DimPlot(seurat_object, group.by = name, reduction = redux_tsne, split.by = 'Experiment', repel=TRUE, label=TRUE)) + NoLegend()
					print(DimPlot(seurat_object, group.by = name, reduction = redux_tsne, split.by = 'Sample', repel=TRUE, label=TRUE)) + NoLegend()
				}
				# ---------------------------------------------------------

				# Proportions ---------------------------------------------
				print('Calculating proportions by experiment...')
				number_perCluster_experiment <- table(seurat_object@meta.data$Experiment, seurat_object@meta.data[[name]])
				number_perCluster_experiment_prop1 <- round(prop.table(number_perCluster_experiment, margin = 1) * 100, 1)

				# combine tables as matrix
				X = as.matrix(number_perCluster_experiment)
				Y = as.matrix(number_perCluster_experiment_prop1)
				Z = matrix(paste0(X, " (", Y, "%", ")"), nrow = nrow(X), dimnames = dimnames(X))

				npce <- reshape2::melt(Z)
				colnames(npce) <- c('expCond', 'cluster', 'numCells')

				proportion_table_name_experiment = paste(report_table_path, '/', project, '_clusterProportions.experiment_', name, '.txt', sep='')
				write.table(npce, file=proportion_table_name_experiment, sep='\t', quote=F, col.names=TRUE, row.names=TRUE, append=FALSE)

				print('Calculating proportions by sample...')
				number_perCluster <- table(seurat_object@meta.data$Sample, seurat_object@meta.data[[name]])

				proportion_table_name = paste(dir_table, '/', project, '_clusterProportions_', name, '.txt', sep='')
				write.table(number_perCluster, file=proportion_table_name, sep='\t', quote=F, col.names=TRUE, row.names=TRUE, append=FALSE)

				# will need to melt table
				df <- reshape2::melt(number_perCluster)
				names(df)[1] <- 'Sample'
				names(df)[2] <- 'Cluster'

				df$Cluster <- as.factor(df$Cluster) 

				print(ggplot(df, aes(x=Sample, y=value, fill=Cluster)) +
					geom_bar(stat='identity', position='fill') + 
					ylab('Cluster Proportion') + xlab('Sample') + guides(fill=guide_legend(title='Cluster')))
				# ---------------------------------------------------------

				# Visualize user-defined genes ----------------------------
				if (is.null(genes) == FALSE)
				{
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

							if (!is.null(visi))
							{
								for (v in visi)
								{
									if (v == 'feature')
									{
										print('making (remove dot, add feature) feature plot...')
										try(print(FeaturePlot(seurat_object, label=TRUE, repel=TRUE, label.size=2.5, features=c(t), reduction=redux_umap)))
									}

									if (v == 'dot')
									{
										print('making dot plot...')
										try(print(DotPlot(seurat_object, features=c(t)) + xlab('Genes') + ylab('Clusters') + RotatedAxis()))
									}

									if (v == 'violin')
									{
										print('making violin plot...')
										try(plot(VlnPlot(seurat_object, features=c(t)) + xlab(label = 'Clusters') + ylab('LogNormalized Expression')))
									}

									if (v == 'ridge')
									{
										print('making ridge plot...')
										try(print(RidgePlot(seurat_object, features=t, stack=TRUE) + ylab(label='Clusters')))
									}
								}
							}
						}
					}

					else
					{
						print('Twelve or few genes were supplied')
						print('Making plots highlighting user-provided genes...')
						print(genes)
						print(visi)

						if (!is.null(visi))
						{
							for (v in visi)
							{
								if (v == 'feature')
								{
									print('making feature plot...')
									try(print(FeaturePlot(seurat_object, label=TRUE, repel=TRUE, label.size=2.5, features=c(genes), reduction=redux_umap)))
								}

								if (v == 'dot')
								{
									print('making dot plot...')
									try(print(DotPlot(seurat_object, features=c(genes)) + xlab('Genes') + ylab('Clusters') + RotatedAxis()))
								}

								if (v == 'violin')
								{
									print('making violin plot...')
									try(plot(VlnPlot(seurat_object, features=c(genes)) + xlab(label = 'Clusters') + ylab('LogNormalized Expression')))
								}

								if (v == 'ridge')
								{
									print('making ridge plot...')
									try(print(RidgePlot(seurat_object, features=genes, stack=TRUE) + ylab(label='Clusters')))
								}
							}
						}
					}
				}

				dev.off() # dimplot + proportions + visualize user-provided genes
				# ---------------------------------------------------------

				# DGEA -------------------------
				if (n == 'sct')
				{
					print('Preping seurat object (sct assay)...')
					set.seed(42)
					seurat_object <- PrepSCTFindMarkers(object = seurat_object)
				}

				print('Finding DGEs...')
				handlers(global = TRUE)
				set.seed(42)
				project.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, min.cells.group = 3, slot='data')
				print(names(project.markers))
				print(summary(project.markers))

				if (length(project.markers$cluster) > 0)
				{
					markers_filename = paste(dir_table, '/', project, '_markers_', name, '.txt', sep='')
					print(markers_filename)
					project.markers %>% group_by(cluster) 
					write.table(project.markers, file=markers_filename, sep='\t', quote=F, col.names=TRUE, row.names=FALSE)

					project.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
					top100 <- project.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)

					markers_filename100 = paste(dir_table, '/', project, '_top100_markers_avg_log2FC_', name, '.txt', sep='')
					print(markers_filename100)
					write.table(top100, file=markers_filename100, sep='\t', quote=F, col.names=TRUE, row.names=FALSE)

					# for shiny app
					markers_filename100_report = paste(report_table_path, '/', project, '_top100_markers_avg_log2FC_', name, '.txt', sep='')
					write.table(top100, file=markers_filename100_report, sep='\t', quote=F, col.names=TRUE, row.names=FALSE)
				}
  			# ---------------------------------------------------------

				# Conserved Genes -----------------------------------------
				if (conserved_genes == 'y' & num_samples > 1)
				{
					conserved_dir=paste0(start_directory, n, '/tables/conserved_genes/')
					dir.create(conserved_dir, showWarnings=FALSE)

					conserved_dir=paste0(start_directory, n, '/tables/conserved_genes/', i)
					dir.create(conserved_dir, showWarnings=FALSE)

					conserved_dir=paste0(start_directory, n, '/tables/conserved_genes/', i, '/', r)
					dir.create(conserved_dir, showWarnings=FALSE)

					print('Finding conserved DGEs...')

					cluster_count <- levels(seurat_object@meta.data[[name]]) # anticipate issues, may need [["name"]]
					print(levels(seurat_object@meta.data[[name]])) 

					for (count in cluster_count)
					{
						count <- as.numeric(count)
						print(count)

						# added 1.17.25 ERR
						# Need to check that at least 3 cells exist in one experimental group to prevent failure
						cluster_subset <- subset(seurat_object, idents = count)
						cells_per_condition <- table(cluster_subset@meta.data$Experiment)
						max_cells_per_condition <- max(cells_per_condition)

						if (max_cells_per_condition > 3)
						{
							set.seed(42)
							project.conserved.markers <- FindConservedMarkers(seurat_object, ident.1=count, ident.2=NULL, min.cells.group = 3, grouping.var='Experiment', verbose=TRUE)

							conserved_markers_filename = paste0(conserved_dir, '/', project, '_conservedMarkers_cluster_', count, '_', name, '.txt')
							write.table(project.conserved.markers, file=conserved_markers_filename, sep='\t', quote=F, col.names=TRUE, row.names=TRUE, append=FALSE)
						}
					}
				}
				# ---------------------------------------------------------
			}

			count_integration = count_integration + 1
		}

		count_normal = count_normal + 1

		if (store == 'rds')
		{
			print('saving object as RDS')
			filename <- paste0('data/endpoints/', project, '/analysis/RDS/', project, '_analyzed_seurat_object.RDS')
			saveRDS(seurat_object, file=filename)
		}
	}
}
# --------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------
# IMPORT DATA
print('importing data...')
S = qread(analyzed_seurat_object)
# ------------------------------------------------

# GET NUMBER OF SAMPLES
num_samples <- length(unique(S@meta.data[['Sample']]))
# ------------------------------------------------

# DETERMINE IF USER HAS SUPPLIED LIST OF GENES
genes <- NULL 
markers <- ''
print(user_gene_file)

if (user_gene_file != 'does_not_exist' & file.exists(user_gene_file) == TRUE)
{
	print('User has supplied a marker file')
	markers <- read.table(user_gene_file, header=FALSE) #this is a list of cell type file locations
	genes = as.vector(unique(markers$V1))
	print('...genes will be visualized...')
	print(genes)


}
# ------------------------------------------------

# MAKE GOODIES
proportions_UMAP_DGE(S, num_samples, visualization_method, genes, markers)
# --------------------------------------------------------------------------------------------
