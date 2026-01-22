#!/usr/bin/env Rscript
# Author:  E. Reichenberger
# Modified by: M. Brown
# Mod note: M. Brown add preprocessing directive, converted positional args to flags with input checks and a default lib path for compatibility with argparse and nextflow integration

# Date:    7.31.2024

# Purpose:  Create cluster UMAP/TsNE images, create cells/cluster proportions plots+table, and find conserved and DE markers (upregulated only)
# for all clusters across all possible resolutions, normalization, and integration methods.

#(r message=FALSE)

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
  make_option(c("--project"), type="character",
              help="Project name, used to create output directories and file names."),
  make_option(c("--storage"), type="character", default='rds',
              help="Storage format for the analyzed Seurat object, either 'rds' or 'qs'."),
  make_option(c("--normalization_method"), type="character",
              help="Normalization method(s) to use, comma-separated. Options are 'sct' and 'standard'."),
  make_option(c("--integration_method"), type="character",
              help="Integration method(s) to use, comma-separated. Options are 'cca', 'harmony', 'rpca', and 'sct'."),
  make_option(c("--resolution"), type="character",
              help="Resolution(s) to use for clustering, comma-separated."),
  make_option(c("--conserved_genes"), type="character", default='n',
              help="Whether to find conserved genes across experiments. Options are 'y' or 'n'. Default is 'n'."),
  make_option(c("--analyzed_seurat_object"), type="character",
              help="Path to the analyzed Seurat object file (RDS or QS format)."),
  make_option(c("--processes"), type="integer", default=1,
              help="Number of processes to use for parallel computation. Default is 1."),
  make_option(c("--tsne_plot"), type="character", default='n',
              help="Whether to create tSNE plots. Options are 'y' or 'n'. Default is 'n'."),
  make_option(c("--report_table_path"), type="character",
              help="Path to save report tables. Default is 'data/endpoints/{project}reports/'."),
  make_option(c("--user_gene_file"), type="character",
              help="Path to a user-defined gene file for visualization. Default is 'does_not_exist'."),
  make_option(c("--visualization"), type="character",
              help="Visualization methods for user-defined genes, comma-separated. Options are 'feature', 'dot', 'violin', and 'ridge'.")
)  

# Since the last args is positional, object sticks the options in a separate key
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments=TRUE, args=args)
project <- if (is.null(opt$options$project)) stop("--project is required. See --help for all opts") else opt$options$project
storage <- opt$options$storage
normalization_method <- if (is.null(opt$options$normalization_method)) stop("--normalization_method is required. See --help for all opts") else opt$options$normalization_method
integration_method <- if (is.null(opt$options$integration_method)) stop("--integration_method is required. See --help for all opts") else opt$options$integration_method
resolution <- if (is.null(opt$options$resolution)) stop("--resolution is required. See --help for all opts") else opt$options$resolution
conserved_genes <- if (is.null(opt$options$conserved_genes)) stop("--conserved_genes is required. See --help for all opts") else opt$options$conserved_genes
analyzed_seurat_object <- if (is.null(opt$options$analyzed_seurat_object) || !file.exists(opt$options$analyzed_seurat_object)) stop("--analyzed_seurat_object is required and must be a valid file path. See --help for all opts")  else opt$options$analyzed_seurat_object
processes <- opt$options$processes
tsne_plot <- if (is.null(opt$options$tsne_plot)) stop("--tsne_plot is required. See --help for all opts") else opt$options$tsne_plot
report_table_path <- if (is.null(opt$options$report_table_path)) stop("--report_table_path is required. See --help for all opts") else opt$options$report_table_path
user_gene_file <- if (is.null(opt$options$user_gene_file) || !file.exists(opt$options$user_gene_file)) 'does_not_exist' else opt$options$user_gene_file
visualization <- if (is.null(opt$options$visualization)) stop("--visualization is required. See --help for all opts") else opt$options$visualization


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
suppressPackageStartupMessages(library(tidyverse, lib.loc=lib_path))

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
base_directory=file.path('data/endpoints', project, 'analysis/normalization')

dir.create(base_directory, showWarnings=FALSE, recursive=TRUE)
dir.create(report_table_path, recursive = TRUE, showWarnings = FALSE)

for (n in normalization_method)
{
  normalization_base = file.path(base_directory, n)
  dir.create(normalization_base, showWarnings=FALSE)

  dir_figures=file.path(normalization_base, 'figures')
  dir_tables=file.path(normalization_base, 'tables')

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

    dir_fig=file.path(start_directory, n, 'figures')
    dir_table=file.path(start_directory, n, 'tables')

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
        z_ae$z.score <- round(z_ae$z.score, 3)

        z_scores_report = paste(report_table_path, '/', project, '_z_scores.', name, '.txt', sep='')
        write.table(z_ae, file=z_scores_report, sep='\t', quote=F, col.names=TRUE, row.names=FALSE)
        # ---------------------------------------------------------

				# phase information ---------------------------------------
				# Experiment
				print('making table of phases across all clusters (by experiment)...')
				cc_phase_cluster <- as.data.frame(table(seurat_object@meta.data$Experiment, seurat_object@meta.data[[name]], seurat_object@meta.data$Phase))
				colnames(cc_phase_cluster) <- c('Experiment', 'Cluster', 'Phase', 'Frequency')

				cc_phase_cluster_prop <- as.data.frame(prop.table(table(seurat_object@meta.data$Experiment, seurat_object@meta.data[[name]], seurat_object@meta.data$Phase), margin = 1) * 100)
				colnames(cc_phase_cluster_prop) <- c('Experiment', 'Cluster', 'Phase', 'Frequency')
				cc_phase_cluster$Cluster <- factor(cc_phase_cluster$Cluster, levels = sort(unique(cc_phase_cluster$Cluster))) #make factor for human-readable ordering 1.22.26

				# consider adding combined frequency output if added to html output
				#cc_phase_cluster$Combined_Frequency <- paste0(cc_phase_cluster$Frequency, ' (', round(cc_phase_cluster_prop$Frequency, 1), '%)')
				#cc_phase_cluster <- cc_phase_cluster[, c('Experiment', 'Cluster', 'Phase', 'Combined_Frequency')]
				#colnames(cc_phase_cluster) <- c('Experiment', 'Cluster', 'Phase', 'Frequency')
				#cc_phase_cluster_sample$Combined_Frequency <- paste0(cc_phase_cluster_sample$Frequency, ' (', round(cc_phase_cluster_sample_prop$Frequency, 1), '%)')

				# write table
				cc_phase_cluster <- cc_phase_cluster %>% arrange(Cluster, Phase)
				phase_table = paste(report_table_path, '/', project, '_phase_experiment_', name, '.txt', sep='')
				write.table(cc_phase_cluster, file=phase_table, sep='\t', quote=F, col.names=TRUE, row.names=FALSE)

				# Sample
				print('making table of phases across all clusters (by sample)...')
				cc_phase_cluster_sample <- as.data.frame(table(seurat_object@meta.data$Sample, seurat_object@meta.data[[name]], seurat_object@meta.data$Phase))
				colnames(cc_phase_cluster_sample) <- c('Sample', 'Cluster', 'Phase', 'Frequency')
				cc_phase_cluster_sample$Cluster <- factor(cc_phase_cluster_sample$Cluster, levels = sort(unique(cc_phase_cluster_sample$Cluster))) # make factor for human-readable ordering 1.22.26

				# write table
				phase_table_sample = paste(report_table_path, '/', project, '_phase_sample_', name, '.txt', sep='')
				write.table(cc_phase_cluster_sample, file=phase_table_sample, sep='\t', quote=F, col.names=TRUE, row.names=FALSE)
				# ---------------------------------------------------------

				# DimPlots ------------------------------------------------
				print('Making dimplots...')
				dimplot_fig_name = paste(dir_fig, '/', project, '_DimPlot_Proportions_', name, '.pdf', sep='')
				pdf(file=dimplot_fig_name, onefile=TRUE, width=11, height=8.5)

				print(DimPlot(seurat_object, group.by = name, reduction = redux_umap, label=TRUE, repel=TRUE))
				print(DimPlot(seurat_object, group.by = name, reduction = redux_umap, split.by = 'Experiment', repel=TRUE, label=TRUE)) + NoLegend()
				print(DimPlot(seurat_object, group.by = name, reduction = redux_umap, split.by = 'Sample', repel=TRUE, label=TRUE)) + NoLegend()
				print(DimPlot(seurat_object, reduction = redux_umap, group.by = 'Phase'))
				print(DimPlot(seurat_object, reduction = redux_umap, group.by = 'Phase', split.by='Experiment'))
				print(DimPlot(seurat_object, reduction = redux_umap, group.by = 'Phase', split.by='Sample'))

				print(ggplot(cc_phase_cluster_sample, aes(x=Cluster, y=Frequency, shape=Phase, color=Sample)) + 
				  geom_point(size=3) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
				  labs(x='Cluster', y='Absolute Frequency', shape='Cell Cycle Phase', color='Sample'))

				print(ggplot(cc_phase_cluster, aes(x=Cluster, y=Frequency, shape=Phase, color=Experiment)) + 
				  geom_point(size=3) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
				  labs(x='Cluster', y='Absolute Frequency', shape='Cell Cycle Phase', color='Experiment'))

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

        proportion_table_name_experiment = file.path(report_table_path, paste0(project, '_clusterProportions.experiment_', name, '.txt'))
        write.table(npce, file=proportion_table_name_experiment, sep='\t', quote=F, col.names=TRUE, row.names=TRUE, append=FALSE)

        print('Calculating proportions by sample...')
        number_perCluster <- table(seurat_object@meta.data$Sample, seurat_object@meta.data[[name]])

        proportion_table_name = file.path(dir_table, paste0(project, '_clusterProportions_', name, '.txt'))
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
          markers_filename = file.path(dir_table, paste0(project, '_markers_', name, '.txt'))
          print(markers_filename)
          project.markers %>% group_by(cluster) 
          write.table(project.markers, file=markers_filename, sep='\t', quote=F, col.names=TRUE, row.names=FALSE)

          project.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
          top100 <- project.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)

          markers_filename100 = file.path(dir_table, paste0(project, '_top100_markers_avg_log2FC_', name, '.txt'))
          print(markers_filename100)
          write.table(top100, file=markers_filename100, sep='\t', quote=F, col.names=TRUE, row.names=FALSE)

          # for shiny app
          markers_filename100_report = file.path(report_table_path, paste0(project, '_top100_markers_avg_log2FC_', name, '.txt'))
          write.table(top100, file=markers_filename100_report, sep='\t', quote=F, col.names=TRUE, row.names=FALSE)
        }
        # ---------------------------------------------------------

        # Conserved Genes -----------------------------------------
        if (conserved_genes == 'y' & num_samples > 1)
        {
          conserved_dir=file.path(paste0(start_directory, n), 'tables/conserved_genes')
          dir.create(conserved_dir, showWarnings=FALSE)

          conserved_dir=file.path(paste0(start_directory, n), 'tables/conserved_genes', i)
          dir.create(conserved_dir, showWarnings=FALSE)

          conserved_dir=file.path(paste0(start_directory, n), 'tables/conserved_genes', i, r)
          dir.create(conserved_dir, showWarnings=FALSE)

          print('Finding conserved DGEs...')

          cluster_count <- levels(seurat_object@meta.data[[name]]) # anticipate issues, may need [["name"]]
          print(levels(seurat_object@meta.data[[name]])) 

          for (count in cluster_count)
          {
            count <- as.numeric(count)
            print(count)

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
      rds_path = file.path('data/endpoints', project, 'analysis/RDS')
      dir.create(rds_path, showWarnings=FALSE, recursive=TRUE)
      filename <- file.path(rds_path, paste0(project, '_analyzed_seurat_object.RDS') )
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
