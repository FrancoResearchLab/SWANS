#!/usr/bin/env Rscript
# Author:	E. Reichenberger
# Date:		2.16.2021

# Purpose: 	Remove contaminant RNA w/ SoupX.

# this script does not work for filtered only output...yet.
#sample_name = basename(sample) #https://stackoverflow.com/questions/9693877/how-do-i-extract-a-file-folder-name-only-from-a-path

sample <- ''
lib_path <- ''
data_type <- ''
project <- ''
soupX_input_path <- ''
soupX_output_path <- ''
starter_data <- ''

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 7) {
  stop('At least six arguments must be supplied.', call.=FALSE)
} else if (length(args)==7) {
	lib_path = args[1]
	sample = args[2] # name of sample
	data_type = args[3] # starting data_type (out, filtered, h5)
	project = args[4] #project name
	soupX_input_path = args[5]
	soupX_output_path = args[6]
	starter_data = args[7]
}

library('SoupX', lib.loc=lib_path)
library('DropletUtils', lib.loc=lib_path)
library('ggplot2', lib.loc=lib_path)

set.seed(42)

# CREATE OUT PATH
soupX_sample_path <- file.path(soupX_output_path, project, sample, "soupX")
dir.create(soupX_sample_path, recursive=TRUE, showWarnings=FALSE)

# contamination plot function
#--------------------------------------------------------------------
contam_plot = function(sc)
{
	fit  <- sc$fit
	dd   <- fit$dd$rhoEst   # all individual contamination estimates
	#post <- fit$posterior   # the averaged posterior density data.frame

	x_vals <- seq(0, 0.5, length.out = length(fit$posterior))
    post <- data.frame(x = x_vals, y = fit$posterior)

	print('fit\n')
	print(summary(fit))
	print('dd\n')
	print(summary(dd))
	print('post\n')
	print(summary(post))

	p_contam <- ggplot(post, aes(x = x, y = y)) +
		geom_line(color = '#440154FF', size = 1.2) +  # blue solid line
		geom_vline(xintercept = fit$rhoEst, color = '#26828EFF', size = 1, linetype = 'solid') +  # red estimate line
		stat_function(fun = function(x) 
		{
			dgamma(x, shape = (fit$priorRho / fit$priorRhoStdDev)^2, rate  = fit$priorRho / (fit$priorRhoStdDev^2))
		},
		color = '#C0C0C0', linetype = 'dashed', size = 1) +  # gray dashed prior
		labs(
			x = expression('Contamination Fraction ('*rho*')'),
			y = 'Posterior Density',
			title = 'Estimated Contamination Fraction',
			subtitle = paste0('Estimated ρ = ', round(fit$rhoEst, 3))
		) +
		theme_classic(base_size = 14) +
		theme(
			plot.title = element_text(face = 'bold', size = 16),
			plot.subtitle = element_text(size = 12, margin = margin(b = 10)),
			axis.title = element_text(face = 'bold')
		)
	
	return(p_contam)
}
#--------------------------------------------------------------------


# CLEAN DATA: outs #
#--------------------------------------------------------------------
soupify_outs <- function(in_path, out_path)
{
	print('loading')
	print(in_path)
	print(out_path)
	sc <- load10X(in_path)

	print('estimating')
	sc=autoEstCont(sc)

	print('adjusting')
	out=adjustCounts(sc, roundToInt=TRUE)
	print('writing')
	write10xCounts(out_path, out, version='3', overwrite=TRUE)

	png(filename = paste0(out_path, '/', project, '_', sample, '_', 'contam_plot.png'), height = 2000, width = 2700, res=300)
	print(contam_plot(sc))
	dev.off()
}

if (data_type == 'outs')
{
	print('outs')

	cellranger_data = paste(soupX_input_path, 'outs/', sep='')
	print(cellranger_data)
	soupify_outs(cellranger_data, soupX_sample_path)
}
#--------------------------------------------------------------------

# CLEAN DATA: no cluster information  #
#--------------------------------------------------------------------
soupify_noclusters <- function(sample, in_path, out_path)
{
	library('Seurat', lib.loc=lib_path)

	print(in_path)
	print(out_path)

	r = in_path
	f = out_path

	if (starter_data == 'cellranger' || starter_data == 'fastq')
	{
		r=file.path(in_path, 'outs/raw_feature_bc_matrix/')
		f=file.path(in_path, 'outs/filtered_feature_bc_matrix/')
	}

	raw <- Read10X(data.dir=r)
	filt <- Read10X(data.dir=f)

	filt2 <- CreateSeuratObject(counts=filt , project=project, min.cells=0, min.features=0)
	filt.object <- NormalizeData(object=filt2)
	filt.object <- FindVariableFeatures(object=filt.object)
	filt.object <- ScaleData(object=filt.object)
	filt.object <- RunPCA(object=filt.object)
	filt.object <- FindNeighbors(object=filt.object)
	filt.object <- FindClusters(object=filt.object)
	filt.object <- RunTSNE(object=filt.object)
	clusters <- filt.object@meta.data$seurat_clusters
	names(clusters)<-names(filt.object@active.ident)

	raw <-Read10X(data.dir=r)
	filt <-Read10X(data.dir=f)
	sc=SoupChannel(raw,filt)
	sc=setClusters(sc,clusters)

	sc=autoEstCont(sc, doPlot=TRUE)

	out=adjustCounts(sc, roundToInt=TRUE)
	write10xCounts(out_path, out, version='3', overwrite = TRUE)

	png(filename = file.path(out_path, paste0('/', project, '_', sample, '_', 'contam_plot.png')), height = 2000, width = 2700, res=300)
	print(contam_plot(sc))
	dev.off()
}

# users need sub-directories (57,58), and no outs dir
if (data_type == 'no_clusters')
{
	soupify_noclusters(sample, soupX_input_path, soupX_sample_path)
}
#--------------------------------------------------------------------

# CLEAN DATA: H5
#--------------------------------------------------------------------
soupify_h5 <- function(in_path, out_path)
{
	library('Seurat', lib.loc=lib_path)

	r=paste(in_path, 'raw_feature_bc_matrix.h5', sep='')
	f=paste(in_path, 'filtered_feature_bc_matrix.h5', sep='')

	raw <-Read10X_h5(r, use.names=TRUE)
	filt <-Read10X_h5(f, use.names=TRUE)

	filt2 <- CreateSeuratObject(counts=filt , project=project, min.cells=0, min.features=0)
	filt.object <- NormalizeData(object=filt2)
	filt.object <- FindVariableFeatures(object=filt.object)
	filt.object <- ScaleData(object=filt.object)
	filt.object <- RunPCA(object=filt.object)
	filt.object <- FindNeighbors(object=filt.object)
	filt.object <- FindClusters(object=filt.object)
	filt.object <- RunTSNE(object=filt.object)
	clusters <- filt.object@meta.data$seurat_clusters
	names(clusters) <- names(filt.object@active.ident)

	sc = SoupChannel(raw,filt)
	sc = setClusters(sc,clusters)
	sc=autoEstCont(sc)
	out = adjustCounts(sc, roundToInt=TRUE)
	write10xCounts(out_path, out, version='3', overwrite = TRUE)

	png(filename = file.path(out_path, paste0('/', project, '_', sample, '_', 'contam_plot.png')), height = 2000, width = 2700, res=300)
	print(contam_plot(sc))
	dev.off()
}

# users need files (98,99), and no outs dir
if (data_type == 'h5') 
{
	print('h5')

	#cellranger_data = paste(soupX_input_path, 'outs/', sep='')
	#soupify_outs(cellranger_data, soupX_sample_path)

	soupify_noclusters(sample, soupX_input_path, soupX_output_path)
}
