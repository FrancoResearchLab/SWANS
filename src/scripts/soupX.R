# Author:	E. Reichenberger
# Date:		2.16.2021

# Purpose: 	Remove contaminant RNA w/ SoupX.

# this script does not work for filtered only output...yet.
#sample_name = basename(sample) #https://stackoverflow.com/questions/9693877/how-do-i-extract-a-file-folder-name-only-from-a-path

sam <- ''
lib_path <- ''
data_type <- ''
project <- ''
soupX_input_path <- ''
soupX_ouput_path <- ''
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 6) {
  stop('At least six arguments must be supplied.', call.=FALSE)
} else if (length(args)==6) {
	lib_path = args[1]
	sam = args[2] #name of sample
	data_type = args[3] # starting data_type (out, filtered, h5)
	project = args[4] #project name
	soupX_input_path = args[5]
	soupX_ouput_path = args[6]
}

library('SoupX', lib.loc=lib_path)
library('DropletUtils', lib.loc=lib_path)

set.seed(42)

# CLEAN DATA: outs #
#--------------------------------------------------------------------
soupify_outs <- function(in_path, out_path)
{
	print('loading')
	sc <- load10X(in_path)
	print('estimating')
	pdf(file = paste0(soupX_ouput_path, project, '_', sam, '_', 'contam_plot.pdf'), height = 8, width = 8)
	sc=autoEstCont(sc)
	dev.off()
	print('adjusting')
	out=adjustCounts(sc, roundToInt=TRUE)
	print('writing')
	write10xCounts(out_path, out, version='3', overwrite=TRUE)
}

if (data_type == 'outs')
{
	print('outs')

	cellranger_data = paste(soupX_input_path, 'outs/', sep='')
	print(cellranger_data)
	soupify_outs(cellranger_data, soupX_ouput_path)
}
#--------------------------------------------------------------------

# CLEAN DATA: no cluster information  #
#--------------------------------------------------------------------
soupify_noclusters <- function(in_path, out_path)
{
	library('Seurat', lib.loc=lib_path)

	print(in_path)
	print(out_path)
	r=paste(in_path, 'raw_feature_bc_matrix/', sep='')
	f=paste(in_path, 'filtered_feature_bc_matrix/', sep='')

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
	pdf(file = paste0(soupX_ouput_path, project, '_', sam, '_', 'contam_plot.pdf'), height = 8, width = 8)
	sc=autoEstCont(sc)
	dev.off()
	out=adjustCounts(sc, roundToInt=TRUE)
	write10xCounts(out_path, out, version='3', overwrite = TRUE)
}

# users need sub-directories (57,58), and no outs dir
if (data_type == 'no_clusters')
{
	soupify_noclusters(soupX_input_path, soupX_ouput_path)
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
	pdf(file = paste0(soupX_ouput_path, project, '_', sam, '_', 'contam_plot.pdf'), height = 8, width = 8)
	sc=autoEstCont(sc)
	dev.off()
	out = adjustCounts(sc, roundToInt=TRUE)
	write10xCounts(out_path, out, version='3', overwrite = TRUE)
}

# users need files (98,99), and no outs dir
if (data_type == 'h5') 
{
	print('h5')

	#cellranger_data = paste(soupX_input_path, 'outs/', sep='')
	#soupify_outs(cellranger_data, soupX_ouput_path)

	soupify_noclusters(soupX_input_path, soupX_ouput_path)
}
