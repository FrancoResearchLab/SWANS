FROM rocker/tidyverse:4.4.0
LABEL maintainer="wafulae@chop.edu"
WORKDIR /rocker-build/

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive

### Install apt-getable packages to start
#########################################

# Installing all apt required packages at once
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
    build-essential \
    bzip2 \
    curl \
    jq \
    libgmp3-dev \
    libgdal-dev \
    libudunits2-dev \
    libpoppler-cpp-dev \
    libglpk-dev \
    libncurses5 \
    libssl-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libreadline-dev \
    libsqlite3-dev \
    libgdbm-dev \
    libdb5.3-dev \
    libbz2-dev \
    libexpat1-dev \
    liblzma-dev \
    libffi-dev \
    libuuid1 \
    gsl-bin \
    libgsl-dev \
    wget \
    xorg \
    zlib1g-dev \
    sendmail \
    sendmail-cf \
    m4 \
    mailutils \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Download and install Python 3.11
RUN cd /usr/src && \
    wget https://www.python.org/ftp/python/3.11.0/Python-3.11.0.tgz && \
    tar xzf Python-3.11.0.tgz && \
    cd Python-3.11.0 && \
    ./configure --enable-optimizations && \
    make altinstall && \
    rm -rf /usr/src/Python-3.11.0.tgz

# Setup the default python commands to use Python 3.11
RUN ln -s /usr/local/bin/python3.11 /usr/local/bin/python3 && \
    ln -s /usr/local/bin/python3.11 /usr/local/bin/python
RUN python3 -m pip install --upgrade pip

# Set working directory
WORKDIR /home/rstudio


# Install python packages
##########################

# Install python3 tools and ALL dependencies
RUN pip3 install \
    "matplotlib==3.7.1" \
    "numpy==1.24.3" \
    "pandas==2.0.1" \
    "scikit-learn==1.2.2" \
    "scipy==1.10.1" \
    "seaborn==0.12.2" \
    "multiqc==1.25" \
    "setuptools==46.3.0" \
    "snakemake==8.11.6" \
    "umap-learn==0.5.7" \
    "utils==1.0.1" \
    "wheel==0.34.2" \
    && rm -rf /root/.cache/pip/wheels

# Standalone tools and libraries
################################

# Add Cell Ranger - download link needs to be update with every build the key pair expires
RUN wget -O cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1749881778&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=aLWS0DZHXLNBHx14vWAvn52VJst8ddZ9rPcvowW5IIRLLlaYRIlF-33WFqH01uEUiJaOG5VAx6OlpyhJX8zlfT-KYskxxjjjEId2MgtOqIPP50syD9zCdLQn657UQVswrUkkYioHdt-yZ~wn0N5BNjbQR2OHQbLZqQsPvGou4Ix7YxbJn7T8ZcOjl0sXNDAdGorqio8Iici65SuMTVZ4dSv8Aj38abUMTtUBsPB3NcHslDhCitcniFcjug-Po0wZY8fYWtNzpU7Wm2Ju-CdVXv6g3bYf9bERiiVv4kXbh9xjp2NrfLZnGxbm4F0ZsGkhfQFrU21dk4yVwELagPpMWw__" && \
    tar -zxvf cellranger-9.0.1.tar.gz && rm -f cellranger-9.0.1.tar.gz && \
    mv cellranger-9.0.1 /usr/local/bin/ && \
    ln -s /usr/local/bin/cellranger-9.0.1/cellranger /usr/local/bin/cellranger

#### R packages
###############

# Set the Bioconductor repository as the primary repository
RUN R -e "options(repos = BiocManager::repositories())"

# Install BiocManager and the desired version of Bioconductor
RUN R -e "install.packages('BiocManager', dependencies=TRUE)"
RUN R -e "BiocManager::install(version = '3.20')"

# Install matrixStats 1.4.1
# dependency required MatrixGenerics CRAN has older version 1.3.0
RUN R -e "remotes::install_github('HenrikBengtsson/matrixStats', ref='develop')"

# Install R packages
RUN R -e 'BiocManager::install(c( \
  "multtest", \
  "DropletUtils", \
  "SingleCellExperiment", \
  "SummarizedExperiment", \
  "Biobase", \
  "GenomicRanges", \
  "GenomeInfoDb", \
  "GenomicAlignments", \
  "glmGamPoi", \
  "IRanges", \
  "S4Vectors", \
  "BiocGenerics", \
  "MatrixGenerics", \
  "DirichletMultinomial", \
  "TFBSTools", \
  "fgsea", \
  "limma", \
  "DelayedArray", \
  "DelayedMatrixStats", \
  "batchelor", \
  "HDF5Array", \
  "biomaRt", \
  "BSgenome.Hsapiens.UCSC.hg38", \
  "EnsDb.Hsapiens.v86", \
  "harmony", \
  "SeuratObject", \
  "sctransform", \
  "msigdbr", \
  "reshape2", \
  "RColorBrewer", \
  "progressr", \
  "presto", \
  "ggrepel", \
  "data.table", \
  "future", \
  "clustree", \
  "ggraph", \
  "yaml", \
  "pdftools", \
  "rsconnect", \
  "lme4", \
  "DT", \
  "xfun", \
  "sp", \
  "qs", \
  "terra", \
  "ggrastr", \
  "vroom", \
  "devtools", \
  "tidyverse", \
  "viridis", \
  "viridisLite", \
  "SoupX", \
  "Seurat", \
  "NMF", \
  "rcartocolor", \
  "gprofiler2", \
  "hdf5r", \
  "jsonlite", \
  "htmltools" \
  ))'

# Install Mononle3
RUN R -e "remotes::install_github('cole-trapnell-lab/monocle3')"

# Install DoubletFinder
RUN R -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')"

# Install Azimuth
RUN R -e "remotes::install_github('satijalab/seurat-data')"
RUN R -e "remotes::install_github('satijalab/azimuth', ref = 'master')"

# Install SeuratWrappers - large download times out
RUN R -e "options(timeout=9999999)"
RUN R -e "remotes::install_github('satijalab/seurat-wrappers')"

# Install CellChat and dependencies addition to NMF  
RUN R -e "remotes::install_github('jokergoo/circlize')"
RUN R -e "remotes::install_github('jokergoo/ComplexHeatmap')"
RUN R -e "remotes::install_github('jinworks/CellChat')"

# Install LoupeR from 10x genomics 
# Will require interactively running "loupeR::setup()" in the constainer 
# to accept the license to install finalize installing executables.
# Otherwise, loading the library, "library(loupeR)" will prompt use to 
# run "loupeR::setup()".
RUN R -e "remotes::install_github('10XGenomics/loupeR')"

# Reset the frontend variable for interactive
ENV DEBIAN_FRONTEND=

WORKDIR /rocker-build/