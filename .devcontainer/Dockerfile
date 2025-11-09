# Use specific R version
FROM rocker/tidyverse:4.0.5

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libgit2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    python3-dev \
    python3-pip \
    wget \
    curl \
    git \
    vim \
    htop \
    nano \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages for single-cell analysis
RUN pip3 install --upgrade pip
RUN pip3 install \
    scanpy==1.9.6 \
    scvelo==0.3.0 \
    cellrank==2.0.2 \
    scrublet==0.2.3 \
    anndata==0.9.2 \
    leidenalg==0.10.1 \
    umap-learn==0.5.4 \
    scikit-learn==1.0.2 \
    scipy==1.10.1 \
    pandas==1.5.3 \
    numpy==1.24.3 \
    jupyterlab==4.0.10 \
    jupyter==1.0.0 \
    radian==0.6.2 \
    plotly==5.17.0

# Install R packages for VS Code integration
RUN R -e "install.packages(c( \
    'languageserver', \
    'httpgd', \
    'lsp', \
    'devtools', \
    'BiocManager' \
    ), dependencies=TRUE)"

# Install core single-cell packages
RUN R -e "install.packages(c( \
    'tidyverse', \
    'Seurat', \
    'ggplot2', \
    'ggpubr', \
    'patchwork', \
    'matrix', \
    'clustree', \
    'harmony', \
    'DoubletFinder' \
    ), dependencies=TRUE)"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c( \
    'SingleCellExperiment', \
    'scater', \
    'scran', \
    'celldex', \
    'SingleR', \
    'ComplexHeatmap', \
    'limma', \
    'DESeq2', \
    'edgeR', \
    'EnsDb.Hsapiens.v86', \
    'org.Hs.eg.db' \
    ), ask=FALSE, update=FALSE)"

# Install GitHub packages
RUN R -e "devtools::install_github('satijalab/seurat-wrappers')"
RUN R -e "devtools::install_github('immunogenomics/harmony')"

# Install visualization packages
RUN R -e "install.packages(c( \
    'viridis', \
    'RColorBrewer', \
    'pheatmap', \
    'ggrepel', \
    'cowplot', \
    'ggridges', \
    'circlize', \
    'ggbeeswarm' \
    ), dependencies=TRUE)"

# Install Seurat explicitly (in case it wasn't installed above)
RUN R -e "if (!requireNamespace('Seurat', quietly = TRUE)) { \
    install.packages('Seurat', dependencies=TRUE) \
}"

# Create workspace
RUN mkdir -p /media/workspace/4f1d6914-c8de-4f37-b4bd-fd311d8ef82d1/SAARLAND UNIVERSITY/Msc course/WS 2025_26/single-cell/single-cell-analysis/{data,results,scripts,Notebooks,references,R}
WORKDIR /media/workspace/4f1d6914-c8de-4f37-b4bd-fd311d8ef82d1/SAARLAND UNIVERSITY/Msc course/WS 2025_26/single-cell/single-cell-analysis

# Set default command
CMD ["/bin/bash"]
