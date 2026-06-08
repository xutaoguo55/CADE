FROM rocker/r-ver:4.5.3

LABEL maintainer="Xutao Guo <gxt827@126.com>" \
      description="CADE: Composition-Aware Differential Expression — reproducible environment" \
      version="1.1.0"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
    libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
    libtiff5-dev libjpeg-dev libgsl-dev libnetcdf-dev \
    python3 python3-pip python3-numpy python3-pandas python3-matplotlib \
    && rm -rf /var/lib/apt/lists/*

# R environment
RUN mkdir -p /usr/local/lib/R/etc/
RUN echo "options(repos = c(CRAN = 'https://cloud.r-project.org'), \
    BioC_mirror = 'https://bioconductor.org')" > /usr/local/lib/R/etc/Rprofile.site

# Install BiocManager
RUN Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")'

# Core packages
RUN Rscript -e 'BiocManager::install(c( \
    "limma", "GSVA", "WGCNA", "GEOquery", "sva", "RUVSeq", \
    "pROC", "quadprog", "EPIC", "fgsea", "msigdbr", \
    "celldex", "dplyr", "tibble", "ggplot2", "patchwork", \
    "pheatmap", "gplots", "affy", "annotate", \
    "hgu133plus2.db", "GSEABase", "MASS", \
    update=FALSE, ask=FALSE)'

# Install from CRAN
RUN Rscript -e 'install.packages(c("cowplot", "ggrepel", "scales", \
    "readxl", "openxlsx"), repos="https://cloud.r-project.org")'

# Python packages
RUN pip3 install --break-system-packages numpy pandas matplotlib scipy scikit-learn

# Create working directory
WORKDIR /cade
COPY . /cade/

# Verify installation
RUN Rscript -e 'cat("=== CADE Environment ===\n"); \
    cat("R version:", as.character(getRversion()), "\n"); \
    for(p in c("limma","GSVA","WGCNA","GEOquery","sva","RUVSeq","pROC","quadprog","EPIC")) { \
        cat(sprintf("  %s: %s\n", p, as.character(packageVersion(p)))) }; \
    cat("========================\n")'

# Default command
CMD ["Rscript", "code/run_all.R", "--skip-benchmarks"]
