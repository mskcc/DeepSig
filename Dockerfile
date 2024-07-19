# To run DeepSig, edit the last command; also need to have ./DeepSig_0.9.8.tar.gz and ./.DeepSig directory (contains the models to run)
FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

ENV CPATH=/usr/include

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        curl libcurl4-openssl-dev \
        gfortran \
        lbzip2 libbz2-dev\
        libblas-dev liblapack-dev \
        liblzma-dev \
        xz-utils \
        openssl \
        libssl-dev \
        libxml2-dev \
        pkgconf \
        python3-dev python3-pip python3-venv python3-virtualenv \
        libgsl-dev \
        libpng-dev \
        libhdf5-dev \
        r-base \
    && rm -rf /var/lib/apt/lists/*

RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"), download.file.method = "libcurl")' >> $(Rscript -e 'cat(R.home())')/etc/Rprofile.site \
    && echo 'options(Ncpus = parallel::detectCores())' >> $(Rscript -e 'cat(R.home())')/etc/Rprofile.site \
    && Rscript -e 'options(warn=2); install.packages("remotes")'
RUN Rscript -e 'remotes::install_version("RColorBrewer", "1.1-3")' \
    && Rscript -e 'remotes::install_version("RCurl", "1.98-1.14")' \
    && Rscript -e 'remotes::install_version("Rcpp", "1.0.12")' \
    && Rscript -e 'remotes::install_version("RcppEigen", "0.3.4.0.0")' \
    && Rscript -e 'remotes::install_version("XML", "3.99-0.16.1")' \
    && Rscript -e 'remotes::install_version("argparse", "2.2.3")' \
    && Rscript -e 'remotes::install_version("clue", "0.3-65")' \
    && Rscript -e 'remotes::install_version("coneproj")' \
    && Rscript -e 'remotes::install_version("cowplot", "1.1.3")' \
    && Rscript -e 'remotes::install_version("data.table", "1.15.4")' \
    && Rscript -e 'remotes::install_version("ggplot2", "3.5.1")' \
    && Rscript -e 'remotes::install_version("ggrepel", "0.9.5")' \
    && Rscript -e 'remotes::install_version("gridGraphics", "0.5-1")' \
    && Rscript -e 'remotes::install_version("gtools", "3.9.5")' \
    && Rscript -e 'remotes::install_version("reticulate", "1.38.0")' \
    && Rscript -e 'remotes::install_version("scales", "1.3.0")' \
    && Rscript -e 'remotes::install_version("tidyr", "1.3.1")' \
    && Rscript -e 'remotes::install_version("httr2")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "BSgenome", \
        "BSgenome.Hsapiens.UCSC.hg19", \
        "Biostrings", \
        "S4Vectors", \
        "GenomeInfoDb", \
        "GenomicAlignments", \
        "GenomicRanges", \
        "Rhtslib", \
        "Rsamtools", \
        "SummarizedExperiment", \
        "rtracklayer" \
    ))' 
RUN Rscript -e 'remotes::install_version("restfulr", "0.0.15")'

COPY . /
RUN python3 -m venv .venv && ./.venv/bin/python -m pip install pandas h5py tensorflow==2.15.0 keras==2.15.0
RUN /bin/bash -c 'source /.venv/bin/activate' && Rscript -e  \
  'install.packages(pkgs = "/DeepSig_0.9.8.tar.gz", repos = NULL);  \
   library(DeepSig);  \
   data <- read.table(system.file("extdata", "tcga-brca_catalog.txt",package="DeepSig"));  \
   z <- DL.call(catalog = t(data), cancer.type = "breast")'
