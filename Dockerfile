ARG UBUNTU_VER=18.04
FROM ubuntu:${UBUNTU_VER}

LABEL authors="Gaofei Zhao, Jun Woo"

# Updates and install necessary packages
RUN apt-get update && \
    apt-get install -y build-essential && \
    apt-get install -y wget && \
    apt-get install -y libgsl-dev unzip python3 python3-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

ENV PATH=$CONDA_DIR/bin:$PATH

# Install R packages
RUN conda install -y -c r r-base r-essentials r-devtools tensorflow

# Install project packages
RUN ln -s /bin/tar /bin/gtar
RUN R -e 'options(unzip = "internal")'
RUN R -e 'devtools::install_github("mskcc/DeepSig")'

# Copy and set work directory
COPY . ./DeepSig
WORKDIR ./DeepSig

RUN mkdir input
RUN mkdir output
