FROM ubuntu:22.04
WORKDIR /WolbachiaSubmission
RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -y \
         locales \
         libtbb2 \
         unzip \
         wget \
         git \
         snakemake \
         curl \
    && rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR=/opt/conda
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" \
    && bash Miniforge3-$(uname)-$(uname -m).sh -bfp ${CONDA_DIR} \
    && rm -rf Miniforge3-$(uname)-$(uname -m).sh

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

# Install mamba
RUN conda install -n base -c conda-forge mamba=1

WORKDIR /WolbachiaSubmission
COPY . .

RUN snakemake -p --cores 1 --use-conda --conda-prefix /opt/conda/envs --conda-create-envs-only -s /WolbachiaSubmission/create-envs.smk
