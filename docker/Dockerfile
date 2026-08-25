FROM ubuntu:24.04

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update --fix-missing && \
    apt-get install -q -y wget curl bzip2 libbz2-dev git build-essential zlib1g-dev locales vim fontconfig

# Set the locale
RUN locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

# Install conda
ARG CONDA=Miniconda3-latest-Linux-x86_64.sh
RUN curl -LO https://repo.anaconda.com/miniconda/$CONDA && \
    bash $CONDA -p /miniconda3 -b && \
    rm $CONDA
ENV PATH=/miniconda3/bin:${PATH}

# Configure conda: zero out default_channels (avoids Anaconda TOS enforcement
# in conda 25+) and use only bioconda + conda-forge. libmamba is the default
# solver in conda 24+, so no separate install needed.
RUN cat > /root/.condarc <<'EOF'
channels:
  - conda-forge
  - bioconda
default_channels: []
channel_priority: flexible
always_yes: true
changeps1: false
EOF

# Install conda dependencies
ADD environment.yaml /
ADD VERSION /

RUN conda env create -f environment.yaml && \
    conda clean -afy

# Make the "riboflow" environment the default. The env provides
# Nextflow, every bioinformatics tool, AND umicollapse (bioconda) on a single
# openjdk>=17 — no separate Java 11 runtime or hand-shipped umicollapse.jar.
ENV PATH=/miniconda3/envs/riboflow/bin:${PATH}
RUN echo "source activate riboflow" >> ~/.bashrc
