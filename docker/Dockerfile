FROM ubuntu:22.04

SHELL ["/bin/bash", "-o", "pipefail", "-c"]


RUN apt-get update --fix-missing && \
  apt-get install -q -y wget curl bzip2 libbz2-dev git build-essential zlib1g-dev locales vim fontconfig 


# Set the locale
RUN locale-gen en_US.UTF-8  
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8     

# Install conda
# This version uses python version 3.7
# So we prefer to use it with ribopy.
ARG CONDA=Miniconda3-latest-Linux-x86_64.sh
RUN curl -LO http://repo.continuum.io/miniconda/$CONDA && \
    bash $CONDA -p /miniconda3 -b && \
    rm $CONDA 
ENV PATH /miniconda3/bin:${PATH}

# # Install conda dependencies
ADD environment.yaml /
ADD VERSION /
RUN pwd
RUN conda config --set always_yes yes --set changeps1 no && \
    conda config --add channels conda-forge && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --get && \
    conda update -q conda && \
    conda info -a && \
    conda env create -f environment.yaml  
 #   conda env update -n root --file environment.yaml && \

# Make the "ribo" environment the default.
ENV PATH /miniconda3/envs/ribo/bin:$PATH
RUN /bin/bash -c "source activate ribo"