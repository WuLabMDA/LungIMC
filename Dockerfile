FROM nvidia/cuda:11.3.0-base-ubuntu20.04
MAINTAINER pingjunchen <pingjunchen@ieee.org>

# Setting environment variables
ENV DEBIAN_FRONTEND noninteractive
ENV HDF5_USE_FILE_LOCKING FALSE

RUN apt-get update && apt-get install -y --no-install-recommends \
  build-essential libgl1-mesa-glx libglib2.0-0 libgeos-dev libvips-tools \
  curl sudo git wget vim ca-certificates python3-openslide \
  && rm -rf /var/lib/apt/lists/*

# Folder preparation
WORKDIR /Data
RUN chmod 777 /Data
WORKDIR /App
RUN chmod 777 /App

# Install Miniconda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-py38_4.8.2-Linux-x86_64.sh \
 && bash Miniconda3-py38_4.8.2-Linux-x86_64.sh -p /App/miniconda -b \
 && rm Miniconda3-py38_4.8.2-Linux-x86_64.sh
ENV PATH=/App/miniconda/bin:$PATH
## Create a Python 3.8.3 environment
RUN /App/miniconda/bin/conda install conda-build \
 && /App/miniconda/bin/conda create -y --name py383 python=3.8.3 \
 && /App/miniconda/bin/conda clean -ya
ENV CONDA_DEFAULT_ENV=py383
ENV CONDA_PREFIX=/App/miniconda/envs/$CONDA_DEFAULT_ENV
ENV PATH=$CONDA_PREFIX/bin:$PATH

# Install python packages
RUN pip install gpustat==0.6.0 setuptools==45 pytz==2021.1
RUN pip install deepdish==0.3.6 pandas==1.1.5 seaborn==0.11.2 matplotlib==3.5.0
RUN pip install scikit-learn==1.0.1 xgboost==1.5.1 statsmodels==0.13.1
RUN pip install opencv-python==4.5.4.60 scikit-image==0.18.0
RUN pip install pyreadr==0.4.4

# Set working directory
WORKDIR /App/LungIMC
