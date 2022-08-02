FROM rocker/r-ver:4.2.0
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
  build-essential curl sudo git wget vim ca-certificates \
  libglib2.0-0 libgeos-dev libxml2 libxslt-dev libglpk-dev \
  libmagick++-6.q16-dev libproj-dev libudunits2-dev \
  libpq-dev gdal-bin libgdal-dev

# Folder preparation
WORKDIR /Data
RUN chmod 777 /Data
WORKDIR /App
RUN chmod 777 /App

# Install packages
RUN R -e "install.packages(c('BiocManager', 'XML', 'proj4'))"
RUN R -e "BiocManager::install(c('pheatmap', 'viridis', \
                       'zoo', 'BiocManager', 'devtools', 'tiff', \
                       'distill', 'openxlsx', 'ggrepel', 'patchwork', 'mclust', \
                       'RColorBrewer', 'uwot', 'Rtsne', 'harmony',  \
                       'SeuratObject', 'cowplot', 'kohonen', 'caret', \
                       'randomForest', 'ggridges', 'gridGraphics', \
                       'scales', 'tidyverse', 'lisaClust', 'spicyR'))"
RUN R -e "BiocManager::install(c('ncdfFlow', 'ggcyto', \
    'scuttle', 'scater', 'dittoSeq', 'batchelor', 'bluster','scran', 'Seurat'))"
RUN R -e "BiocManager::install(c('CytoML', 'FlowSOM', 'CATALYST'))"
RUN R -e "BiocManager::install(c('EBImage', 'SpatialExperiment', 'imcRtools', 'cytomapper'))"
RUN R -e "devtools::install_github('i-cyto/Rphenograph')"

# Set working directory
WORKDIR /App/MouseIMC

CMD ["/bin/bash"]
