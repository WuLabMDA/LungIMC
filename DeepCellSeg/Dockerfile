FROM vanvalenlab/deepcell-tf:0.11.0-gpu
MAINTAINER pingjunchen <pingjunchen@ieee.org>

# Folder preparation
WORKDIR /.local
RUN chmod 777 /.local
WORKDIR /.keras
RUN chmod 777 /.keras
WORKDIR /Data
RUN chmod 777 /Data
WORKDIR /App
RUN chmod 777 /App

# Set working directory
WORKDIR /App/LungIMC/DeepCellSeg

CMD ["/bin/bash"]