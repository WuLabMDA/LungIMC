#!/usr/bin/env bash

BASEDIR="/rsrch1/ip/pchen6/LungIMCData/LungROIProcessing/Steinbock"
cd ${BASEDIR}

# setup steinbock alias
shopt -s expand_aliases
alias steinbock="docker run -v ${BASEDIR}:/data -u $(id -u):$(id -g) ghcr.io/bodenmillergroup/steinbock:0.14.1"

# measurement
steinbock measure intensities --masks masks_deepcell
steinbock measure regionprops --masks masks_deepcell
steinbock measure neighbors --masks masks_deepcell --type expansion --dmax 4