# -*- coding: utf-8 -*-

import os, sys
import numpy as np
from datetime import datetime
import argparse, shutil, pickle, pytz
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pyreadr
from sklearn.manifold import TSNE

# from pheno_utils import antibody_names, random_colors
from pheno_utils import random_colors


def set_args():
    parser = argparse.ArgumentParser(description = "IMC Community Detection")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--batchcorrection_dir",    type=str,       default="BatchCorrection")
    parser.add_argument("--vis_dir",                type=str,       default="Vis")
    parser.add_argument("--fea_option",             type=str,       default="Transform", choices = ["Transform", "SelfCorrect", "ControlCorrect"])
    parser.add_argument("--sample",                 default=True,   action="store_false")
    parser.add_argument("--sample_cell_size",       type=int,       default=100000)
    parser.add_argument("--seed",                   type=int,       default=1234)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()
    # load data
    cellfea_dir = os.path.join(args.data_root, args.batchcorrection_dir, "RData")
    # load metadata
    metadata_path = os.path.join(cellfea_dir, "StudyMeta.RData")
    metadata_dict = pyreadr.read_r(metadata_path)
    fea_filenames = metadata_dict["fea_filenames"]["fea_filenames"].to_list()
    # markers = metadata_dict["markers"]["markers"].to_list()
    antibody_names = metadata_dict["markers"]["markers"].to_list()
    roi_nrows = metadata_dict["roi_nrows"]["roi_nrows"].to_list()

    # load ROI diagnosis information
    study_roi_info_path = os.path.join(args.data_root, "Metadata", "StudyROI_Info.xlsx")
    slide_roi_info = pd.read_excel(study_roi_info_path, sheet_name = "Sheet1", header= 0, index_col=None)
    name_loc_dict = {key: value for (key, value) in zip(slide_roi_info["ROI_ID"], slide_roi_info["ROI_Location"])}
    name_diag_dict = {key: value for (key, value) in zip(slide_roi_info["ROI_ID"], slide_roi_info["ROI_Diag"])}
    roi_diag_dict = {}
    for ele in fea_filenames:
        roi_loc = name_loc_dict[ele]
        if roi_loc in ["AdjacentNormal", "DistantNormal", "Normal"]:
            roi_diag_dict[ele] = "Normal"
        elif roi_loc == "Tumor":
            roi_diag_dict[ele] = name_diag_dict[ele]
        else:
            print("Unknow location: {}".format(roi_loc))
    import pdb; pdb.set_trace()
    a = 1

    # load cell features and communites
    cell_feas, communities = None, None
    cellfea_dir = os.path.join(args.data_root, args.batchcorrection_dir, "RData")
    # obtain cell features
    fea_path = os.path.join(cellfea_dir, args.fea_option + "Feas.csv")
    cell_feas = pd.read_csv(fea_path)
    cell_feas = cell_feas[antibody_names].to_numpy()
    # load communites
    community_path = os.path.join(cellfea_dir, args.fea_option + "CommunitiesSOM.RData")
    community_rdata = pyreadr.read_r(community_path)
    communities = community_rdata["communities"]
    communities = communities.to_numpy().astype(int)
    communities = np.squeeze(communities)
