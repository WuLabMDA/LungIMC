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
    # prepare directory
    vis_dir = os.path.join(args.data_root, args.vis_dir, args.fea_option)
    if not os.path.exists(vis_dir):
        os.makedirs(vis_dir)

    # load metadata
    cellfea_dir = os.path.join(args.data_root, args.batchcorrection_dir, "RData")
    metadata_path = os.path.join(cellfea_dir, "StudyMeta.RData")
    metadata_dict = pyreadr.read_r(metadata_path)
    fea_filenames = metadata_dict["fea_filenames"]["fea_filenames"].to_list()
    # markers = metadata_dict["markers"]["markers"].to_list()
    antibody_names = metadata_dict["markers"]["markers"].to_list()
    roi_nrows = metadata_dict["roi_nrows"]["roi_nrows"].to_list()
    # load cell features and communites
    cell_feas, communities = None, None
    cellfea_dir = os.path.join(args.data_root, args.batchcorrection_dir, "RData")
    fea_path = os.path.join(cellfea_dir, args.fea_option + "Feas.csv")
    cell_feas = pd.read_csv(fea_path)
    cell_feas = cell_feas[antibody_names].to_numpy()
    community_path = os.path.join(cellfea_dir, args.fea_option + "CommunitiesSOM.RData")
    community_rdata = pyreadr.read_r(community_path)
    communities = community_rdata["communities"]
    communities = communities.to_numpy().astype(int)
    communities = np.squeeze(communities)

    # load ROI diagnosis
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
    # build cell source list
    cell_diag_list = []
    for ind, ele in enumerate(fea_filenames):
        cell_num = roi_nrows[ind]
        roi_diag = roi_diag_dict[ele]
        cell_diag_list.extend([roi_diag] * cell_num)

    # Cell sampling
    if args.sample:
        cell_indices = np.arange(0, cell_feas.shape[0]).astype(int)
        np.random.shuffle(cell_indices)
        sample_indices = cell_indices[:args.sample_cell_size]
        cell_feas = cell_feas[sample_indices, :]
        communities = communities[sample_indices]
        cell_diag_list = [cell_diag_list[ind] for ind in sample_indices]
    else:
        sample_indices = np.arange(0, cell_feas.shape[0]).astype(int)
    communities = communities.tolist()

    # print cell information
    cell_num, fea_num = cell_feas.shape
    print("Input data has {} cells and {} features.".format(cell_num, fea_num))
    unique_ids = np.unique(communities)
    community_num = len(unique_ids)
    print("Number of communities is: {}".format(community_num))
    # Draw t-SNE embbeded map
    print("Start t-SNE @ ", datetime.now().strftime("%H:%M:%S"))
    tsne = TSNE(n_components=2)
    embed_feas = tsne.fit_transform(cell_feas)
    community_colors = random_colors(community_num)
    color_dict = {unique_ids[ind]: (np.array(cur_color) * 255.0).astype(np.uint8) for ind, cur_color in enumerate(community_colors)}
    cell_colors = [color_dict[val] for val in communities]
    hex_colors = ["#{:02x}{:02x}{:02x}".format(ele[0], ele[1], ele[2]) for ele in cell_colors]
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(7, 5))
    axes.scatter(embed_feas[:, 0], embed_feas[:, 1], c=hex_colors, s=0.1)
    axes.set_title("Cell Community Detection")
    s_sne_name = "SOM{}CommunitiesTSNE{}Cells{}Markers.png".format(community_num, cell_num, fea_num)
    t_sne_path = os.path.join(vis_dir, s_sne_name)
    plt.savefig(t_sne_path, dpi=300)
    plt.close()
    print("t-SNE Done @ ", datetime.now().strftime("%H:%M:%S"))

    # feature normalization
    min_feas = np.percentile(cell_feas, q=1, axis=0)
    max_feas = np.percentile(cell_feas, q=99, axis=0)
    range_feas = max_feas - min_feas
    norm_feas = (cell_feas - min_feas) / range_feas
    norm_feas[norm_feas < 0] = 0
    norm_feas[norm_feas > 1] = 1.0

    # Draw antibody heatmap
    interested_antibodies = ["CK", "CD45"]
    for antibody_name in interested_antibodies:
        ck_embed_dir = os.path.join(vis_dir, "{}_SOM_Embed".format(antibody_name))
        if not os.path.exists(ck_embed_dir):
            os.makedirs(ck_embed_dir)
        antibody_ind = antibody_names.index(antibody_name)
        # draw overall figure
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(7, 5))
        axes.scatter(embed_feas[:, 0], embed_feas[:, 1], c=norm_feas[:, antibody_ind], s=0.1, cmap=plt.cm.jet)
        axes.set_title("Feature: {}".format(antibody_name))
        plt.savefig(os.path.join(ck_embed_dir, "{}_heatmap.png".format(antibody_name)), dpi=300)
        plt.close()
        # draw seperate figure
        diagnosis_list = list(set(cell_diag_list))
        for cur_diag in diagnosis_list:
            cell_idx = [ind for ind, ele in enumerate(cell_diag_list) if ele == cur_diag]
            fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(7, 5))
            axes.scatter(embed_feas[cell_idx, 0], embed_feas[cell_idx, 1], c=norm_feas[cell_idx, antibody_ind], s=0.1, cmap=plt.cm.jet)
            axes.set_title("Feature: {}-{}".format(antibody_name, cur_diag))
            plt.savefig(os.path.join(ck_embed_dir, "{}-{}_heatmap.png".format(antibody_name, cur_diag)), dpi=300)
            plt.close()
