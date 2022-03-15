# -*- coding: utf-8 -*-

import os, sys
import numpy as np
from datetime import datetime
import argparse, shutil, pickle, pytz
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pyreadr


def set_args():
    parser = argparse.ArgumentParser(description = "IMC Community Detection")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--batchcorrection_dir",    type=str,       default="BatchCorrection")
    parser.add_argument("--phenotype_dir",          type=str,       default="Phenotype")
    parser.add_argument("--control_option",         type=str,       default="NoControl", choices = ["NoControl", "WithControl"])
    parser.add_argument("--fea_option",             type=str,       default="Corrected", choices = ["Transformed", "Corrected"])
    parser.add_argument("--sample",                 default=False,  action="store_true")
    parser.add_argument("--sample_cell_size",       type=int,       default=100000)
    parser.add_argument("--seed",                   type=int,       default=1234)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()
    # prepare directory
    phenotype_dir = os.path.join(args.data_root, args.phenotype_dir, args.control_option)
    if not os.path.exists(phenotype_dir):
        os.makedirs(phenotype_dir)
    cellfea_dir = os.path.join(args.data_root, args.batchcorrection_dir, args.control_option)
    community_path = os.path.join(cellfea_dir, args.fea_option + "CommunitiesSOM.RData")
    community_rdata = pyreadr.read_r(community_path)
    cell_feas, communities = None, None
    if args.fea_option == "Corrected":
        cell_feas = community_rdata["cell_feas"]
        communities = community_rdata["correct_communities"]
    else:
        cell_feas = community_rdata["cell_feas"]
        communities = community_rdata["transform_communities"]
    # Obtain cell features
    antibody_names = list(cell_feas)
    cell_feas = cell_feas.to_numpy()
    # convert communites to list (integer)
    communities = communities.to_numpy().astype(int)
    communities = np.squeeze(communities).tolist()
    # Cell sampling
    if args.sample:
        cell_indices = np.arange(0, cell_feas.shape[0])
        np.random.shuffle(cell_indices)
        sample_indices = cell_indices[:args.sample_cell_size].astype(int)
        cell_feas = cell_feas[sample_indices, :]
        communities = communities[sample_indices]
    else:
        sample_indices = np.arange(0, cell_feas.shape[0])
    # Print cell information
    cell_num, fea_num = cell_feas.shape
    print("Input data has {} cells and {} features.".format(cell_num, fea_num))
    unique_ids = np.unique(communities)
    community_num = len(unique_ids)
    print("Number of communities is: {}".format(community_num))

    # Draw clustered staining heatmap
    heat_mat = np.zeros((community_num, fea_num), dtype=np.float32)
    cluster_ids = []
    for community_id in unique_ids:
        cur_community_fea = cell_feas[communities==community_id, :]
        cur_mean_fea = np.mean(cur_community_fea, axis=0)
        heat_mat[community_id-1, :] = cur_mean_fea
        cluster_ids.append(str(community_id) + "-" + "{:.3f}".format(np.sum(communities==community_id) * 1.0 / len(communities)))
    min_antibody = np.min(heat_mat, axis=0)
    max_antibody = np.max(heat_mat, axis=0)
    range_antibody = max_antibody - min_antibody
    norm_heat = (heat_mat - min_antibody) / range_antibody
    heat_df = pd.DataFrame(norm_heat, columns=antibody_names)
    heat_df.index = cluster_ids
    heat_g = sns.clustermap(data=heat_df, figsize=(5, 15), metric="euclidean", method="ward", col_cluster=False, cmap="jet")
    heatmap_name = "SOM_{}{}CommunitiesHeatmap{}Cells{}Markers.png".format(args.fea_option, community_num, cell_num, fea_num)
    fea_heatmap = os.path.join(phenotype_dir, heatmap_name)
    heat_g.savefig(fea_heatmap, dpi=300)
