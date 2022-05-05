# -*- coding: utf-8 -*-

import os, sys
import numpy as np
from datetime import datetime
import argparse, shutil, pickle, pytz
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pyreadr

from phenotype_utils import interested_immune_antibodies

def set_args():
    parser = argparse.ArgumentParser(description = "IMC Community Detection")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--batchcorrection_dir",    type=str,       default="BatchCorrection")
    parser.add_argument("--phenotype_dir",          type=str,       default="Phenotype")
    parser.add_argument("--fea_option",             type=str,       default="Raw", choices = ["Raw", "Correct"])
    parser.add_argument("--seed",                   type=int,       default=1234)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()
    # prepare directory
    phenotype_dir = os.path.join(args.data_root, args.phenotype_dir, args.fea_option)
    if not os.path.exists(phenotype_dir):
        os.makedirs(phenotype_dir)
    cell_feas, communities = None, None

    # load communites
    cellfea_dir = os.path.join(args.data_root, args.batchcorrection_dir, "RData")
    community_path = os.path.join(cellfea_dir, args.fea_option + "FeaCommunities.RData")
    community_rdata = pyreadr.read_r(community_path)
    cell_feas = community_rdata["cell_feas"]
    communities = community_rdata["communities"]

    # Obtain cell features
    antibody_names = list(cell_feas)
    cell_feas = cell_feas.to_numpy()
    # reorder antibodies
    antibody_orders = [antibody_names.index(cur_antibody) for cur_antibody in interested_immune_antibodies]
    cell_feas = cell_feas[:, antibody_orders]
    # convert communites to list (integer)
    communities = communities.to_numpy().astype(int)
    communities = np.squeeze(communities).tolist()
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
    heat_df = pd.DataFrame(norm_heat, columns=interested_immune_antibodies)
    heat_df.index = cluster_ids
    # heat_g = sns.clustermap(data=heat_df, figsize=(20, 25), metric="euclidean", method="ward", col_cluster=False, cmap="jet", annot=True, fmt='.3g')
    heat_g = sns.clustermap(data=heat_df, figsize=(12, 28), metric="euclidean", method="ward", col_cluster=False, cmap="jet")
    # heat_g.cax.set_visible(False)
    heatmap_name = "Heatmap{}Cells{}Communities{}Markers.png".format(cell_num, community_num, fea_num)
    fea_heatmap = os.path.join(phenotype_dir, heatmap_name)
    heat_g.savefig(fea_heatmap, dpi=300)
