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
# import umap

from pheno_utils import antibody_names, random_colors

def set_args():
    parser = argparse.ArgumentParser(description = "IMC Community Detection")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--batchcorrection_dir",    type=str,       default="BatchCorrection")
    parser.add_argument("--phenotype_dir",          type=str,       default="Phenotype")
    parser.add_argument("--control_option",         type=str,       default="NoControl", choices = ["NoControl", "WithControl"])
    parser.add_argument("--fea_option",             type=str,       default="Corrected", choices = ["Transformed", "Corrected"])
    parser.add_argument("--sample",                 default=True,  action="store_false")
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
    fea_path = os.path.join(cellfea_dir, args.fea_option + "Feas.RData")
    fea_rdata = pyreadr.read_r(fea_path)
    community_path = os.path.join(cellfea_dir, args.fea_option + "CommunitiesSOM.RData")
    community_rdata = pyreadr.read_r(community_path)
    cell_feas, communities = None, None
    if args.fea_option == "Corrected":
        cell_feas = fea_rdata["corrected"]
        communities = community_rdata["correct_communities"]
    else:
        cell_feas = fea_rdata["uncorrected"]
        communities = community_rdata["transform_communities"]
    # Obtain cell features
    cell_feas = cell_feas[antibody_names].to_numpy()
    # convert communites to list (integer)
    communities = communities.to_numpy().astype(int)
    communities = np.squeeze(communities)
    # Cell sampling
    if args.sample:
        cell_indices = np.arange(0, cell_feas.shape[0])
        np.random.shuffle(cell_indices)
        sample_indices = cell_indices[:args.sample_cell_size].astype(int)
        cell_feas = cell_feas[sample_indices, :]
        communities = communities[sample_indices]
    else:
        sample_indices = np.arange(0, cell_feas.shape[0])
    communities = communities.tolist()
    # Print cell information
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
    s_sne_name = "SOM_{}{}CommunitiesTSNE{}Cells{}Markers.png".format(args.fea_option, community_num, cell_num, fea_num)
    t_sne_path = os.path.join(phenotype_dir, s_sne_name)
    plt.savefig(t_sne_path, dpi=300)
    plt.close()
    print("t-SNE Done @ ", datetime.now().strftime("%H:%M:%S"))

    # Draw antibody heatmap
    fea_heatmap_dir = os.path.join(phenotype_dir, "SOM_{}{}CellsStainsHeatmap".format(args.fea_option, cell_num))
    if os.path.exists(fea_heatmap_dir):
        shutil.rmtree(fea_heatmap_dir)
    os.makedirs(fea_heatmap_dir)
    min_feas = np.percentile(cell_feas, q=5, axis=0)
    max_feas = np.percentile(cell_feas, q=95, axis=0)
    range_feas = max_feas - min_feas
    norm_feas = (cell_feas - min_feas) / range_feas
    norm_feas[norm_feas < 0] = 0
    norm_feas[norm_feas > 1] = 1.0
    for ind, fea_name in enumerate(antibody_names):
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(7, 5))
        axes.scatter(embed_feas[:, 0], embed_feas[:, 1], c=norm_feas[:, ind], s=0.1, cmap=plt.cm.jet)
        axes.set_title("Feature: {}".format(fea_name))
        plt.savefig(os.path.join(fea_heatmap_dir, "{}_heatmap.png".format(fea_name)), dpi=300)
        plt.close()