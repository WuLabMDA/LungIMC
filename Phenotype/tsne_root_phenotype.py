# -*- coding: utf-8 -*-

import os, sys
import random
import numpy as np
from datetime import datetime
import argparse, shutil, pickle, pytz
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pyreadr
from sklearn.manifold import TSNE

from phenotype_utils import antibody_names, random_colors
from phenotype_utils import interested_immune_antibodies, interested_nonimmnue_antibodies


def set_args():
    parser = argparse.ArgumentParser(description = "IMC Community Detection")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--batchcorrection_dir",    type=str,       default="BatchCorrection")
    parser.add_argument("--phenotype_dir",          type=str,       default="Phenotype")
    parser.add_argument("--fea_option",             type=str,       default="Transform", choices = ["Transform", "SelfCorrect"])
    parser.add_argument("--sample",                 default=True,   action="store_false")
    parser.add_argument("--sample_cell_size",       type=int,       default=100000)
    parser.add_argument("--seed",                   type=int,       default=1234)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()
    random.seed(args.seed)
    np.random.seed(args.seed)
    print("Start @ {}".format(datetime.now(pytz.timezone('America/Chicago')).strftime("%m/%d/%Y, %H:%M:%S")))
    # prepare directory
    phenotype_dir = os.path.join(args.data_root, args.phenotype_dir, args.fea_option)
    if not os.path.exists(phenotype_dir):
        os.makedirs(phenotype_dir)
    # cluster info
    identified_cluster_path = os.path.join(args.data_root, args.phenotype_dir, "IdentifiedClusters.xlsx")
    cluster_info = pd.read_excel(identified_cluster_path, sheet_name=args.fea_option, header=0, index_col=None)
    cluster_dict = {cluster_id: category for cluster_id, category in zip(cluster_info["ClusterID"], cluster_info["Category"])}
    # cluster_dict = {cluster_id: category for cluster_id, category in zip(cluster_info["ClusterID"], cluster_info["VisConfirm"])}

    # features & communites
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

    # obtain tsne embedded features
    embed_feas = None
    embed_feas_path = os.path.join(phenotype_dir, "FeaEmbedTSNE{}Cells.npy".format(cell_num))
    if os.path.exists(embed_feas_path):
        embed_feas = np.load(embed_feas_path)
    else:
        tsne = TSNE(n_components=2)
        embed_feas = tsne.fit_transform(cell_feas)
        np.save(embed_feas_path, embed_feas)

    # # Merge tsne Plotting
    # fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(7, 5))
    # identified_communites = [cluster_dict[val] for val in communities]
    # merge_unique_ids = ["CK", "CD45", "Stroma"]
    # community_colors = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
    # color_dict = {merge_unique_ids[ind]: (np.array(cur_color) * 255.0).astype(np.uint8) for ind, cur_color in enumerate(community_colors)}
    # cell_colors = [color_dict[val] for val in identified_communites]
    # hex_colors = ["#{:02x}{:02x}{:02x}".format(ele[0], ele[1], ele[2]) for ele in cell_colors]
    # axes.scatter(embed_feas[:, 0], embed_feas[:, 1], c=hex_colors, s=0.1)
    # axes.set_title("Merge Cell Distribution")
    # axes.set_xlim([-50, 50])
    # axes.set_ylim([-50, 50])
    # s_sne_name = "TSNE{}CellsMergeCommunities.png".format(cell_num)
    # t_sne_path = os.path.join(phenotype_dir, s_sne_name)
    # plt.savefig(t_sne_path, dpi=300)
    # plt.close()

    # # nterested Clusters
    # batch_num, cluster_num = 8, 8
    # interested_clusters = list(range(1+(batch_num-1)*cluster_num, 1+batch_num*cluster_num))
    # community_colors = random_colors(len(interested_clusters))
    # color_dict = {}
    # for ind, cur_color in enumerate(community_colors):
    #     ind_color = (np.array(cur_color) * 255.0).astype(np.uint8)
    #     color_dict[interested_clusters[ind]] = "#{:02x}{:02x}{:02x}".format(ind_color[0], ind_color[1], ind_color[2])
    # fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(7, 5))
    # for cur_cluster in interested_clusters:
    #     interested_inds = [i for i, val in enumerate(communities) if val == cur_cluster]
    #     intereted_feas = embed_feas[interested_inds]
    #     intereted_communites = [communities[ind] for ind in interested_inds]
    #     axes.scatter(intereted_feas[:, 0], intereted_feas[:, 1], c=color_dict[cur_cluster], label=str(cur_cluster), s=0.1)
    # axes.set_title("Distribution of Interested Clusters")
    # axes.legend(loc="upper left", markerscale=5, scatterpoints=1, fontsize=10)
    # axes.set_xlim([-60, 50])
    # axes.set_ylim([-50, 50])
    # s_sne_name = "TSNE{}CellsInterestedCommunitiesBatch{}.png".format(cell_num, batch_num)
    # t_sne_path = os.path.join(phenotype_dir, s_sne_name)
    # plt.savefig(t_sne_path, dpi=300)
    # plt.close()

    # normalize features
    min_feas = np.percentile(cell_feas, q=1, axis=0)
    max_feas = np.percentile(cell_feas, q=99, axis=0)
    range_feas = max_feas - min_feas
    norm_feas = (cell_feas - min_feas) / range_feas
    norm_feas[norm_feas < 0] = 0
    norm_feas[norm_feas < 0.1] = 0 # ignore less than 0.1
    norm_feas[norm_feas > 1] = 1.0
    # Draw heatmaps of interested antibodies
    interested_antibody_dict = {"Immnue": interested_immune_antibodies, "NonImmnue": interested_nonimmnue_antibodies}
    for category, interested_antibodies in interested_antibody_dict.items():
        fea_heatmap_dir = os.path.join(phenotype_dir, "{}-StainsHeatmap".format(category))
        if os.path.exists(fea_heatmap_dir):
            shutil.rmtree(fea_heatmap_dir)
        os.makedirs(fea_heatmap_dir)
        intereted_antibodies = []
        for fea_name in interested_antibodies:
            fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(7, 5))
            ind = antibody_names.index(fea_name)
            axes.scatter(embed_feas[:, 0], embed_feas[:, 1], c=norm_feas[:, ind], s=0.1, cmap=plt.cm.jet)
            axes.set_title("Feature: {}".format(fea_name))
            axes.set_xlim([-50, 50])
            axes.set_ylim([-50, 50])
            plt.savefig(os.path.join(fea_heatmap_dir, "{}_heatmap.png".format(fea_name)), dpi=300)
            plt.close()
        # draw merge features
        indices = [antibody_names.index(fea_name) for fea_name in interested_antibodies]
        max_feas = np.max(norm_feas[:, indices], axis=1)
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(7, 5))
        axes.scatter(embed_feas[:, 0], embed_feas[:, 1], c=max_feas, s=0.1, cmap=plt.cm.jet)
        axes.set_title("Merge Feature {}".format(category))
        axes.set_xlim([-50, 50])
        axes.set_ylim([-50, 50])
        plt.savefig(os.path.join(fea_heatmap_dir, "{}_heatmap.png".format(category)), dpi=300)
        plt.close()
    print("Finish @ {}".format(datetime.now(pytz.timezone('America/Chicago')).strftime("%m/%d/%Y, %H:%M:%S")))
