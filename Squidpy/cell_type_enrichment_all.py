# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import scanpy as sc
import anndata
import squidpy as sq


def set_args():
    parser = argparse.ArgumentParser(description = "Load data into squidpy image container")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--meta_dir",               type=str,       default="Metadata")    
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--steinbock_dir",          type=str,       default="Steinbock")
    parser.add_argument("--result_dir",             type=str,       default="Results")
    parser.add_argument("--plot_format",            type=str,       default=".png", choices=[".png", ".pdf"])        

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    # directory setting
    steinbock_root_dir = os.path.join(args.data_root, args.data_type, args.steinbock_dir)
    metadata_dir = os.path.join(args.data_root, args.meta_dir)
    enrich_result_dir = os.path.join(args.data_root, args.result_dir, "Enrichment")
    if not os.path.exists(enrich_result_dir):
        os.makedirs(enrich_result_dir)

    # load steinbock exported data
    cur_time_str = datetime.now(pytz.timezone('America/Chicago')).strftime("%m/%d/%Y, %H:%M:%S")
    print("----Load steinbock data @ {}".format(cur_time_str))
    rois_h5ad_path = os.path.join(steinbock_root_dir, "rois.h5ad")
    rois_adata = sc.read_h5ad(rois_h5ad_path)
    # add spatial coordinates
    spatial_coords = np.transpose(np.asarray([rois_adata.obs["centroid-1"].tolist(), rois_adata.obs["centroid-0"].tolist()]))
    rois_adata.obsm["spatial"] = spatial_coords
    # obtain roi list based on h5ad file
    roi_list = []
    for cur_roi in rois_adata.obs["image"]:
        cur_roi_name = cur_roi[:-5]
        if cur_roi_name not in roi_list:
            roi_list.append(cur_roi_name)

    # load cell phenotype & library ids information
    cur_time_str = datetime.now(pytz.timezone('America/Chicago')).strftime("%m/%d/%Y, %H:%M:%S")
    print("----Load cell phenotypes @ {}".format(cur_time_str))
    cell_phenotype_dict = None
    cell_phenotype_path = os.path.join(steinbock_root_dir, "cell_phenotypes.json")
    with open(cell_phenotype_path) as fp:
        cell_phenotype_dict = json.load(fp)
    cell_ids = [ele for ele in cell_phenotype_dict.keys()]
    cell_phenotypes = []
    library_ids = []
    cell_nums = []
    # collect roi list
    roi_lst = []
    for cell_id in cell_ids:
        roi_name = cell_id[:cell_id.rfind("_")]
        if roi_name not in roi_lst:
            roi_lst.append(roi_name)
    # collect roi cell numbers
    roi_cell_dict = {}
    for cell_id in cell_ids:
        roi_name = cell_id[:cell_id.rfind("_")]
        if roi_name not in roi_cell_dict:
            roi_cell_dict[roi_name] = 1
        else:   
            roi_cell_dict[roi_name] += 1
    # collect cell_phenotypes
    for roi_name in roi_lst:
        cell_num = roi_cell_dict[roi_name]
        cell_phenotypes.extend([cell_phenotype_dict[roi_name + "_" + str(num)] for num in range(1, cell_num + 1)])
        library_ids.extend([roi_name, ] * cell_num)
        cell_nums.extend([ele for ele in range(1, cell_num + 1)])
    # convert to Categorical
    rois_cell_phenotypes = pd.Categorical(cell_phenotypes) 
    rois_library_ids = pd.Categorical(library_ids) 
    rois_adata.obs["cell_type"] = rois_cell_phenotypes
    rois_adata.obs["library_id"] = rois_library_ids
    rois_adata.obs["cell_id"] = cell_nums

    # spaital cell type enrichment analysis
    cur_time_str = datetime.now(pytz.timezone('America/Chicago')).strftime("%m/%d/%Y, %H:%M:%S")
    print("----Spaital cell type enrichment analysis @ {}".format(cur_time_str))
    sq.gr.spatial_neighbors(rois_adata, spatial_key="spatial", library_key="library_id", coord_type="generic", radius=16.0)
    sq.gr.nhood_enrichment(rois_adata, cluster_key="cell_type")
    sq.pl.nhood_enrichment(rois_adata, cluster_key="cell_type")    
    plot_name = "cell_type_enrichment_all"
    plot_path = os.path.join(enrich_result_dir, plot_name + args.plot_format)
    plt.savefig(plot_path, transparent=False, dpi=300, bbox_inches="tight")  
    cur_time_str = datetime.now(pytz.timezone('America/Chicago')).strftime("%m/%d/%Y, %H:%M:%S")
    print("----Finish @ {}".format(cur_time_str))