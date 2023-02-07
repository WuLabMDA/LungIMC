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
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--meta_dir",               type=str,       default="Metadata")
    parser.add_argument("--steinbock_dir",          type=str,       default="Steinbock")
    parser.add_argument("--enrichment_dir",         type=str,       default="Enrichment")
    parser.add_argument("--result_dir",             type=str,       default="Results")
    parser.add_argument("--enrichment_mode",        type=str,       default="zscore", choices=["zscore", "count"])      
    parser.add_argument("--plot_format",            type=str,       default=".png", choices=[".png", ".pdf"])    
    parser.add_argument("--path_stage",             type=str,       default="Normal", choices=["Normal", "AAH", "AIS", "MIA", "ADC"])      

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    # directory setting
    roi_root_dir = os.path.join(args.data_root, args.data_type)
    metadata_dir = os.path.join(args.data_root, args.meta_dir)
    enrichment_dir = os.path.join(roi_root_dir, args.enrichment_dir, args.path_stage)
    enrich_result_dir = os.path.join(args.data_root, args.result_dir, "Enrichment")
    if not os.path.exists(enrich_result_dir):
        os.makedirs(enrich_result_dir)

    # load roi stage information
    roi_info_path = os.path.join(metadata_dir, "StudyROI_Info.xlsx")
    study_roi_df = pd.read_excel(roi_info_path)
    study_roi_lst, roi_diag_lst = study_roi_df["ROI_ID"].tolist(), study_roi_df["ROI_Diag"].tolist()
    enrichment_roi_lst = [roi for roi, diag in zip(study_roi_lst, roi_diag_lst) if diag == args.path_stage]
    # print("{} has {} ROIs.".format(args.path_stage, len(enrichment_roi_lst)))


    # load steinbock exported data
    cur_time_str = datetime.now(pytz.timezone('America/Chicago')).strftime("%m/%d/%Y, %H:%M:%S")
    print("----Load steinbock data @ {}".format(cur_time_str))
    rois_h5ad_path = os.path.join(enrichment_dir, "rois.h5ad")
    rois_adata = sc.read_h5ad(rois_h5ad_path)
    # add spatial coordinates
    spatial_coords = np.transpose(np.asarray([rois_adata.obs["centroid-1"].tolist(), rois_adata.obs["centroid-0"].tolist()]))
    rois_adata.obsm["spatial"] = spatial_coords
    # obtain roi list based on h5ad file
    roi_lst = []
    for cur_roi in rois_adata.obs["image"]:
        cur_roi_name = cur_roi[:-5]
        if cur_roi_name not in roi_lst:
            roi_lst.append(cur_roi_name)
    assert len(roi_lst) == len(enrichment_roi_lst), "{} ROI number not match.".format(args.path_stage)
    # load cell phenotype & library ids information
    cur_time_str = datetime.now(pytz.timezone('America/Chicago')).strftime("%m/%d/%Y, %H:%M:%S")
    print("----Load cell phenotypes @ {}".format(cur_time_str))
    cell_phenotype_dict = None
    cell_phenotype_path = os.path.join(roi_root_dir, args.steinbock_dir, "cell_phenotypes.json")
    with open(cell_phenotype_path) as fp:
        cell_phenotype_dict = json.load(fp)
    cell_ids = [ele for ele in cell_phenotype_dict.keys()]
    cell_phenotypes = []
    library_ids = []
    cell_nums = []
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
    sq.gr.spatial_neighbors(rois_adata, spatial_key="spatial", library_key="library_id", coord_type="generic", radius=20.0)
    sq.gr.nhood_enrichment(rois_adata, cluster_key="cell_type", show_progress_bar=False)
    
    # save the plot
    plot_name = "cell_type_enrichment_{}_{}".format(args.path_stage, args.enrichment_mode)
    plot_path = os.path.join(enrich_result_dir, plot_name + args.plot_format)
    plot_title = "{} - Neighborhood Enrichment Analysis".format(args.path_stage)
    sq.pl.nhood_enrichment(rois_adata, cluster_key="cell_type", mode=args.enrichment_mode,
        title=plot_title, save=plot_path, dpi=300)
    cur_time_str = datetime.now(pytz.timezone('America/Chicago')).strftime("%m/%d/%Y, %H:%M:%S")
    print("----Finish @ {}".format(cur_time_str))    