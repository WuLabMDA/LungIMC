# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, pickle
from datetime import datetime
import pandas as pd
import numpy as np


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])  
    parser.add_argument("--feature_dir",            type=str,       default="FeatureAnalysis")  
    parser.add_argument("--aggregation_dir",        type=str,       default="Aggregation")
    parser.add_argument("--dist_dir",               type=str,       default="Distance")
    parser.add_argument("--dist_prog_dir",          type=str,       default="DistProgression")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)
    slide_distprog_dir = os.path.join(slide_agg_dir, args.dist_prog_dir)
    if not os.path.exists(slide_distprog_dir):
        os.makedirs(slide_distprog_dir)

    # collect roi percentage to core
    distance_dir = os.path.join(slide_agg_dir, args.dist_dir)
    lesion_roi_dist_path = os.path.join(distance_dir, "lesion_roi_dist.json")
    lesion_roi_dist_dict = None
    with open(lesion_roi_dist_path) as handle:
        lesion_roi_dist_dict = json.load(handle)
    roi_distance_dict = {}
    lesion_names = lesion_roi_dist_dict.keys()
    lesion_roi_dist_dict.pop("2571-1D")
    lesion_names = sorted(lesion_names)
    print("There are {} lesions.".format(len(lesion_names)))


    # collect all ROI's distance to core
    for cur_lesion in lesion_roi_dist_dict.keys():
        cur_lesion_dict = lesion_roi_dist_dict[cur_lesion]
        for cur_roi in cur_lesion_dict.keys():
            roi_distance_dict[cur_roi] = cur_lesion_dict[cur_roi]["CoreDist"]
    dist_roi_names = [ele for ele in roi_distance_dict.keys()]
    # load lesion roi features
    lesion_roi_fea_path = os.path.join(slide_agg_dir, "lesion_roi_feas.csv")
    roi_fea_df = pd.read_csv(lesion_roi_fea_path)
    roi_fea_df = roi_fea_df[roi_fea_df["ROI_ID"].isin(dist_roi_names)]
    roi_core_dist_lst = [roi_distance_dict[ele] for ele in roi_fea_df["ROI_ID"]]
    roi_fea_df.insert(loc = 3, column = "CoreDist", value = roi_core_dist_lst)
    # interested feature list
    study_fea_lst = ["CD8-T-Cell-Proportion", "CD8-T-Cell-Density", "Macrophage-Proportion", "Macrophage-Density", 
                     "B-Cell-Proportion", "B-Cell-Density", "Endothelial-Cell-Proportion", "Endothelial-Cell-Density",
                     "CD8-T-Cell-Epithelial-Cell", "Macrophage-Epithelial-Cell", "Endothelial-Cell-Fibroblast", "Neutrophil-T-Reg-Cell"]
    # filter interested study features
    filter_col_names = [ele for ele in roi_fea_df.columns.tolist()[:4]] + study_fea_lst
    roi_fea_df = roi_fea_df[filter_col_names]
    print("There are {} ROIs with features.".format(roi_fea_df.shape[0]))
    roi_fea_lst = [ele for ele in roi_fea_df["ROI_ID"].tolist()]
    # save features with distance information
    lesion_roi_core_dist_fea_path = os.path.join(slide_distprog_dir, "lesion_roi_coredist_feas.csv")
    roi_fea_df.to_csv(lesion_roi_core_dist_fea_path, index = False) 

    core_distprog_dict = {}
    ttl_roi_num = 0
    for lesion_name in lesion_roi_dist_dict.keys():
        lesion_info = lesion_roi_dist_dict[lesion_name]
        lesion_rois = [ele for ele in lesion_info.keys() if ele in roi_fea_lst]
        ttl_roi_num += len(lesion_rois)
        lesion_fea_df = roi_fea_df[roi_fea_df["ROI_ID"].isin(lesion_rois)]
        lesion_stage_lst = [ele for ele in lesion_fea_df["ROI_Stage"]]
        if len(set(lesion_stage_lst)) != 1:
            print(" has different roi stages.".format(lesion_name))
        lesion_smoke_lst = [ele for ele in lesion_fea_df["SmokeStatus"]]
        if len(set(lesion_smoke_lst)) != 1:
            print(" has different roi smoke status.".format(lesion_name))
        lesion_dict = {}
        lesion_dict["LesionROIs"] = sorted(lesion_rois)
        lesion_dict["LesionStage"] = lesion_stage_lst[0]
        lesion_dict["SmokeStatus"] = lesion_smoke_lst[0]
        core_distprog_dict[lesion_name] = lesion_dict
    print("There are {} ROIs with distances to core.".format(ttl_roi_num))
    
    # save lesion information
    lesion_core_distprog_path = os.path.join(slide_distprog_dir, "lesion_roi_info.json")
    with open(lesion_core_distprog_path, "w") as fp:
        json.dump(core_distprog_dict, fp)