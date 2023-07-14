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
    parser.add_argument("--margin_core_dir",        type=str,       default="MarginCore")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)
    margin_core_dir = os.path.join(slide_agg_dir, args.margin_core_dir)
    if not os.path.exists(margin_core_dir):
        os.makedirs(margin_core_dir)

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
    # print("There are {} lesions.".format(len(lesion_names)))

    # collect all ROI's distance to core
    for cur_lesion in lesion_roi_dist_dict.keys():
        cur_lesion_dict = lesion_roi_dist_dict[cur_lesion]
        for cur_roi in cur_lesion_dict.keys():
            if cur_lesion_dict[cur_roi]["BorderDist"] > 1000:
                roi_distance_dict[cur_roi] = "TumorCore"
            else:
                roi_distance_dict[cur_roi] = "TumorMargin"
    dist_roi_names = [ele for ele in roi_distance_dict.keys()]
    # load lesion roi features
    lesion_roi_fea_path = os.path.join(slide_agg_dir, "lesion_roi_feas_raw.csv")
    roi_fea_df = pd.read_csv(lesion_roi_fea_path)
    roi_fea_df = roi_fea_df[roi_fea_df["ROI_ID"].isin(dist_roi_names)]
    roi_core_dist_lst = [roi_distance_dict[ele] for ele in roi_fea_df["ROI_ID"]]
    roi_fea_df.insert(loc = 3, column = "CoreMargin", value = roi_core_dist_lst)
    # interested feature list
    study_fea_lst = ["CD8-T-Cell-Proportion", "CD8-T-Cell-Density", "Macrophage-Proportion", "Macrophage-Density", 
                     "B-Cell-Proportion", "B-Cell-Density", "Endothelial-Cell-Proportion", "Endothelial-Cell-Density",
                     "CD8-T-Cell-Epithelial-Cell", "Macrophage-Epithelial-Cell", "Endothelial-Cell-Fibroblast", "Neutrophil-T-Reg-Cell"]
    # filter interested study features
    filter_col_names = [ele for ele in roi_fea_df.columns.tolist()[:4]] + study_fea_lst
    roi_fea_df = roi_fea_df[filter_col_names]
    # print("There are {} ROIs with features.".format(roi_fea_df.shape[0]))
    roi_fea_lst = [ele for ele in roi_fea_df["ROI_ID"].tolist()]

    # collect lesion features
    margin_core_cols = ["Lesion", "Loc",] + study_fea_lst
    margin_core_df = pd.DataFrame(columns = margin_core_cols)
    for lesion_name in lesion_roi_dist_dict.keys():
        lesion_info = lesion_roi_dist_dict[lesion_name]
        lesion_rois = [ele for ele in lesion_info.keys() if ele in roi_fea_lst]
        lesion_fea_df = roi_fea_df[roi_fea_df["ROI_ID"].isin(lesion_rois)]
        lesion_stage_lst = [ele for ele in lesion_fea_df["ROI_Stage"]]
        if lesion_stage_lst[0] != "ADC":
            continue
        lesion_core_df = lesion_fea_df[lesion_fea_df["CoreMargin"] == "TumorCore"]
        lesion_margin_df = lesion_fea_df[lesion_fea_df["CoreMargin"] == "TumorMargin"]
        if len(lesion_core_df) == 0 or len(lesion_margin_df) == 0:
            continue
        lesion_core_lst = [lesion_name, "TumorCore", ]
        lesion_margin_lst = [lesion_name, "TumorMargin", ]
        for cur_fea in study_fea_lst:
            lesion_core_lst.append(np.mean([ele for ele in lesion_core_df[cur_fea].tolist()]))
            lesion_margin_lst.append(np.mean([ele for ele in lesion_margin_df[cur_fea].tolist()]))
        margin_core_df.loc[len(margin_core_df)] = lesion_core_lst
        margin_core_df.loc[len(margin_core_df)] = lesion_margin_lst

    # save lesion information
    adc_lesion_loc_path = os.path.join(margin_core_dir, "adc_lesion_info.csv")
    margin_core_df.to_csv(adc_lesion_loc_path, index=False)