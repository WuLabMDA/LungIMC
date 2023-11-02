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

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)
    
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
    for cur_lesion in lesion_roi_dist_dict.keys():
        cur_lesion_dict = lesion_roi_dist_dict[cur_lesion]
        for cur_roi in cur_lesion_dict.keys():
            roi_core_dist = cur_lesion_dict[cur_roi]["CoreDist"]
            roi_border_dist = cur_lesion_dict[cur_roi]["BorderDist"]
            core_dist_ratio = roi_core_dist * 1.0 / (roi_core_dist + roi_border_dist)
            roi_distance_dict[cur_roi] = core_dist_ratio
    dist_roi_names = [ele for ele in roi_distance_dict.keys()]
    # print("There are {} rois.".format(len(dist_roi_names)))

    # load lesion roi features
    lesion_roi_fea_path = os.path.join(slide_agg_dir, "lesion_roi_feas.csv")
    roi_fea_df = pd.read_csv(lesion_roi_fea_path)
    # fea_roi_lst = [ele for ele in roi_fea_df["ROI_ID"]]
    # for ind, cur_roi in enumerate(fea_roi_lst):
    #     if cur_roi not in dist_roi_names and roi_fea_df.loc[ind, "ROI_Stage"] != "DistantNormal":
    #         print(cur_roi)
    roi_fea_df = roi_fea_df[roi_fea_df["ROI_ID"].isin(dist_roi_names)]
    roi_core_dist_lst = [roi_distance_dict[ele] for ele in roi_fea_df["ROI_ID"]]
    roi_fea_df.insert(loc = 3, column = "CoreDist", value = roi_core_dist_lst)
    # save features with distance information
    lesion_roi_dist_fea_path = os.path.join(slide_agg_dir, "lesion_roi_dist_feas.csv")
    roi_fea_df.to_csv(lesion_roi_dist_fea_path, index = False)    