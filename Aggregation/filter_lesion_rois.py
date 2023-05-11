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

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    # load aggregated feature
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)
    agg_roi_fea_path = os.path.join(slide_agg_dir, "roi_fea_aggregation.csv")
    roi_fea_df = pd.read_csv(agg_roi_fea_path)

    # Add smoking information
    slide_smoke_path = os.path.join(slide_agg_dir, "lesion_smoke_info.xlsx")
    lesion_df = pd.read_excel(slide_smoke_path)
    lesion_slide_lst = [ele for ele in lesion_df["Slide_ID"].tolist()]
    smoke_slide_lst = [ele for ele in lesion_df["SmokeLevel"].tolist()]
    lesion_smoke_dict = {lesion: smoke for lesion, smoke in zip(lesion_slide_lst, smoke_slide_lst)}
    roi_smoke_lst = []
    lesion_roi = [ele for ele in roi_fea_df["ROI_ID"]]
    for cur_roi in lesion_roi:
        cur_lesion = cur_roi[:-7]
        cur_status = lesion_smoke_dict[cur_lesion]
        roi_smoke_lst.append(cur_status)
    roi_fea_df.insert(loc = 2, column = "SmokeStatus", value = roi_smoke_lst)
    roi_fea_df.loc[roi_fea_df["SmokeStatus"] == "Light", "SmokeStatus"] = "Never"
    roi_fea_df.loc[roi_fea_df["ROI_Stage"] == "DistantNormal", "ROI_Stage"] = "Normal"
    print("{} ROIs inside lesion.".format(len(roi_fea_df)))    
    lesion_roi_fea_path = os.path.join(slide_agg_dir, "lesion_roi_feas.xlsx")
    roi_fea_df.to_excel(lesion_roi_fea_path, index = False)