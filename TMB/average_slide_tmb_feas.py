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
    parser.add_argument("--tmb_dir",                type=str,       default="TMB")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    # set directory
    dataset_dir = os.path.join(args.data_root, args.data_set)
    slide_tmb_dir = os.path.join(dataset_dir, args.tmb_dir)

    # load aggregated feature
    lesion_roi_fea_path = os.path.join(slide_tmb_dir, "lesion_roi_feas.csv")
    roi_fea_df = pd.read_csv(lesion_roi_fea_path)

    # create empty slide dataframe
    slide_column_lst = ["LesionID", "LesionStage", "TMB"]
    roi_fea_columns = [ele for ele in roi_fea_df.columns.tolist()]
    roi_fea_columns = roi_fea_columns[3:]
    slide_column_lst.extend(roi_fea_columns)
    slide_df = pd.DataFrame(columns=slide_column_lst)

    roi_lst = roi_fea_df["ROI_ID"].tolist()
    lesion_lst = list(set([ele[:-7] for ele in roi_lst]))    
    for cur_lesion in lesion_lst:
        roi_names = [roi_lst[ind] for ind in np.arange(len(roi_lst)) if roi_lst[ind].startswith(cur_lesion)]
        cur_lesion_df = roi_fea_df[roi_fea_df["ROI_ID"].isin(roi_names)]
        cur_tmbs = [ele for ele in cur_lesion_df["TMB"].tolist()]
        if len(set(cur_tmbs)) != 1:
            print("Multiple TMBs in {}".format(cur_lesion))
            continue
        stage_lst = cur_lesion_df["ROI_Stage"].tolist()
        distinct_stages = list(set(stage_lst))
        if len(distinct_stages) > 1:
            print("{} has mutliple stages.".format(cur_lesion))
            continue
        for cur_stage in distinct_stages:
            row_val_lst = [cur_lesion, cur_stage, cur_tmbs[0]]
            stage_rois = []
            for cur_roi, roi_stage in zip(roi_names, stage_lst):
                if roi_stage == cur_stage:
                    stage_rois.append(cur_roi)
            cur_lesion_df = roi_fea_df[roi_fea_df["ROI_ID"].isin(stage_rois)]
            lesion_fea_df = cur_lesion_df.iloc[:, 3:]
            for cur_fea in roi_fea_columns:
                row_val_lst.append(np.mean(cur_lesion_df[cur_fea].tolist()))
            slide_df.loc[len(slide_df.index)] = row_val_lst

    # Check slide dataframe information 
    lesion_fea_path = os.path.join(slide_tmb_dir, "lesion_avg_feas.csv")
    slide_df.to_csv(lesion_fea_path, index = False)