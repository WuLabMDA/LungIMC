# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, pickle
from datetime import datetime
import pandas as pd
import numpy as np
from scipy import stats
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import cosine


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])  
    parser.add_argument("--feature_dir",            type=str,       default="FeatureAnalysis")  
    parser.add_argument("--aggregation_dir",        type=str,       default="Aggregation")
    parser.add_argument("--ith_dir",                type=str,       default="ITH")
    parser.add_argument("--roi_num",                type=int,       default=10)
    parser.add_argument("--rand_seed",              type=int,       default=1234)

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()
    np.random.seed(args.rand_seed)

    # set directory
    dataset_dir = os.path.join(args.data_root, args.data_set)
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)
    lesion_ith_dir = os.path.join(slide_agg_dir, args.ith_dir)
    if not os.path.exists(lesion_ith_dir):
        os.makedirs(lesion_ith_dir)

    # load aggregated feature
    lesion_roi_fea_path = os.path.join(slide_agg_dir, "lesion_roi_feas.csv")
    roi_fea_df = pd.read_csv(lesion_roi_fea_path)
    # remove Normal ROIs
    roi_fea_df = roi_fea_df[roi_fea_df["ROI_Stage"] != "Normal"]
    roi_names = roi_fea_df["ROI_ID"].tolist()

    # locate all lesions    
    lesion_roi_dict = {}
    for cur_roi in roi_names:
        cur_leion = cur_roi[:-7]
        if cur_leion not in lesion_roi_dict:
            lesion_roi_dict[cur_leion] = [cur_roi, ]
        else:
            lesion_roi_dict[cur_leion].append(cur_roi)
    # remove lesion with multiple stages
    lesion_roi_dict.pop("2571-1D")

    # create empty slide dataframe
    slide_column_lst = ["LesionID", "SmokeStatus", "RawITH", "RandomITH"]
    slide_ith_df = pd.DataFrame(columns=slide_column_lst)

    for cur_lesion in lesion_roi_dict.keys():
        roi_inds = [roi_names.index(cur_roi) for cur_roi in lesion_roi_dict[cur_lesion]]
        cur_lesion_df = roi_fea_df.iloc[roi_inds]
        stage_lst = cur_lesion_df["ROI_Stage"].tolist()
        if stage_lst[0] != "ADC":
            continue
        smoke_lst = cur_lesion_df["SmokeStatus"].tolist()
        row_val_lst = [cur_lesion, smoke_lst[0]]
        lesion_fea_df = cur_lesion_df.iloc[:, 3:]
        lesion_fea_np = lesion_fea_df.to_numpy()
        # RawITH
        kernel_mat = 1.0 - pairwise_distances(lesion_fea_np, metric="cosine") / 2.0
        kernel_triu = kernel_mat[np.triu_indices(len(roi_inds), k = 1)]
        raw_ith = np.median(kernel_triu)
        row_val_lst.append(raw_ith)
        # RandomITH
        if lesion_fea_np.shape[0] <= args.roi_num:
            row_val_lst.append(raw_ith)
        else:
            np.random.shuffle(lesion_fea_np)
            random_fea_np = lesion_fea_np[:args.roi_num, :]
            random_kernel_mat = 1.0 - pairwise_distances(random_fea_np, metric="cosine") / 2.0
            random_kernel_triu = random_kernel_mat[np.triu_indices(args.roi_num, k = 1)]
            random_ith = np.median(random_kernel_triu)
            row_val_lst.append(random_ith)        
        slide_ith_df.loc[len(slide_ith_df.index)] = row_val_lst
    # Check slide dataframe information 
    lesion_ith_path = os.path.join(lesion_ith_dir, "lesion_ith_adc_{}.csv".format(args.roi_num))
    slide_ith_df.to_csv(lesion_ith_path, index = False)