# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, pickle
import math, random
from datetime import datetime
import pandas as pd
import numpy as np
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import cosine


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])  
    parser.add_argument("--feature_dir",            type=str,       default="FeatureAnalysis")  
    parser.add_argument("--tmb_dir",                type=str,       default="TMB")
    parser.add_argument("--rand_roi_num",           type=int,       default=50)
    parser.add_argument("--rand_times",             type=int,       default=3)
    parser.add_argument("--rand_seed",              type=int,       default=1234)

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()
    random.seed(args.rand_seed)

    # set directory
    dataset_dir = os.path.join(args.data_root, args.data_set)
    slide_tmb_dir = os.path.join(dataset_dir, args.tmb_dir)

    # load aggregated feature
    lesion_roi_fea_path = os.path.join(slide_tmb_dir, "lesion_roi_feas.csv")
    roi_fea_df = pd.read_csv(lesion_roi_fea_path)
    roi_names = roi_fea_df["ROI_ID"].tolist()

    # locate all lesions    
    lesion_roi_dict = {}
    for cur_roi in roi_names:
        cur_leion = cur_roi[:-7]
        if cur_leion not in lesion_roi_dict:
            lesion_roi_dict[cur_leion] = [cur_roi, ]
        else:
            lesion_roi_dict[cur_leion].append(cur_roi)
    # print("There are {} different lesions.".format(len(lesion_roi_dict)))

    # create empty slide dataframe
    slide_column_lst = ["LesionID", "LesionStage", "TMB"]
    roi_fea_columns = [ele for ele in roi_fea_df.columns.tolist()]
    roi_fea_columns = roi_fea_columns[3:]
    slide_column_lst.extend(roi_fea_columns)
    slide_df = pd.DataFrame(columns=slide_column_lst)

    for cur_lesion in lesion_roi_dict.keys():
        lesion_roi_lst = lesion_roi_dict[cur_lesion]
        for rind in np.arange(args.rand_times):
            rand_rois = random.choices(lesion_roi_lst, k=args.rand_roi_num)
            rand_inds = []
            for ele in rand_rois:
                if ele in roi_names:
                    rand_inds.append(roi_names.index(ele))
            if len(rand_inds) < (args.rand_roi_num - 20):
                print("{} sampling {}".format(cur_lesion, len(rand_inds)))
            rand_lesion_df = roi_fea_df.iloc[rand_inds]
            tmb_lst = rand_lesion_df["TMB"].tolist()
            if len(set(tmb_lst)) != 1:
                print("Multiple TMB levels in {}".format(cur_lesion))
                continue
            stage_lst = rand_lesion_df["ROI_Stage"].tolist()
            if len(set(stage_lst)) != 1:
                print("Multiple stages in {}".format(cur_lesion))
                continue
            row_val_lst = [cur_lesion, stage_lst[0], tmb_lst[0]]
            # # average ROIs
            # for cur_fea in roi_fea_columns:
            #     row_val_lst.append(np.mean(rand_lesion_df[cur_fea].tolist()))
            lesion_fea_df = rand_lesion_df.iloc[:, 3:]
            lesion_fea_np = lesion_fea_df.to_numpy()
            kernel_mat = (1.0 - pairwise_distances(lesion_fea_np, metric="cosine") + 1.0) / 2.0
            roi_weights = 1.0 / np.sum(kernel_mat, axis=0)
            roi_weights = roi_weights /np.sum(roi_weights) # normalize
            weigh_fea = np.matmul(roi_weights, lesion_fea_np)
            row_val_lst.extend(weigh_fea.tolist())
            # add lesions     
            slide_df.loc[len(slide_df.index)] = row_val_lst
    # Check slide dataframe information 
    bootstrap_lesion_fea_path = os.path.join(slide_tmb_dir, "lesion_bootstrap_feas.csv")
    slide_df.to_csv(bootstrap_lesion_fea_path, index = False)    