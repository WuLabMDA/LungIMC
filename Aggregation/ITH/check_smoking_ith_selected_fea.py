# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, pickle, math
from datetime import datetime
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats import multitest
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import cosine


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])  
    parser.add_argument("--feature_dir",            type=str,       default="FeatureAnalysis")  
    parser.add_argument("--aggregation_dir",        type=str,       default="Aggregation")
    parser.add_argument("--ith_dir",                type=str,       default="ITH")
    parser.add_argument("--pval_thresh",            type=float,     default=0.05)
    parser.add_argument("--dominant",               type=str,       default="Heavy", choices=["Heavy", "Never"])
    parser.add_argument("--path_stage",             type=str,       default="AAH", choices=["AAH", "AIS", "MIA", "ADC"])
    
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    # set directory
    dataset_dir = os.path.join(args.data_root, args.data_set)
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)
    lesion_ith_dir = os.path.join(slide_agg_dir, args.ith_dir)
    lesion_ith_path = os.path.join(lesion_ith_dir, "lesion_ith_per_fea.csv")
    slide_ith_df = pd.read_csv(lesion_ith_path)

    # obtain feature names
    lesion_fea_columns = [ele for ele in slide_ith_df.columns.tolist()]
    lesion_fea_names = lesion_fea_columns[3:]

    # filtering stage
    print("--Stage {}".format(args.path_stage))
    stage_fea_df = slide_ith_df[slide_ith_df["LesionStage"] == args.path_stage]
    never_inds, heavy_inds = [], []
    smoke_lst = [ele for ele in stage_fea_df["SmokeStatus"]]
    for ind in np.arange(len(smoke_lst)):
        if smoke_lst[ind] == "Never":
            never_inds.append(ind)
        else:
            heavy_inds.append(ind)
    print("--{} never smokers & {} heavy smokers".format(len(never_inds), len(heavy_inds)))
    # select features with significance
    p_fea_names = []
    for cur_fea_name in lesion_fea_names:
        cur_fea_lst = [ele for ele in stage_fea_df[cur_fea_name].tolist()]
        never_feas = [cur_fea_lst[ele] for ele in never_inds]
        heavy_feas = [cur_fea_lst[ele] for ele in heavy_inds]
        if args.dominant == "Heavy":
            if np.mean(never_feas) > np.mean(heavy_feas):
                continue
        elif args.dominant == "Never":
            if np.mean(heavy_feas) > np.mean(never_feas):
                continue            
        else:
            assert("unknow category: {}".format(args.dominant))
        fea_ttest = stats.ttest_ind(heavy_feas, never_feas)
        fea_pval = fea_ttest.pvalue
        if math.isnan(fea_pval):
            fea_pval = 1.0
        if fea_pval < args.pval_thresh:
            p_fea_names.append(cur_fea_name)
    if len(p_fea_names) == 0:
        print("No features selected.")
        sys.exit(1)
    else:
        print("--{} features selected for ITH analysis:".format(len(p_fea_names)))
        print(p_fea_names)

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
    slide_column_lst = ["LesionID", "LesionStage", "SmokeStatus", "LesionITH"]
    slide_ith_df = pd.DataFrame(columns=slide_column_lst)

    for cur_lesion in lesion_roi_dict.keys():
        roi_inds = [roi_names.index(cur_roi) for cur_roi in lesion_roi_dict[cur_lesion]]
        cur_lesion_df = roi_fea_df.iloc[roi_inds]
        stage_lst = cur_lesion_df["ROI_Stage"].tolist()
        smoke_lst = cur_lesion_df["SmokeStatus"].tolist()
        row_val_lst = [cur_lesion, stage_lst[0], smoke_lst[0]]
        # lesion_fea_df = cur_lesion_df.iloc[:, 3:] 
        lesion_select_df = cur_lesion_df[p_fea_names]
        lesion_fea_np = lesion_select_df.to_numpy()
        kernel_mat = 1.0 - pairwise_distances(lesion_fea_np, metric="cosine") / 2.0
        kernel_triu = kernel_mat[np.triu_indices(len(roi_inds), k = 1)]
        row_val_lst.append(np.median(kernel_triu))
        slide_ith_df.loc[len(slide_ith_df.index)] = row_val_lst
    # # Check slide dataframe information 
    # lesion_ith_path = os.path.join(lesion_ith_dir, "lesion_ith_select_fea.csv")
    # slide_ith_df.to_csv(lesion_ith_path, index = False)            

    # filtering stage
    # print("--Stage {}".format(args.path_stage))
    stage_fea_df = slide_ith_df[slide_ith_df["LesionStage"] == args.path_stage]
    never_inds, heavy_inds = [], []
    smoke_lst = [ele for ele in stage_fea_df["SmokeStatus"]]
    for ind in np.arange(len(smoke_lst)):
        if smoke_lst[ind] == "Never":
            never_inds.append(ind)
        else:
            heavy_inds.append(ind)
    # print("{} never smokers".format(len(never_inds)))
    # print("{} heavy smokers".format(len(heavy_inds)))

    # measure signficance
    cur_fea_lst = [ele for ele in stage_fea_df["LesionITH"].tolist()]
    never_feas = [cur_fea_lst[ele] for ele in never_inds]
    heavy_feas = [cur_fea_lst[ele] for ele in heavy_inds]
    fea_ttest = stats.ttest_ind(heavy_feas, never_feas)
    fea_pval = fea_ttest.pvalue
    if math.isnan(fea_pval):
        fea_pval = 1.0  
    print("--{} Dominant ITH signficance: {}".format(args.dominant, fea_pval))