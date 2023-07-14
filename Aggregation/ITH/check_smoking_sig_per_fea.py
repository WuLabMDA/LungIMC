# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, pickle, math
from datetime import datetime
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats import multitest


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])  
    parser.add_argument("--feature_dir",            type=str,       default="FeatureAnalysis")  
    parser.add_argument("--aggregation_dir",        type=str,       default="Aggregation")
    parser.add_argument("--ith_dir",                type=str,       default="ITH")
    parser.add_argument("--pval_thresh",            type=float,     default=0.05)
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
    print("{} never smokers".format(len(never_inds)))
    print("{} heavy smokers".format(len(heavy_inds)))

    # collect features 
    fea_pval_lst = []
    p_fea_names = []
    for cur_fea_name in lesion_fea_names:
        cur_fea_lst = [ele for ele in stage_fea_df[cur_fea_name].tolist()]
        never_feas = [cur_fea_lst[ele] for ele in never_inds]
        heavy_feas = [cur_fea_lst[ele] for ele in heavy_inds]
        fea_ttest = stats.ttest_ind(heavy_feas, never_feas)
        fea_pval = fea_ttest.pvalue
        if math.isnan(fea_pval):
            fea_pval = 1.0  
        fea_pval_lst.append(fea_pval)
        if fea_pval < args.pval_thresh:
            p_fea_names.append(cur_fea_name)


    # p-value FDR adjustment 
    pval_rejects, adjusted_pvals = multitest.fdrcorrection(fea_pval_lst)    
    sig_fea_num = np.sum([ele < args.pval_thresh for ele in adjusted_pvals])
    print("{} features with signficance ITH".format(sig_fea_num))
    print("Before FDR, {} significant features:".format(len(p_fea_names)))
    print(p_fea_names)