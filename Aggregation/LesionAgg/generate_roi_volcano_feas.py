# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, pickle, math
from datetime import datetime
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats import multitest
from bioinfokit import visuz
import matplotlib.pyplot as plt
plt.rcParams['font.serif'] = ['Times New Roman']


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])                        
    parser.add_argument("--aggregation_dir",        type=str,       default="Aggregation")
    parser.add_argument("--rank_dir",               type=str,       default="Rank")
    parser.add_argument("--path_stage",             type=str,       default="AAH", choices=["AAH", "AIS", "MIA", "ADC"])

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)

    # ROI feature analysis
    lesion_roi_fea_path = os.path.join(slide_agg_dir, "lesion_roi_feas.csv")
    roi_fea_df = pd.read_csv(lesion_roi_fea_path)
    # normalize all features between 0.0 to 1.0
    roi_fea_columns = [ele for ele in roi_fea_df.columns.tolist()]
    roi_fea_names = roi_fea_columns[3:]
    for fea in roi_fea_names:
        roi_fea_df[fea] = (roi_fea_df[fea] + 3.0) / 6.0
    # filtering stage
    stage_fea_df = roi_fea_df[roi_fea_df["ROI_Stage"] == args.path_stage]
    never_inds, heavy_inds = [], []
    smoke_lst = [ele for ele in stage_fea_df["SmokeStatus"]]
    for ind in np.arange(len(smoke_lst)):
        if smoke_lst[ind] == "Never":
            never_inds.append(ind)
        else:
            heavy_inds.append(ind)
    # collect features 
    volcano_fea_lst = ["Feature", "Log2FC", "Pvalue"]
    fea_pval_lst, fea_fc_lst = [], []
    for cur_fea_name in roi_fea_names:
        cur_fea_lst = [ele for ele in stage_fea_df[cur_fea_name].tolist()]
        never_feas = [cur_fea_lst[ele] for ele in never_inds]
        heavy_feas = [cur_fea_lst[ele] for ele in heavy_inds]
        fea_log_fc = np.log2(np.mean(heavy_feas) / np.mean(never_feas))
        fea_ttest = stats.ttest_ind(heavy_feas, never_feas)
        fea_pval = fea_ttest.pvalue
        if math.isnan(fea_pval):
            fea_pval = 1.0
        fea_pval_lst.append(fea_pval)
        fea_fc_lst.append(fea_log_fc)
    # p-value FDR adjustment 
    pval_rejects, adjusted_pvals = multitest.fdrcorrection(fea_pval_lst)

    # collect signficant features
    stage_vol_df = pd.DataFrame(columns=volcano_fea_lst)
    for ind in np.arange(len(roi_fea_names)):
        cur_fea_name = roi_fea_names[ind]
        fea_log_fc = fea_fc_lst[ind]
        fea_pval = adjusted_pvals[ind]
        stage_vol_df.loc[len(stage_vol_df.index)] = [cur_fea_name, fea_log_fc, fea_pval]

    # Save ROI volcano features
    roi_rank_dir = os.path.join(slide_agg_dir, args.rank_dir, "ROI")
    if not os.path.exists(roi_rank_dir):
        os.makedirs(roi_rank_dir)
    roi_rank_volcano_path = os.path.join(roi_rank_dir, "{}_ROI_VolcanoFeatures.csv".format(args.path_stage))
    stage_vol_df.to_csv(roi_rank_volcano_path, index = False)