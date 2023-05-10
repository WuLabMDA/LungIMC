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
import matplotlib.font_manager


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])                        
    parser.add_argument("--aggregation_dir",        type=str,       default="Aggregation")
    parser.add_argument("--pval_thresh",            type=float,     default=0.05)       
    parser.add_argument("--path_stage",             type=str,       default="AAH")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)
    lesion_roi_fea_path = os.path.join(slide_agg_dir, "lesion_roi_feas.xlsx")
    roi_fea_df = pd.read_excel(lesion_roi_fea_path)

    # normalize all features between 0.0 to 1.0
    roi_fea_columns = [ele for ele in roi_fea_df.columns.tolist()]
    roi_fea_names = roi_fea_columns[3:]
    for fea in roi_fea_names:
        roi_fea_df[fea] = (roi_fea_df[fea] - roi_fea_df[fea].mean()) / roi_fea_df[fea].std(ddof=0)
        roi_fea_df[fea] = (roi_fea_df[fea].clip(-2.0, 2.0) + 2.0) / 4.0

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
    fea_pval_lst = []
    fea_fc_lst = []
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
    print("Minimum p-val before adjustment: {}".format(np.min(fea_pval_lst)))    
    pval_rejects, adjusted_pvals = multitest.fdrcorrection(fea_pval_lst)
    print("Minimum p-val after adjustment: {}".format(np.min(adjusted_pvals)))

    heavy_ones = []
    never_ones = []
    stage_vol_df = pd.DataFrame(columns=volcano_fea_lst)
    for ind in np.arange(len(roi_fea_names)):
        cur_fea_name = roi_fea_names[ind]
        fea_log_fc = fea_fc_lst[ind]
        fea_pval = adjusted_pvals[ind]
        if fea_pval < args.pval_thresh:
            if fea_log_fc > 0.0:
                heavy_ones.append(cur_fea_name)
            else:
                never_ones.append(cur_fea_name)
        stage_vol_df.loc[len(stage_vol_df.index)] = [cur_fea_name, fea_log_fc, fea_pval]

    volcano_dir = os.path.join(slide_agg_dir, "ROI-Volcano")
    if not os.path.exists(volcano_dir):
        os.makedirs(volcano_dir)
    # plot volcano 
    plot_name = "{}_roi_volcano_plot".format(args.path_stage)
    fig_path = os.path.join(volcano_dir, plot_name)
    visuz.GeneExpression.volcano(df=stage_vol_df, lfc="Log2FC", pv="Pvalue", geneid="Feature", 
        lfc_thr=(0.0, 0.0), pv_thr=(args.pval_thresh, args.pval_thresh), sign_line=True, 
        xlm=(-0.7, 0.8, 0.1), ylm=(0, 10, 2),
        gstyle=2, axtickfontsize=10,
        plotlegend=True, legendlabels=["Smoker significant up", "No signficance", "Smoker significant down"],
        figname=fig_path, figtype="pdf")
    
    print("ROI on stage {}".format(args.path_stage))
    print("--Heavy smokers dominant features:")
    print(heavy_ones)
    print("Never smokers dominant features:")
    print(never_ones)