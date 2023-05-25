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


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])                        
    parser.add_argument("--recur_dir",              type=str,       default="Recurrence")
    parser.add_argument("--pval_thresh",            type=float,     default=0.05) 

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    slide_recur_dir = os.path.join(dataset_dir, args.recur_dir)
    lesion_roi_fea_path = os.path.join(slide_recur_dir, "lesion_roi_feas.csv")
    roi_fea_df = pd.read_csv(lesion_roi_fea_path)

    # normalize all features between 0.0 to 1.0
    roi_fea_columns = [ele for ele in roi_fea_df.columns.tolist()]
    roi_fea_names = roi_fea_columns[3:]
    for fea in roi_fea_names:
        roi_fea_df[fea] = (roi_fea_df[fea] + 3.0) / 6.0

    # filtering stage
    nonrecur_inds, recur_inds = [], []
    recur_lst = [ele for ele in roi_fea_df["Recur"]]
    for ind in np.arange(len(recur_lst)):
        if recur_lst[ind] == "Recurrence":
            recur_inds.append(ind)
        elif recur_lst[ind] == "NonRecurrence":
            nonrecur_inds.append(ind)
        else:
            print("Unknow status {}".format(recur_lst[ind]))

    # collect features 
    volcano_fea_lst = ["Feature", "Log2FC", "Pvalue"]
    fea_pval_lst = []
    fea_fc_lst = []
    for cur_fea_name in roi_fea_names:
        cur_fea_lst = [ele for ele in roi_fea_df[cur_fea_name].tolist()]
        nonrecur_feas = [cur_fea_lst[ele] for ele in nonrecur_inds]
        recur_feas = [cur_fea_lst[ele] for ele in recur_inds]
        fea_log_fc = np.log2(np.mean(recur_feas) / np.mean(nonrecur_feas))
        fea_ttest = stats.ttest_ind(recur_feas, nonrecur_feas)
        fea_pval = fea_ttest.pvalue
        if math.isnan(fea_pval):
            fea_pval = 1.0
        fea_pval_lst.append(fea_pval)
        fea_fc_lst.append(fea_log_fc)
        
    # p-value FDR adjustment 
    pval_rejects, adjusted_pvals = multitest.fdrcorrection(fea_pval_lst)

    # collect signficant features
    recur_ones = []
    nonrecur_ones = []
    stage_vol_df = pd.DataFrame(columns=volcano_fea_lst)
    for ind in np.arange(len(roi_fea_names)):
        cur_fea_name = roi_fea_names[ind]
        fea_log_fc = fea_fc_lst[ind]
        fea_pval = adjusted_pvals[ind]
        if fea_pval < args.pval_thresh:
            if fea_log_fc > 0.0:
                recur_ones.append(cur_fea_name)
            else:
                nonrecur_ones.append(cur_fea_name)
        stage_vol_df.loc[len(stage_vol_df.index)] = [cur_fea_name, fea_log_fc, fea_pval]

    volcano_dir = os.path.join(slide_recur_dir, "ROI-Volcano")
    if not os.path.exists(volcano_dir):
        os.makedirs(volcano_dir)
        
    # plot volcano 
    fig_path = os.path.join(volcano_dir, "recurrence_roi_volcano_plot")
    visuz.GeneExpression.volcano(df=stage_vol_df, lfc="Log2FC", pv="Pvalue", geneid="Feature", 
        lfc_thr=(0.0, 0.0), pv_thr=(args.pval_thresh, args.pval_thresh), sign_line=True, 
        xlm=(-0.7, 0.8, 0.1), ylm=(0, 10, 2),
        gstyle=2, axtickfontsize=10,
        plotlegend=True, legendlabels=["Smoker significant up", "No signficance", "Smoker significant down"],
        figname=fig_path, figtype="pdf")

    print("No. NonRecurrence ROIs: {}".format(len(nonrecur_inds)))
    print("No. Recurrence ROIs: {}".format(len(recur_inds)))
    print("Minimum p-val before adjustment: {}".format(np.min(fea_pval_lst)))    
    print("Minimum p-val after adjustment: {}".format(np.min(adjusted_pvals)))
    print("--Recurrence dominant features:")
    print(recur_ones)
    print("--NonRecurrence dominant features:")
    print(nonrecur_ones)