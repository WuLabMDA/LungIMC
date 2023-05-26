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
    parser.add_argument("--volcano_dir",            type=str,       default="Volcano")    
    parser.add_argument("--pval_thresh",            type=float,     default=0.05)         
    parser.add_argument("--path_stage",             type=str,       default="AAH")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)
    lesion_fea_path = os.path.join(slide_agg_dir, "lesion_bmean_feas.csv")
    lesion_fea_df = pd.read_csv(lesion_fea_path)

    # normalize all features between 0.0 to 1.0
    lesion_fea_columns = [ele for ele in lesion_fea_df.columns.tolist()]
    lesion_fea_names = lesion_fea_columns[3:]
    for fea in lesion_fea_names:
        lesion_fea_df[fea] = (lesion_fea_df[fea] + 3.0) / 6.0

    # filtering stage
    stage_fea_df = lesion_fea_df[lesion_fea_df["LesionStage"] == args.path_stage]
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
    for cur_fea_name in lesion_fea_names:
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

    heavy_ones = []
    never_ones = []
    stage_vol_df = pd.DataFrame(columns=volcano_fea_lst)
    for ind in np.arange(len(lesion_fea_names)):
        cur_fea_name = lesion_fea_names[ind]
        fea_log_fc = fea_fc_lst[ind]
        fea_pval = adjusted_pvals[ind]
        if fea_pval < args.pval_thresh:
            if fea_log_fc > 0.0:
                heavy_ones.append(cur_fea_name)
            else:
                never_ones.append(cur_fea_name)
        stage_vol_df.loc[len(stage_vol_df.index)] = [cur_fea_name, fea_log_fc, fea_pval]

    volcano_dir = os.path.join(slide_agg_dir, args.volcano_dir, "Bmean")
    if not os.path.exists(volcano_dir):
        os.makedirs(volcano_dir)

    # save feature information
    fea_txt_name = "{}_lesion_feas.txt".format(args.path_stage)
    fea_tx_path = os.path.join(volcano_dir, fea_txt_name)
    with open(fea_tx_path, 'w') as file:
        file.write("Slide on stage {}\n".format(args.path_stage))
        file.write("No. Heavy Smoker: {}\n".format(len(heavy_inds)))
        file.write("No. Never Smoker: {}\n".format(len(never_inds)))
        file.write("Minimum p-val before adjustment: {}\n".format(np.min(fea_pval_lst)))
        file.write("Minimum p-val after adjustment: {}\n".format(np.min(adjusted_pvals)))
        file.write("--Heavy smokers ({}) dominant features:\n {}\n".format(len(heavy_ones), heavy_ones))
        file.write("--Never smokers ({}) dominant features:\n {}\n".format(len(never_ones), never_ones))        
        
    # plot volcano 
    plot_name = "{}_lesion_volcano_plot".format(args.path_stage)
    fig_path = os.path.join(volcano_dir, plot_name)
    visuz.GeneExpression.volcano(df=stage_vol_df, lfc="Log2FC", pv="Pvalue", geneid="Feature", 
        lfc_thr=(0.0, 0.0), pv_thr=(0.05, 0.05), sign_line=True, 
        xlm=(-1.4, 1.5, 0.1), ylm=(0, 5, 1),
        gstyle=2, axtickfontsize=10,
        plotlegend=True, legendlabels=["Smoker Significant Up", "No Signficance", "Smoker Significant Down"],
        figname=fig_path, figtype="pdf")