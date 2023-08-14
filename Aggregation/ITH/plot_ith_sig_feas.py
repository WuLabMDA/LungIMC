# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, pickle, math
from datetime import datetime
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats import multitest
import matplotlib.pyplot as plt


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])  
    parser.add_argument("--aggregation_dir",        type=str,       default="Aggregation")
    parser.add_argument("--ith_dir",                type=str,       default="ITH")
    parser.add_argument("--pval_thresh",            type=float,     default=0.05)

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
    
    # iterate each feature
    sig_ith_fea_lst = []
    lesion_stage_lst = [ele for ele in slide_ith_df["LesionStage"].tolist()]
    for cur_fea in lesion_fea_names:
        aah_lst, ais_lst, mia_lst, adc_lst = [], [], [], []
        fea_ith_lst = [ele for ele in slide_ith_df[cur_fea].tolist()]
        for cur_stage, cur_fea_val in zip(lesion_stage_lst, fea_ith_lst):
            if cur_stage == "AAH":
                aah_lst.append(cur_fea_val)
            elif cur_stage == "AIS":
                ais_lst.append(cur_fea_val)
            elif cur_stage == "MIA":
                mia_lst.append(cur_fea_val)
            elif cur_stage == "ADC":
                adc_lst.append(cur_fea_val)
            else:
                sys.exit("Unknow stage: {}".format(cur_stage))
        # Kruskal Wallis Test
        kw_test = stats.kruskal(aah_lst, ais_lst, mia_lst, adc_lst)
        if kw_test.pvalue < 0.05:
            sig_ith_fea_lst.append(cur_fea)
            # plot
            fig, ax = plt.subplots(figsize =(10, 7))
            bp = ax.boxplot([aah_lst, ais_lst, mia_lst, adc_lst])
            ax.set_xticklabels(["AAH", "AIS", "MIA", "ADC"])
            plt.title("{} ITH KW Test P-val: {}".format(cur_fea, kw_test.pvalue))

            # prepare folder
            fea_ith_plot_dir = os.path.join(lesion_ith_dir, "FeaITH")
            if "CN" in cur_fea:
                fea_ith_plot_dir = os.path.join(fea_ith_plot_dir, "CN")
            else:
                fea_ith_plot_dir = os.path.join(fea_ith_plot_dir, "CT")
            if "Proportion" in cur_fea:
                fea_ith_plot_dir = os.path.join(fea_ith_plot_dir, "Proportion")
            elif "Density" in cur_fea:
                fea_ith_plot_dir = os.path.join(fea_ith_plot_dir, "Density")
            else:
                fea_ith_plot_dir = os.path.join(fea_ith_plot_dir, "Others")
            if not os.path.exists(fea_ith_plot_dir):
                os.makedirs(fea_ith_plot_dir)
            fea_ith_path = os.path.join(fea_ith_plot_dir, "{}.pdf".format(cur_fea))
            plt.savefig(fea_ith_path, transparent=False, dpi=300)
