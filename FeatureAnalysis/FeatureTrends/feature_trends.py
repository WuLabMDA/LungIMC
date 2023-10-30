# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, itertools
from datetime import datetime
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from skimage import io
import seaborn as sns


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])  
    parser.add_argument("--feature_dir",            type=str,       default="FeatureAnalysis") 

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    feature_root_dir = os.path.join(dataset_dir, args.feature_dir)
    # load features
    roi_fea_path = os.path.join(feature_root_dir, "ROI_Fea_Aggregation.csv")
    roi_fea_df = pd.read_csv(roi_fea_path)
    # stage indices
    roi_stages = [ele for ele in roi_fea_df["ROI_Stage"].tolist()]
    normal_inds, aah_inds, ais_inds, mia_inds, adc_inds = [], [], [], [], []
    for ind, value in enumerate(roi_stages):
        if value == "Normal":
            normal_inds.append(ind)
        elif value == "AAH":
            aah_inds.append(ind)
        elif value == "AIS":
            ais_inds.append(ind)
        elif value == "MIA":
            mia_inds.append(ind)
        elif value == "ADC":
            adc_inds.append(ind)
        else:
            print("Unknown stage: {}".format(value))
    # plot directory
    feature_trend_dir = os.path.join(feature_root_dir, "FeatureTrends")
    if not os.path.exists(feature_trend_dir):
        os.makedirs(feature_trend_dir)
    # plot features
    feature_lst = ["Epithelial-Cell-Eccentricity",  "Endothelial-Cell-Proportion",  "Fibroblast-Eccentricity",  "Epithelial-Cell-Proportion",
                   "Macrophage-TIM3", "Monocyte-TIM3",  "CD4-T-Cell-TIM3",  "CD8-T-Cell-TIM3", "MDSC-TIM3",  "Dendritic-Cell-TIM3",  
                   "CD4-T-Cell-Area", "CD8-T-Cell-Area",  "Endothelial-Cell-Area",  "Endothelial-Cell-MajorAxisLength", "Fibroblast-NK-Cell",
                   "NK-Cell-Macrophage", "NK-Cell-NK-Cell", "NK-Cell-CD8-T-Cell", "NK-Cell-Fibroblast", "Proliferating-Cell-Density"]
    for cur_fea in feature_lst:
        # obtain kruskal wallis
        fea_vals = [ele for ele in roi_fea_df[cur_fea].tolist()]
        normal_feas = [fea_vals[ind] for ind in normal_inds]
        aah_feas = [fea_vals[ind] for ind in normal_inds]
        ais_feas = [fea_vals[ind] for ind in ais_inds]
        mia_feas = [fea_vals[ind] for ind in mia_inds]
        adc_feas = [fea_vals[ind] for ind in adc_inds]
        kw_test = stats.kruskal(normal_feas, aah_feas, ais_feas, mia_feas, adc_feas)
        # plot
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))
        roi_stages = roi_fea_df["ROI_Stage"]
        stage_orders = ["Normal", "AAH", "AIS", "MIA", "ADC"]
        sns.boxplot(x=roi_stages, y=roi_fea_df[cur_fea], orient='v', order = stage_orders, ax=ax)
        ax.set_title("{}-Pvalue:{}".format(cur_fea, kw_test.pvalue))
        fea_box_path = os.path.join(feature_trend_dir, "{}_fea_boxplot.pdf".format(cur_fea))
        plt.savefig(fea_box_path, transparent=False, dpi=300)