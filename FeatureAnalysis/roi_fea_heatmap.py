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
    roi_fea_columns = [ele for ele in roi_fea_df.columns.tolist()]
    roi_fea_names = roi_fea_columns[2:] 

    # Ordering ROIs based on their stages
    stage_order_lst = ["Normal", "AAH", "AIS", "MIA", "ADC"]
    roi_fea_df["ROI_Stage"] = pd.Categorical(roi_fea_df["ROI_Stage"], stage_order_lst)
    roi_fea_df = roi_fea_df.sort_values("ROI_Stage")

    # min-max rescaling features
    roi_fea_df[roi_fea_names] -= roi_fea_df[roi_fea_names].min()
    roi_fea_df[roi_fea_names] /= roi_fea_df[roi_fea_names].max()

    # plot feature heatmap
    roi_fea_arr = (np.transpose(roi_fea_df[roi_fea_names].to_numpy()) * 255.0).astype(np.uint8)
    cmap = plt.get_cmap('OrRd')
    roi_fea_arr = cmap(roi_fea_arr)
    fea_heatmap_path = os.path.join(feature_root_dir, "fea_heatmap_stage.pdf")
    fig, ax = plt.subplots(figsize=(6, 8))
    plt.imshow(roi_fea_arr)
    plt.savefig(fea_heatmap_path, transparent=False, dpi=300)

    # stage_mean_df = keep_fea_df.groupby("ROI_Stage")[dominant_features].mean()
    # stage_mean_df = stage_mean_df.T
    # fig, ax = plt.subplots(figsize=(3, 18))
    # fig = plt.figure()
    # sns.set(font_scale=0.2)
    # sns.heatmap(stage_mean_df, cmap='RdYlGn_r', linewidths=0.3, annot=False)
    # dominant_mean_map_path = os.path.join(feature_root_dir, "dominant_fea_stage_mean_heatmap.pdf")
    # plt.savefig(dominant_mean_map_path, transparent=False, dpi=600)    