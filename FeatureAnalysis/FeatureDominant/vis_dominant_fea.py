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
    roi_fea_names = roi_fea_columns[2:] # exclude ROI_ID & ROI_Stage
    stage_order_lst = ["Normal", "AAH", "AIS", "MIA", "ADC"]

    normal_features = ['Epithelial-Cell-Epithelial-Cell', '5-Eccentricity', 'Proliferating-Cell-Eccentricity', 'CN-Richness', 'Macrophage-Macrophage', 'Macrophage-Eccentricity', '1-Eccentricity', 'Undefined-Eccentricity', '7-Eccentricity', '5-MajorAxisLength']
    aah_features = ['Epithelial-Cell-Epithelial-Cell', 'Fibroblast-Fibroblast', 'Macrophage-Macrophage', 'CT-Richness', 'CN-ShannonEntropy', 'Proliferating-Cell-Eccentricity', '2-MinorAxisLength', '6-MinorAxisLength', '1-Eccentricity', '5-MajorAxisLength']
    ais_features = ['CT-Richness', 'Proliferating-Cell-MajorAxisLength', 'Proliferating-Cell-MinorAxisLength', '1-MinorAxisLength', 'Proliferating-Cell-Undefined', 'Endothelial-Cell-Neutrophil', '1-MajorAxisLength', 'Fibroblast-Neutrophil', 'Proliferating-Cell-Area', 'Epithelial-Cell-Undefined']
    mia_features = ['Undefined-Undefined', 'CT-Richness', 'Proliferating-Cell-Proliferating-Cell', '5-MinorAxisLength', '5-MajorAxisLength', 'Dendritic-Cell-Dendritic-Cell', '1-MinorAxisLength', '5-Area', 'Macrophage-Dendritic-Cell', '1-MajorAxisLength']
    adc_features = ['CN-Richness', 'CT-Richness', 'Proliferating-Cell-Proliferating-Cell', '2-MinorAxisLength', '6-MinorAxisLength', 'CD8-T-Cell-CD8-T-Cell', 'Undefined-Proliferating-Cell', 'Undefined-Fibroblast', 'Endothelial-Cell-Fibroblast', 'Undefined-Endothelial-Cell']
    dominant_features = normal_features + aah_features + ais_features + mia_features + adc_features

    keep_features = ["ROI_ID", "ROI_Stage"]
    keep_features.extend(dominant_features)
    keep_fea_df = roi_fea_df.loc[:, roi_fea_df.columns.isin(keep_features)]
    keep_fea_df["ROI_Stage"] = pd.Categorical(keep_fea_df["ROI_Stage"], stage_order_lst)
    keep_fea_df = keep_fea_df.sort_values("ROI_Stage")


    # keep_fea_df[dominant_features] -= keep_fea_df[dominant_features].min()
    # keep_fea_df[dominant_features] /= keep_fea_df[dominant_features].max()

    # dominant_fea_arr = (np.transpose(keep_fea_df[dominant_features].to_numpy()) * 255.0).astype(np.uint8)
    # cmap = plt.get_cmap('OrRd')
    # dominant_fea_arr = cmap(dominant_fea_arr)
    # dominant_fea_map_path = os.path.join(feature_root_dir, "dominant_fea_stage_roi_heatmap.pdf")
    # fig, ax = plt.subplots(figsize=(12, 6))
    # plt.imshow(dominant_fea_arr)
    # plt.savefig(dominant_fea_map_path, transparent=False, dpi=600)

    # stage_mean_df = keep_fea_df.groupby("ROI_Stage")[dominant_features].mean()
    # stage_mean_df = stage_mean_df.T
    # fig, ax = plt.subplots(figsize=(3, 18))
    # fig = plt.figure()
    # sns.set(font_scale=0.2)
    # sns.heatmap(stage_mean_df, cmap='RdYlGn_r', linewidths=0.3, annot=False)
    # dominant_mean_map_path = os.path.join(feature_root_dir, "dominant_fea_stage_mean_heatmap.pdf")
    # plt.savefig(dominant_mean_map_path, transparent=False, dpi=600)    