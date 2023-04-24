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
    stage_order_lst = ["Normal", "AAH", "AIS_MIA", "ADC"]

    normal_features = ['4-Eccentricity', 'Epithelial-Cell-Eccentricity', 'CT-SimpsonIndex', 'Fibroblast-Eccentricity', 'Endothelial-Cell-Proportion', 'Undefined-Eccentricity', 'CN4-Proportion', 'Epithelial-Cell-Proportion', 'CN2-Proportion', '9-VISTA']
    aah_features = ['3-TIM3', 'Proliferating-Cell-TIM3', 'MDSC-TIM3', '10-TIM3', 'Monocyte-TIM3', 'CD8-T-Cell-TIM3', '5-TIM3', 'CD4-T-Cell-TIM3', '4-TIM3', '1-TIM3']
    ais_features = ['4-HLADR', '2-HLADR', 'Endothelial-Cell-HLADR', '6-HLADR', 'Proliferating-Cell-HLADR', 'Epithelial-Cell-HLADR', 'MDSC-Proportion', '3-HLADR', 'CN6-Density', 'MDSC-Density']
    mia_features = ['CD8-T-Cell-Area', 'Undefined-MajorAxisLength', 'Undefined-MinorAxisLength', 'Endothelial-Cell-MinorAxisLength', 'Undefined-Area', 'Endothelial-Cell-Area', 'Endothelial-Cell-MajorAxisLength', 'CD4-T-Cell-Area', 'CN1-Density', '5-MajorAxisLength']
    adc_features = ['NK-Cell-Macrophage', 'NK-Cell-Fibroblast', 'NK-Cell-NK-Cell', 'Proliferating-Cell-NK-Cell', 'NK-Cell-CD8-T-Cell', 'NK-Cell-Proliferating-Cell', 'Proliferating-Cell-Proliferating-Cell', 'Proliferating-Cell-Density', 'NK-Cell-Endothelial-Cell', 'Fibroblast-NK-Cell']]
    dominant_features = normal_features + aah_features + ais_features + mia_features + adc_features

    keep_features = ["ROI_ID", "ROI_Stage"]
    keep_features.extend(dominant_features)
    keep_fea_df = roi_fea_df.loc[:, roi_fea_df.columns.isin(keep_features)]
    keep_fea_df["ROI_Stage"] = pd.Categorical(keep_fea_df["ROI_Stage"], stage_order_lst)
    keep_fea_df = keep_fea_df.sort_values("ROI_Stage")
    keep_fea_df[dominant_features] -= keep_fea_df[dominant_features].min()
    keep_fea_df[dominant_features] /= keep_fea_df[dominant_features].max()

    dominant_fea_arr = (np.transpose(keep_fea_df[dominant_features].to_numpy()) * 255.0).astype(np.uint8)
    cmap = plt.get_cmap('OrRd')
    dominant_fea_arr = cmap(dominant_fea_arr)
    dominant_fea_map_path = os.path.join(feature_root_dir, "dominant_fea_stage_roi_heatmap.pdf")
    fig, ax = plt.subplots(figsize=(12, 6))
    plt.imshow(dominant_fea_arr)
    plt.savefig(dominant_fea_map_path, transparent=False, dpi=600)