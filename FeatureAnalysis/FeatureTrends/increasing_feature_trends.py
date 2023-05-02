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
    # increase_fea_lst = ["T-Reg-Cell-Area", "Epithelial-Cell-Area", "T-Reg-Cell-MajorAxisLength"]
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    roi_stages = roi_fea_df["ROI_Stage"]
    stage_orders = ["Normal", "AAH", "AIS", "MIA", "ADC"]
    sns.boxplot(x=roi_stages, y=roi_fea_df["T-Reg-Cell-Area"], orient='v', order = stage_orders, ax=axes[0])
    sns.boxplot(x=roi_stages, y=roi_fea_df["Epithelial-Cell-Area"], orient='v', order = stage_orders, ax=axes[1])
    sns.boxplot(x=roi_stages, y=roi_fea_df["T-Reg-Cell-MajorAxisLength"], orient='v', order = stage_orders, ax=axes[2])
    increase_fea_box_path = os.path.join(feature_root_dir, "increasing_fea_boxplot.pdf")
    plt.savefig(increase_fea_box_path, transparent=False, dpi=300)