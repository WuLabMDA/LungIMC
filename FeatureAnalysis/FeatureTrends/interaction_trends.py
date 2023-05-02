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
    interaction_fea_lst = ["ROI_Stage", "Epithelial-Cell-T-Reg-Cell", "Epithelial-Cell-Macrophage", "Epithelial-Cell-B-Cell"]
    interaction_df = roi_fea_df[interaction_fea_lst]

    sns.lineplot(x="ROI_Stage", y="value", hue="variable", data=pd.melt(interaction_df, ["ROI_Stage"]))
    interaction_trends_path = os.path.join(feature_root_dir, "interaction_trends.pdf")
    plt.savefig(interaction_trends_path, transparent=False, dpi=300)