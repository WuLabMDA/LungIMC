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
    print("-Normal has {} ROIs.".format(len(normal_inds)))
    print("-AAH has {} ROIs.".format(len(aah_inds)))
    print("-AIS has {} ROIs.".format(len(ais_inds)))
    print("-MIA has {} ROIs.".format(len(mia_inds)))
    print("-ADC has {} ROIs.".format(len(adc_inds)))

    all_fea_df = roi_fea_df.iloc[:, 2:]
    fea_names = [ele for ele in all_fea_df.columns.tolist()]
    print("There are {} features in total.".format(len(fea_names)))    

    # # Kruskal Wallis Test
    # kw_fea_inds = []
    # for ind, fea_name in enumerate(fea_names):
    #     fea_vals = [ele for ele in all_fea_df[fea_name].tolist()]
    #     normal_feas = [fea_vals[ind] for ind in normal_inds]
    #     aah_feas = [fea_vals[ind] for ind in normal_inds]
    #     ais_feas = [fea_vals[ind] for ind in ais_inds]
    #     mia_feas = [fea_vals[ind] for ind in mia_inds]
    #     adc_feas = [fea_vals[ind] for ind in adc_inds]
    #     kw_test = stats.kruskal(normal_feas, aah_feas, ais_feas, mia_feas, adc_feas)
    #     if kw_test.pvalue < 0.05:
    #         kw_fea_inds.append(ind)
    # print("{}/{} pass Kruskal Wallis test.".format(len(kw_fea_inds), len(fea_names)))
    # all_fea_df = all_fea_df.iloc[:, kw_fea_inds]
    # print("KW feature shape:")
    # print(all_fea_df.shape)

    all_fea_df.index = roi_stages
    all_fea_arr = all_fea_df.to_numpy().astype("float").T
    import pdb; pdb.set_trace()
    plt.imshow(all_fea_arr, cmap ="bwr", interpolation='none')
    plt.colorbar()
    fea_heatmap_path = os.path.join(feature_root_dir, "fea_heatmap_stage-new.pdf")
    plt.savefig(fea_heatmap_path, transparent=False, dpi=300)