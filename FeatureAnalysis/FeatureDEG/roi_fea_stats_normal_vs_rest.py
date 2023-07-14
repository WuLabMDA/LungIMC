# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, itertools
from operator import itemgetter
from datetime import datetime
import pandas as pd
import numpy as np
from scipy import stats
from collections import OrderedDict
import matplotlib.pyplot as plt


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])  
    parser.add_argument("--feature_dir",            type=str,       default="FeatureAnalysis")
    parser.add_argument("--pval_thresh",            type=float,     default=0.05)     

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
    print("Of overall {} features:".format(len(roi_fea_names)))

    feature_lst = []
    # iterate each feature
    for cur_fea_name in roi_fea_names:
        stage_dict = {}
        for cur_stage in roi_fea_df["ROI_Stage"].unique():
            stage_dict[cur_stage] = roi_fea_df[cur_fea_name][roi_fea_df["ROI_Stage"]==cur_stage].values
        normal_feas = stage_dict["Normal"]
        aah_feas = stage_dict["AAH"]
        aismia_feas = stage_dict["AIS_MIA"]
        adc_feas = stage_dict["ADC"]
        # calculate mean values of each stage
        mean_normal = np.mean(normal_feas)
        mean_aah = np.mean(aah_feas)
        mean_aismia = np.mean(aismia_feas)
        mean_adc = np.mean(adc_feas)
        max_mean_rest = np.max([mean_aah, mean_aismia, mean_adc])
        if mean_normal > max_mean_rest:
            normal_aah_p = stats.ttest_ind(normal_feas, aah_feas)
            normal_aismia_p = stats.ttest_ind(normal_feas, aismia_feas)
            normal_adc_p = stats.ttest_ind(normal_feas, adc_feas)
            max_pval = np.max([normal_aah_p.pvalue, normal_aismia_p.pvalue, normal_adc_p.pvalue])
            if max_pval < args.pval_thresh:
                feature_lst.append(cur_fea_name)
    print("Normal has {} DEG significant features.".format(len(feature_lst)))
    print(feature_lst)
    
