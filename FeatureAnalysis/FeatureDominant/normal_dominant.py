# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, itertools
import math, operator
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
    parser.add_argument("--top_n",                  type=int,       default=10)       

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

    # iterate each feature
    feature_fc_dict = {}
    for cur_fea_name in roi_fea_names:
        stage_dict = {}
        for cur_stage in roi_fea_df["ROI_Stage"].unique():
            stage_dict[cur_stage] = roi_fea_df[cur_fea_name][roi_fea_df["ROI_Stage"]==cur_stage].values
        normal_feas = stage_dict["Normal"]
        aah_feas = stage_dict["AAH"]
        ais_feas = stage_dict["AIS"]
        mia_feas = stage_dict["MIA"]
        adc_feas = stage_dict["ADC"]

        # calculate mean values of each stage
        mean_normal = np.mean(normal_feas)
        mean_aah = np.mean(aah_feas)
        mean_ais = np.mean(ais_feas)
        mean_mia = np.mean(mia_feas)
        mean_adc = np.mean(adc_feas)
        max_mean_rest = np.max([mean_aah, mean_ais, mean_mia, mean_adc])
        if mean_normal > max_mean_rest:
            normal_aah_p = stats.ttest_ind(normal_feas, aah_feas)
            normal_ais_p = stats.ttest_ind(normal_feas, ais_feas)
            normal_mia_p = stats.ttest_ind(normal_feas, mia_feas)
            normal_adc_p = stats.ttest_ind(normal_feas, adc_feas)
            max_pval = np.max([normal_aah_p.pvalue, normal_ais_p.pvalue, normal_mia_p.pvalue, normal_adc_p.pvalue])
            if max_pval < args.pval_thresh:
                # stage_fc_lst = [mean_normal*1.0/mean_aah, mean_normal*1.0/mean_ais, mean_normal*1.0/mean_mia, mean_normal*1.0/mean_adc]
                # avg_log2fc = np.mean([math.log2(ele) for ele in stage_fc_lst])
                feature_fc_dict[cur_fea_name] = max_pval
    sorted_fea_dict = dict(sorted(feature_fc_dict.items(), key=operator.itemgetter(1),reverse=False))
    print("Normal has {} significant features.".format(len(sorted_fea_dict)))
    print("Top {} features are: ".format(args.top_n))
    feature_keys = [fea for fea in sorted_fea_dict.keys()]
    top_feas = [feature_keys[ind] for ind in np.arange(args.top_n)]
    print(top_feas)
    for feature_name in top_feas:
        fea_pval = sorted_fea_dict[feature_name]
        print("Feature {} with p-value: {}".format(feature_name, fea_pval))