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

    normal_pval_dict = {}
    aah_pval_dict = {}
    aismia_pval_dict = {}
    adc_pval_dict = {}

    # iterate each feature
    for cur_fea_name in roi_fea_names:
        stage_dict = {}
        for cur_stage in roi_fea_df["ROI_Stage"].unique():
            stage_dict[cur_stage] = roi_fea_df[cur_fea_name][roi_fea_df["ROI_Stage"]==cur_stage].values
        normal_feas = stage_dict["Normal"]
        aah_feas = stage_dict["AAH"]
        aismia_feas = stage_dict["AIS_MIA"]
        adc_feas = stage_dict["ADC"]
        rest_normal_feas = list(itertools.chain(aah_feas, aismia_feas, adc_feas))
        rest_aah_feas = list(itertools.chain(normal_feas, aismia_feas, adc_feas))
        rest_aismia_feas = list(itertools.chain(normal_feas, aah_feas, adc_feas))
        rest_adc_feas = list(itertools.chain(normal_feas, aah_feas, aismia_feas))

        if np.mean(normal_feas) > np.mean(rest_normal_feas):
            normal_test = stats.ttest_ind(normal_feas, rest_normal_feas)
            normal_pval_dict[cur_fea_name] = normal_test.pvalue    
        if np.mean(aah_feas) > np.mean(rest_aah_feas):        
            aah_test = stats.ttest_ind(aah_feas, rest_aah_feas)
            aah_pval_dict[cur_fea_name] = aah_test.pvalue
        if np.mean(aismia_feas) > np.mean(rest_aismia_feas):
            aismia_test = stats.ttest_ind(aismia_feas, rest_aismia_feas)
            aismia_pval_dict[cur_fea_name] = aismia_test.pvalue
        if np.mean(adc_feas) > np.mean(rest_adc_feas):
            adc_test = stats.ttest_ind(adc_feas, rest_adc_feas)
            adc_pval_dict[cur_fea_name] = adc_test.pvalue
    # merge
    top_num = 10
    top_normal = dict(sorted(normal_pval_dict.items(), key=itemgetter(1))[:top_num])
    top_normal_feas = [ele for ele in top_normal.keys()]
    top_aah = dict(sorted(aah_pval_dict.items(), key=itemgetter(1))[:top_num])
    top_aah_feas = [ele for ele in top_aah.keys()]
    top_aismia = dict(sorted(aismia_pval_dict.items(), key=itemgetter(1))[:top_num])
    top_aismia_feas = [ele for ele in top_aismia.keys()]    
    top_adc = dict(sorted(adc_pval_dict.items(), key=itemgetter(1))[:top_num])
    top_adc_feas = [ele for ele in top_adc.keys()]
    merge_fea_lst = list(itertools.chain(top_normal_feas, top_aah_feas, top_aismia_feas, top_adc_feas))
    merge_fea_lst = list(dict.fromkeys(merge_fea_lst))

    keep_features = ["ROI_ID", "ROI_Stage"]
    keep_features.extend(merge_fea_lst)     
    keep_fea_df = roi_fea_df.loc[:, roi_fea_df.columns.isin(keep_features)]
    stage_order_lst = ["Normal", "AAH", "AIS_MIA", "ADC"]
    keep_fea_df["ROI_Stage"] = pd.Categorical(keep_fea_df["ROI_Stage"], stage_order_lst)
    keep_fea_df = keep_fea_df.sort_values("ROI_Stage")
    keep_fea_df[merge_fea_lst] -= keep_fea_df[merge_fea_lst].min()
    keep_fea_df[merge_fea_lst] /= keep_fea_df[merge_fea_lst].max()

    keep_fea_arr = (np.transpose(keep_fea_df[merge_fea_lst].to_numpy()) * 255.0).astype(np.uint8)
    cmap = plt.get_cmap('OrRd')
    keep_fea_arr = cmap(keep_fea_arr)
    keep_fea_map_path = os.path.join(feature_root_dir, "one_rest_fea_heatmap.pdf")
    fig, ax = plt.subplots(figsize=(12, 6))
    plt.imshow(keep_fea_arr)
    plt.savefig(keep_fea_map_path, transparent=False, dpi=600)    