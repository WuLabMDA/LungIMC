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
    stage_order_lst = ["Normal", "AAH", "AIS", "MIA", "ADC"]
    # print ROI number counts
    roi_stages = [ele for ele in roi_fea_df["ROI_Stage"].tolist()]
    for cur_stage in stage_order_lst:
        print("{} has {} ROIs".format(cur_stage, sum([ele==cur_stage for ele in roi_stages])))    

    normal_aah_pval_lst = []
    aah_ais_pval_lst = []
    ais_mia_pval_lst = []
    mia_adc_pval_lst = []
    
    # iterate each feature
    increase_features = []
    decrease_features = []
    for cur_fea_name in roi_fea_names:
        stage_dict = {}
        for cur_stage in roi_fea_df["ROI_Stage"].unique():
            stage_dict[cur_stage] = roi_fea_df[cur_fea_name][roi_fea_df["ROI_Stage"]==cur_stage].values
        normal_feas = stage_dict["Normal"]
        aah_feas = stage_dict["AAH"]
        ais_feas = stage_dict["AIS"]
        mia_feas = stage_dict["MIA"]
        adc_feas = stage_dict["ADC"]

        normal_aah_p = stats.ttest_ind(normal_feas, aah_feas)
        normal_aah_pval_lst.append(normal_aah_p.pvalue)
        normal_aah_d = 1 if np.mean(aah_feas) > np.mean(normal_feas) else -1

        aah_ais_p = stats.ttest_ind(aah_feas, ais_feas)
        aah_ais_pval_lst.append(aah_ais_p.pvalue)
        aah_ais_d = 1 if np.mean(ais_feas) > np.mean(aah_feas) else -1
        
        ais_mia_p = stats.ttest_ind(ais_feas, mia_feas)
        ais_mia_pval_lst.append(ais_mia_p.pvalue)
        ais_mia_d = 1 if np.mean(mia_feas) > np.mean(ais_feas) else -1

        mia_adc_p = stats.ttest_ind(mia_feas, adc_feas)
        mia_adc_pval_lst.append(mia_adc_p.pvalue)
        mia_adc_d = 1 if np.mean(adc_feas) > np.mean(mia_feas) else -1
        
        # find the maximum pval
        max_pval = max(normal_aah_p.pvalue, aah_ais_p.pvalue, ais_mia_p.pvalue, mia_adc_p.pvalue)
        sum_d = normal_aah_d + aah_ais_d + ais_mia_d + mia_adc_d
        if max_pval < args.pval_thresh:
            if sum_d == 4:
                increase_features.append(cur_fea_name)
            elif sum_d == -4:
                decrease_features.append(cur_fea_name)
            
    print("There are {} features significance increasing between featuers.".format(len(increase_features)))
    print(increase_features)
    print("There are {} features significance decreasing between featuers.".format(len(decrease_features)))
    print(decrease_features)

    monotonous_features = list(itertools.chain(decrease_features, increase_features))
    keep_features = ["ROI_ID", "ROI_Stage"]
    keep_features.extend(monotonous_features)
    keep_fea_df = roi_fea_df.loc[:, roi_fea_df.columns.isin(keep_features)]
    keep_fea_df["ROI_Stage"] = pd.Categorical(keep_fea_df["ROI_Stage"], stage_order_lst)
    keep_fea_df = keep_fea_df.sort_values("ROI_Stage")
    keep_fea_df[monotonous_features] -= keep_fea_df[monotonous_features].min()
    keep_fea_df[monotonous_features] /= keep_fea_df[monotonous_features].max()

    monotonous_fea_arr = (np.transpose(keep_fea_df[monotonous_features].to_numpy()) * 255.0).astype(np.uint8)
    cmap = plt.get_cmap('OrRd')
    monotonous_fea_arr = cmap(monotonous_fea_arr)
    monotonous_fea_map_path = os.path.join(feature_root_dir, "monotonous_fea_stage_roi_heatmap.pdf")
    fig, ax = plt.subplots(figsize=(12, 6))
    plt.imshow(monotonous_fea_arr)
    plt.savefig(monotonous_fea_map_path, transparent=False, dpi=600)
    # io.imsave(monotonous_fea_map_path, monotonous_fea_arr)

    stage_mean_df = keep_fea_df.groupby("ROI_Stage")[monotonous_features].mean()
    stage_mean_df = stage_mean_df.T
    fig, ax = plt.subplots(figsize=(3, 18))
    fig = plt.figure()
    sns.set(font_scale=0.2)
    sns.heatmap(stage_mean_df, cmap='RdYlGn_r', linewidths=0.3, annot=False)
    monotonous_mean_map_path = os.path.join(feature_root_dir, "monotonous_fea_stage_mean_heatmap.pdf")
    plt.savefig(monotonous_mean_map_path, transparent=False, dpi=600)