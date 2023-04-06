# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz
from datetime import datetime
import pandas as pd
import numpy as np
from scipy import stats


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
    roi_fea_path = os.path.join(feature_root_dir, "ROI_Fea_Merge.csv")
    roi_fea_df = pd.read_csv(roi_fea_path)
    roi_fea_columns = [ele for ele in roi_fea_df.columns.tolist()]
    roi_fea_names = roi_fea_columns[2:] # exclude ROI_ID & ROI_Stage

    normal_aah_pval_lst = []
    aah_ais_pval_lst = []
    ais_mia_pval_lst = []
    mia_adc_pval_lst = []
    
    # iterate each feature
    common_features = []
    for cur_fea_name in roi_fea_names:
        stage_dict = {}
        for cur_stage in roi_fea_df["ROI_Stage"].unique():
            stage_dict[cur_stage] = roi_fea_df[cur_fea_name][roi_fea_df["ROI_Stage"]==cur_stage].values
        normal_feas = stage_dict["Normal"]
        aah_feas = stage_dict["AAH"]
        ais_feas = stage_dict["AIS"]
        mia_feas = stage_dict["MIA"]
        adc_feas = stage_dict["ADC"]
        import pdb; pdb.set_trace()
        # rest_normal_feas = 

        p_normal_aah = stats.ttest_ind(normal_feas, aah_feas)
        normal_aah_pval_lst.append(p_normal_aah.pvalue)
        p_aah_ais = stats.ttest_ind(aah_feas, ais_feas)
        aah_ais_pval_lst.append(p_aah_ais.pvalue)
        p_ais_mia = stats.ttest_ind(ais_feas, mia_feas)
        ais_mia_pval_lst.append(p_ais_mia.pvalue)
        p_mia_adc = stats.ttest_ind(mia_feas, adc_feas)
        mia_adc_pval_lst.append(p_mia_adc.pvalue)
        # find the maximum pval
        max_pval = max(p_normal_aah.pvalue, p_aah_ais.pvalue, p_ais_mia.pvalue, p_mia_adc.pvalue)
        if max_pval < 0.01:
            common_features.append(cur_fea_name)
    print("There are {} features with significance between featuers.".format(len(common_features)))

    
