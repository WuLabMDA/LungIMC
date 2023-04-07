# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, itertools
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

    normal_pval_lst = []
    aah_pval_lst = []
    ais_pval_lst = []
    mia_pval_lst = []
    adc_pval_lst = []

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
        rest_normal_feas = list(itertools.chain(aah_feas, ais_feas, mia_feas, adc_feas))
        rest_aah_feas = list(itertools.chain(normal_feas, ais_feas, mia_feas, adc_feas))
        rest_ais_feas = list(itertools.chain(normal_feas, aah_feas, mia_feas, adc_feas))
        rest_mia_feas = list(itertools.chain(normal_feas, aah_feas, ais_feas, adc_feas))
        rest_adc_feas = list(itertools.chain(normal_feas, aah_feas, ais_feas, mia_feas))

        normal_test = stats.ttest_ind(normal_feas, rest_normal_feas)
        normal_pval_lst.append(normal_test.pvalue)
        aah_test = stats.ttest_ind(aah_feas, rest_aah_feas)
        aah_pval_lst.append(aah_test.pvalue)
        ais_test = stats.ttest_ind(ais_feas, rest_ais_feas)
        ais_pval_lst.append(ais_test.pvalue)
        mia_test = stats.ttest_ind(mia_feas, rest_mia_feas)
        mia_pval_lst.append(mia_test.pvalue)
        adc_test = stats.ttest_ind(adc_feas, rest_adc_feas)
        adc_pval_lst.append(adc_test.pvalue)

    # 
    print("Normal has {} features with significance.".format(sum([(ele < 0.01) for ele in normal_pval_lst])))
    print("AAH has {} features with significance.".format(sum([(ele < 0.01) for ele in aah_pval_lst])))
    print("AIS has {} features with significance.".format(sum([(ele < 0.01) for ele in ais_pval_lst])))
    print("MIA has {} features with significance.".format(sum([(ele < 0.01) for ele in mia_pval_lst])))
    print("ADC has {} features with significance.".format(sum([(ele < 0.01) for ele in adc_pval_lst])))                