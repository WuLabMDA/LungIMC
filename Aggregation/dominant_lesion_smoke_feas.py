# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, pickle
from datetime import datetime
import pandas as pd
import numpy as np
from scipy import stats


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])                        
    parser.add_argument("--aggregation_dir",        type=str,       default="Aggregation")
    parser.add_argument("--pval_thresh",            type=float,     default=0.05)   
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)
    lesion_avg_fea_path = os.path.join(slide_agg_dir, "lesion_avg_feas.csv")
    lesion_fea_df = pd.read_csv(lesion_avg_fea_path)

    # normalize all features between 0.0 to 1.0
    lesion_fea_columns = [ele for ele in lesion_fea_df.columns.tolist()]
    lesion_fea_names = lesion_fea_columns[3:]
    for fea in lesion_fea_names:
        lesion_fea_df[fea] = (lesion_fea_df[fea] + 3.0) / 6.0
            
    # collect stage & smoke
    lesion_stages = lesion_fea_df["LesionStage"]
    lesion_smokes = lesion_fea_df["SmokeStatus"]
    aah_never, ais_never, mia_never, adc_never = [], [], [], []
    aah_heavy, ais_heavy, mia_heavy, adc_heavy = [], [], [], []
    for ind in np.arange(len(lesion_stages)):
        if lesion_stages[ind] == "AAH" and lesion_smokes[ind] == "Never":
            aah_never.append(ind)
        elif lesion_stages[ind] == "AIS" and lesion_smokes[ind] == "Never":
            ais_never.append(ind)            
        elif lesion_stages[ind] == "MIA" and lesion_smokes[ind] == "Never":
            mia_never.append(ind)
        elif lesion_stages[ind] == "ADC" and lesion_smokes[ind] == "Never":
            adc_never.append(ind) 
        elif lesion_stages[ind] == "AAH" and lesion_smokes[ind] == "Heavy":
            aah_heavy.append(ind)
        elif lesion_stages[ind] == "AIS" and lesion_smokes[ind] == "Heavy":
            ais_heavy.append(ind)            
        elif lesion_stages[ind] == "MIA" and lesion_smokes[ind] == "Heavy":
            mia_heavy.append(ind)
        elif lesion_stages[ind] == "ADC" and lesion_smokes[ind] == "Heavy":
            adc_heavy.append(ind)
        else:
            print("Unkown stage: {} and smoke {}".format(lesion_stages[ind], lesion_smokes[ind]))
    print("AAH Never: {}".format(len(aah_never)))
    print("AIS Never: {}".format(len(ais_never)))
    print("MIA Never: {}".format(len(mia_never)))
    print("ADC Never: {}".format(len(adc_never)))
    print("AAH Heavy: {}".format(len(aah_heavy)))
    print("AIS Heavy: {}".format(len(ais_heavy)))
    print("MIA Heavy: {}".format(len(mia_heavy)))
    print("ADC Heavy: {}".format(len(adc_heavy)))

    # feature analysis
    sig0_feas, sig1_feas, sig2_feas, sig3_feas, sig4_feas = [], [], [], [], []
    sig0_directions, sig1_directions, sig2_directions, sig3_directions, sig4_directions = [], [], [], [], []
    for cur_fea_name in lesion_fea_names:
        cur_fea_lst = [ele for ele in lesion_fea_df[cur_fea_name].tolist()]
        sig_num = 0
        sig_directions = []
        aah_never_feas = [cur_fea_lst[ele] for ele in aah_never]
        aah_heavy_feas = [cur_fea_lst[ele] for ele in aah_heavy]
        aah_p = stats.ttest_ind(aah_never_feas, aah_heavy_feas)
        if aah_p.pvalue < args.pval_thresh:
            sig_num += 1
            if aah_p.statistic > 0:
                sig_directions.append("-AAH")
            else:
                sig_directions.append("+AAH")
        ais_never_feas = [cur_fea_lst[ele] for ele in ais_never]
        ais_heavy_feas = [cur_fea_lst[ele] for ele in ais_heavy]
        ais_p = stats.ttest_ind(ais_never_feas, ais_heavy_feas)
        if ais_p.pvalue < args.pval_thresh:
            sig_num += 1
            if ais_p.statistic > 0:
                sig_directions.append("-AIS")
            else:
                sig_directions.append("+AIS")
        mia_never_feas = [cur_fea_lst[ele] for ele in mia_never]
        mia_heavy_feas = [cur_fea_lst[ele] for ele in mia_heavy]
        mia_p = stats.ttest_ind(mia_never_feas, mia_heavy_feas)
        if mia_p.pvalue < args.pval_thresh:
            sig_num += 1
            if mia_p.statistic > 0:
                sig_directions.append("-MIA")
            else:
                sig_directions.append("+MIA")
        adc_never_feas = [cur_fea_lst[ele] for ele in adc_never]
        adc_heavy_feas = [cur_fea_lst[ele] for ele in adc_heavy]
        adc_p = stats.ttest_ind(adc_never_feas, adc_heavy_feas)
        if adc_p.pvalue < args.pval_thresh:
            sig_num += 1
            if adc_p.statistic > 0:
                sig_directions.append("-ADC")
            else:
                sig_directions.append("+ADC")

        # stratify feature
        if sig_num == 0:
            sig0_feas.append(cur_fea_name)
            sig0_directions.append(sig_directions)
        if sig_num == 1:
            sig1_feas.append(cur_fea_name)
            sig1_directions.append(sig_directions)
        if sig_num == 2:
            sig2_feas.append(cur_fea_name)
            sig2_directions.append(sig_directions)
        if sig_num == 3:
            sig3_feas.append(cur_fea_name)
            sig3_directions.append(sig_directions)
        if sig_num == 4:
            sig4_feas.append(cur_fea_name) 
            sig4_directions.append(sig_directions)

    print("Sig0 has {}".format(len(sig0_feas)))     
    print("Sig1 has {}".format(len(sig1_feas)))     
    print("Sig2 has {}".format(len(sig2_feas)))
    print("Sig3 has {}".format(len(sig3_feas)))   
    print("Sig4 has {}".format(len(sig4_feas)))      

    print("----Sig3 features:") 
    for fea, direction in zip(sig3_feas, sig3_directions):
        print("{}:{}".format(fea, direction))            