# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, pickle, math
from datetime import datetime
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats import multitest


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])  
    parser.add_argument("--aggregation_dir",        type=str,       default="Aggregation")
    parser.add_argument("--ith_dir",                type=str,       default="ITH")
    parser.add_argument("--pval_thresh",            type=float,     default=0.05)

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    # set directory
    dataset_dir = os.path.join(args.data_root, args.data_set)
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)
    lesion_ith_dir = os.path.join(slide_agg_dir, args.ith_dir)
    lesion_ith_path = os.path.join(lesion_ith_dir, "lesion_ith_per_fea.csv")
    slide_ith_df = pd.read_csv(lesion_ith_path)

    # obtain feature names
    lesion_fea_columns = [ele for ele in slide_ith_df.columns.tolist()]
    lesion_fea_names = lesion_fea_columns[3:]
    
    # iterate each feature
    lesion_stage_lst = [ele for ele in slide_ith_df["LesionStage"].tolist()]
    for cur_fea in lesion_fea_names:
        aah_lst, ais_lst, mia_lst, adc_lst = [], [], [], []
        fea_ith_lst = [ele for ele in slide_ith_df[cur_fea].tolist()]
        for cur_stage, cur_fea in zip(lesion_stage_lst, fea_ith_lst):
            if cur_stage == "AAH":
                aah_lst.append(cur_fea)
            elif cur_stage == "AIS":
                ais_lst.append(cur_fea)
            elif cur_stage == "MIA":
                ais_lst.append(cur_fea)
            elif cur_stage == "ADC":
                ais_lst.append(cur_fea)
            else:
                sys.exit("Unknow stage: {}".format(cur_stage))
        import pdb; pdb.set_trace()
        