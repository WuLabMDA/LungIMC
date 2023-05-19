# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, pickle
from datetime import datetime
import pandas as pd
import numpy as np


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])  
    parser.add_argument("--feature_dir",            type=str,       default="FeatureAnalysis")  
    parser.add_argument("--aggregation_dir",        type=str,       default="Aggregation")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    # load location information
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)

    # save information to json
    lesion_roi_dist_path = os.path.join(slide_agg_dir, "lesion_roi_dist.json")
    lesion_roi_weight_dict =  None
    with open(lesion_roi_dist_path) as fp:
        lesion_roi_weight_dict = json.load(fp)
    
    # load slide dataframe information 
    lesion_fea_path = os.path.join(slide_agg_dir, "lesion_avg_feas.csv")
    slide_df = pd.read_csv(lesion_fea_path)
    slide_df = slide_df[slide_df["LesionStage"] != "Normal"]
    lesion_id_lst = slide_df["LesionID"].tolist()
    lesion_stage_lst = slide_df["LesionStage"].tolist()
    stage_lesion_dict = {}
    for cur_id, cur_stage in zip(lesion_id_lst, lesion_stage_lst):
        if cur_stage not in stage_lesion_dict:
            stage_lesion_dict[cur_stage] = [cur_id, ]
        else:
            stage_lesion_dict[cur_stage].append(cur_id)

    stage_lst = ["AAH", "AIS", "MIA", "ADC"]
    for cur_stage in stage_lst:
        lesion_lst = stage_lesion_dict[cur_stage]
        roi_num_lst = []
        for cur_lesion in lesion_lst:
            if cur_lesion in lesion_roi_weight_dict:
                roi_num_lst.append(len(lesion_roi_weight_dict[cur_lesion]))
        print("Stage {}".format(cur_stage))
        print("--Average  No. of ROIs of lesion: {:.2f}".format(np.mean(roi_num_lst)))
        print("--Minmum   No. of ROIs of lesion: {}".format(np.min(roi_num_lst)))
        print("--Maxmium  No. of ROIs of lesion: {}".format(np.max(roi_num_lst)))