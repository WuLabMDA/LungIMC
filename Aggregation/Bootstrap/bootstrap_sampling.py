# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse, pytz, pickle
import json, math, random
from datetime import datetime
import pandas as pd
import numpy as np


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])  
    parser.add_argument("--feature_dir",            type=str,       default="FeatureAnalysis")  
    parser.add_argument("--aggregation_dir",        type=str,       default="Aggregation")
    parser.add_argument("--rand_roi_num",           type=int,       default=50)
    parser.add_argument("--rand_times",             type=int,       default=5)
    parser.add_argument("--rand_seed",              type=int,       default=1234)

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()
    random.seed(args.rand_seed)

    dataset_dir = os.path.join(args.data_root, args.data_set)
    # load location information
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)

    # load aggregated feature
    lesion_roi_fea_path = os.path.join(slide_agg_dir, "lesion_roi_feas.xlsx")
    roi_fea_df = pd.read_excel(lesion_roi_fea_path)
    roi_names = roi_fea_df["ROI_ID"].tolist()

    # create empty slide dataframe
    slide_column_lst = ["LesionID", "LesionStage", "SmokeStatus"]
    roi_fea_columns = [ele for ele in roi_fea_df.columns.tolist()]
    roi_fea_columns = roi_fea_columns[3:]
    slide_column_lst.extend(roi_fea_columns)
    slide_df = pd.DataFrame(columns=slide_column_lst)    

    # save information to json
    lesion_roi_dist_path = os.path.join(slide_agg_dir, "lesion_roi_dist.json")
    lesion_roi_weight_dict =  None
    with open(lesion_roi_dist_path) as fp:
        lesion_roi_weight_dict = json.load(fp)
    
    lesion_lst = [ele for ele in lesion_roi_weight_dict.keys()]
    for cur_lesion in lesion_lst[2:]:
        lesion_roi_lst = [ele for ele in lesion_roi_weight_dict[cur_lesion].keys()]
        for rind in np.arange(args.rand_times):
            rand_rois = random.choices(lesion_roi_lst, k=args.rand_roi_num)
            rand_inds = []
            for ele in rand_rois:
                if ele in roi_names:
                    rand_inds.append(roi_names.index(ele))
            if len(rand_inds) < (args.rand_roi_num - 10):
                print("{} sampling {}".format(cur_lesion, len(rand_inds)))
            rand_lesion_df = roi_fea_df.iloc[rand_inds]
            smoke_lst = rand_lesion_df["SmokeStatus"].tolist()
            if len(set(smoke_lst)) != 1:
                print("Multiple smokes in {}".format(cur_lesion))
                continue  
            cur_smoke = smoke_lst[0]                      
            stage_lst = rand_lesion_df["ROI_Stage"].tolist()
            if len(set(stage_lst)) != 1:
                print("Multiple stages in {}".format(cur_lesion))
                continue
            cur_stage = stage_lst[0]                
            row_val_lst = [cur_lesion, cur_stage, cur_smoke]
            for cur_fea in roi_fea_columns:
                row_val_lst.append(np.mean(rand_lesion_df[cur_fea].tolist()))
            slide_df.loc[len(slide_df.index)] = row_val_lst
    # Check slide dataframe information 
    bootstrap_lesion_fea_path = os.path.join(slide_agg_dir, "lesion_bootstrap_feas.xlsx")
    slide_df.to_excel(bootstrap_lesion_fea_path, index = False)