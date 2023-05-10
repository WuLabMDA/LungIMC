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

    # set directory
    dataset_dir = os.path.join(args.data_root, args.data_set)
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)

    # load aggregated feature
    lesion_roi_fea_path = os.path.join(slide_agg_dir, "lesion_roi_feas.xlsx")
    roi_fea_df = pd.read_excel(lesion_roi_fea_path)

    # create empty slide dataframe
    slide_column_lst = ["LesionID", "LesionStage", "SmokeStatus"]
    roi_fea_columns = [ele for ele in roi_fea_df.columns.tolist()]
    roi_fea_columns = roi_fea_columns[3:]
    slide_column_lst.extend(roi_fea_columns)
    slide_df = pd.DataFrame(columns=slide_column_lst)    

    
    # for cur_lesion in lesion_roi_dict.keys():
    #     row_val_lst = [cur_lesion, ]
    #     roi_names = [ele for ele in lesion_roi_dict[cur_lesion].keys()]
    #     cur_lesion_df = roi_fea_df[roi_fea_df["ROI_ID"].isin(roi_names)]
    #     if len(cur_lesion_df) == 0:
    #         print("No lesion detected in {}".format(cur_lesion))
    #         continue
    #     # add stage
    #     cur_stages = cur_lesion_df["ROI_Stage"].tolist()
    #     if len(set(cur_stages)) != 1:
    #         print("Multiple stages in {}".format(cur_lesion))
    #     else:
    #         row_val_lst.append(cur_stages[0])
    #     # add smoke
    #     cur_smokes = cur_lesion_df["SmokeStatus"].tolist()
    #     if len(set(cur_smokes)) != 1:
    #         print("Multiple smokes in {}".format(cur_lesion))
    #     else:
    #         row_val_lst.append(cur_smokes[0])
    #     lesion_fea_df = cur_lesion_df.iloc[:, 3:]
    #     for cur_fea in roi_fea_columns:
    #         row_val_lst.append(np.mean(cur_lesion_df[cur_fea].tolist()))
    #     slide_df.loc[len(slide_df.index)] = row_val_lst
    
    # # Check slide dataframe information 
    # lesion_fea_path = os.path.join(slide_agg_dir, "lesion_avg_feas.xlsx")
    # slide_df.to_excel(lesion_fea_path, index = False)