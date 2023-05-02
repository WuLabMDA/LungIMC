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
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    # load aggregated feature
    feature_root_dir = os.path.join(dataset_dir, args.feature_dir)    
    merge_roi_fea_path = os.path.join(feature_root_dir, "ROI_Fea_Aggregation.csv")
    roi_fea_df = pd.read_csv(merge_roi_fea_path)
    all_roi_ids = [ele for ele in roi_fea_df["ROI_ID"].tolist()]

    # load location information
    demo_root_dir = os.path.join(dataset_dir, "NatureFigures", "Fig03")    
    imitation_dir = os.path.join(demo_root_dir, "Imitations")
    lesion_roi_loc_path = os.path.join(imitation_dir, "lesion_roi_loc.pkl")
    lesion_roi_dict = None
    with open(lesion_roi_loc_path, 'rb') as handle:
        lesion_roi_dict = pickle.load(handle)        
    lesion_names = [ele for ele in lesion_roi_dict.keys()]
    print("There are {} annotated leision.".format(len(lesion_names)))
    # , 'ROI007':