# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, pickle
from datetime import datetime
import pandas as pd
import numpy as np
import cv2
from pycontour import cv2_transform 

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
    lesion_roi_loc_path = os.path.join(slide_agg_dir, "lesion_roi_loc.pkl")
    lesion_roi_dict = None
    with open(lesion_roi_loc_path, 'rb') as handle:
        lesion_roi_dict = pickle.load(handle)        
    lesion_names = [ele for ele in lesion_roi_dict.keys()]
    print("There are {} annotated leision.".format(len(lesion_names)))
    lesion_roi_weight_dict = {}
    ttl_roi_num = 0
    inside_roi_num = 0
    for cur_lesion in lesion_names:
        roi_weight_dict = {}
        lesion_dict = lesion_roi_dict[cur_lesion]
        lesion_np_cnt = lesion_dict["lesion_cnt"]
        lesion_cv_cnt = cv2_transform.np_arr_to_cv_cnt(lesion_np_cnt)
        lesion_M = cv2.moments(lesion_cv_cnt)
        lesion_cx = int(lesion_M["m10"] / lesion_M["m00"])
        lesion_cy = int(lesion_M["m01"] / lesion_M["m00"])
        cur_roi_lst = [ele for ele in lesion_dict.keys() if ele.startswith("ROI")]
        for cur_roi in cur_roi_lst:
            cur_roi_name = cur_lesion + "-" + cur_roi
            ttl_roi_num += 1
            roi_cx = lesion_dict[cur_roi]["X"]
            roi_cy = lesion_dict[cur_roi]["Y"]
            border_dist = cv2.pointPolygonTest(lesion_cv_cnt, (roi_cx, roi_cy), measureDist=True)
            if border_dist > 0:
                inside_roi_num += 1
                core_dist = np.sqrt(np.power(lesion_cx-roi_cx, 2) + np.power(lesion_cy-roi_cy, 2))
                roi_weight_dict[cur_roi_name] = {"BorderDist": int(border_dist), "CoreDist": int(core_dist)}
        lesion_roi_weight_dict[cur_lesion] = roi_weight_dict
    print("ROI total number is: {}".format(ttl_roi_num))
    print("ROI inside number is: {}".format(inside_roi_num))
    # save information to json
    lesion_roi_dist_path = os.path.join(slide_agg_dir, "lesion_roi_dist.json")
    with open(lesion_roi_dist_path, 'w') as fp:
        json.dump(lesion_roi_weight_dict, fp)