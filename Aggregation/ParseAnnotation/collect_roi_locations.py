# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, pickle
from datetime import datetime
import openslide
import pandas as pd
import numpy as np
from skimage import io
from pycontour import coor_transform, rela 

from annotation_utils import parse_imagescope_annotations


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])  
    parser.add_argument("--aggregation_dir",        type=str,       default="Aggregation")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    demo_root_dir = os.path.join(dataset_dir, "NatureFigures", "Fig03")
    annotation_dir = os.path.join(demo_root_dir, "Slidesannotation")
    imitation_dir = os.path.join(demo_root_dir, "Imitations")
    if not os.path.exists(imitation_dir):
        os.makedirs(imitation_dir)   
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)
    if not os.path.exists(slide_agg_dir):
        os.makedirs(slide_agg_dir) 

    slide_lst = [os.path.splitext(ele)[0] for ele in os.listdir(annotation_dir) if ele.endswith(".svs")]
    print("Number of slides: {}".format(len(slide_lst)))
    ann_slide_lst = [ele for ele in slide_lst if os.path.exists(os.path.join(annotation_dir, ele + ".xml"))]

    lesion_roi_dict = {}
    for ind, cur_slide in enumerate(ann_slide_lst):
        print("Collect on {} {}/{}".format(cur_slide, ind+1, len(ann_slide_lst)))
        # load annotations
        cur_annotation_path = os.path.join(annotation_dir, cur_slide + ".xml")
        lesion_vertices, roi_anno_dict = parse_imagescope_annotations(cur_annotation_path)
        lesion_status = len(lesion_vertices) > 0
        if lesion_status:
            # draw lesion
            lesion_dict = {}
            lesion_cnt_arr = np.transpose(lesion_vertices).astype(np.int32)
            lesion_cnt_arr = coor_transform.swap_wh(lesion_cnt_arr)
            lesion_dict["lesion_cnt"] = lesion_cnt_arr
            for roi_name, roi_cnt in roi_anno_dict.items():
                reduction_cnt = roi_cnt
                mean_x, mean_y = np.mean(reduction_cnt, axis=0)
                mean_x, mean_y = int(mean_x), int(mean_y)
                in_lesion_status = rela.point_in_contour(lesion_cnt_arr, mean_y, mean_x)
                lesion_dict[roi_name] = {"X":mean_x, "Y": mean_y}
            lesion_roi_dict[cur_slide] = lesion_dict
    print("There are {} lesions with annotated ROIs.".format(len(lesion_roi_dict)))
    # save information to pkl
    lesion_roi_loc_path = os.path.join(slide_agg_dir, "lesion_roi_loc.pkl")
    with open(lesion_roi_loc_path, 'wb') as handle:
        pickle.dump(lesion_roi_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)