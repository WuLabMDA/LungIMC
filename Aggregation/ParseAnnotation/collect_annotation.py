# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz
from datetime import datetime
import openslide
import pandas as pd
import numpy as np
from skimage import io
import cv2
from pycontour import coor_transform, rela 

from annotation_utils import parse_imagescope_annotations


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])  
    parser.add_argument("--reduction_ratio",        type=float,     default=10.0) 
                    
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

    slide_lst = [os.path.splitext(ele)[0] for ele in os.listdir(annotation_dir) if ele.endswith(".svs")]
    print("Number of slides: {}".format(len(slide_lst)))
    ann_slide_lst = [ele for ele in slide_lst if os.path.exists(os.path.join(annotation_dir, ele + ".xml"))]

    slide_num = 0
    for ind, cur_slide in enumerate(ann_slide_lst):
        print("Collect on {} {}/{}".format(cur_slide, ind+1, len(ann_slide_lst)))
        # load annotations
        cur_annotation_path = os.path.join(annotation_dir, cur_slide + ".xml")
        lesion_vertices, roi_anno_dict = parse_imagescope_annotations(cur_annotation_path)
        lesion_status = len(lesion_vertices) > 0
        if lesion_status:
            # draw lesion
            lesion_cnt_arr = lesion_cnt_arr.swap_wh(lesion_vertices).astype(np.int32)
            for roi_name, roi_cnt in roi_anno_dict.items():
                reduction_cnt = roi_cnt
                mean_x, mean_y = np.mean(reduction_cnt, axis=0)
                in_lesion_status = rela.point_in_contour(lesion_cnt_arr, mean_y, mean_x)
                import pdb; pdb.set_trace()

