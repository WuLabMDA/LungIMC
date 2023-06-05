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

from annotation_utils import parse_imagescope_annotations


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])  
    parser.add_argument("--aggregation_dir",        type=str,       default="Aggregation")
    parser.add_argument("--ith_dir",                type=str,       default="ITH")
    parser.add_argument("--dlevel",                 type=int,       default=2)
                    
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    demo_root_dir = os.path.join(dataset_dir, "NatureFigures", "Fig03")
    annotation_dir = os.path.join(demo_root_dir, "Slidesannotation")
    
    # prepare folders for annotation imitation
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)
    ith_dir = os.path.join(slide_agg_dir, args.ith_dir)
    imitation_dir = os.path.join(ith_dir, "Imitations")
    if not os.path.exists(imitation_dir):
        os.makedirs(imitation_dir)
    slide_lst = [os.path.splitext(ele)[0] for ele in os.listdir(annotation_dir) if ele.endswith(".svs")]
    ann_slide_lst = [ele for ele in slide_lst if os.path.exists(os.path.join(annotation_dir, ele + ".xml"))]

    for ind, cur_slide in enumerate(ann_slide_lst):
        if cur_slide != "2166-1B":
            continue
        print("Imitate on {} {}/{}".format(cur_slide, ind+1, len(ann_slide_lst)))
        # create background image
        cur_slide_path = os.path.join(annotation_dir, cur_slide + ".svs")
        cur_slide_head = openslide.OpenSlide(cur_slide_path)
        slide_width, slide_height = cur_slide_head.dimensions 
        # read image
        cur_slide_img = cur_slide_head.read_region(location=(0,0), level=args.dlevel, size=cur_slide_head.level_dimensions[args.dlevel])
        cur_slide_img = np.asarray(cur_slide_img)[:,:,:3]
        # contour image
        reduction_ratio = cur_slide_head.level_downsamples[args.dlevel]
        imitate_width = int(slide_width / reduction_ratio)
        imitate_height = int(slide_height / reduction_ratio)
        # load annotations
        cur_annotation_path = os.path.join(annotation_dir, cur_slide + ".xml")
        lesion_vertices, roi_anno_dict = parse_imagescope_annotations(cur_annotation_path)
        if len(lesion_vertices) > 0:
            lesion_dict = {}
            # print("--- with {} ROIs".format(len(roi_anno_dict)))        
            for roi_name, roi_cnt in roi_anno_dict.items():
                reduction_cnt = roi_cnt / reduction_ratio
                mean_x, mean_y = np.mean(reduction_cnt, axis=0)
                lesion_dict[roi_name] = {"xcoor": int(mean_x), "ycoor": int(mean_y)}
            lesion_dict["cnt"] = [[int(ele[0]), int(ele[1])] for ele in (lesion_vertices / reduction_ratio).tolist()]
            # save imitation image
            slide_img_path = os.path.join(imitation_dir, cur_slide + ".png")
            io.imsave(slide_img_path, cur_slide_img)
            slide_anno_path = os.path.join(imitation_dir, cur_slide + ".json")
            with open(slide_anno_path, 'w') as fp:
                json.dump(lesion_dict, fp)            