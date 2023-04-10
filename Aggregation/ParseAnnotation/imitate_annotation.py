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
    ann_slide_lst = [ele for ele in slide_lst if os.path.exists(os.path.join(annotation_dir, ele + ".xml"))]
    # print("There are {} annotated slides.".format(len(ann_slide_lst)))

    for ind, cur_slide in enumerate(ann_slide_lst):
        print("Imitate on {} {}/{}".format(cur_slide, ind+1, len(ann_slide_lst)))
        # create background image
        cur_slide_path = os.path.join(annotation_dir, cur_slide + ".svs")
        cur_slide_head = openslide.OpenSlide(cur_slide_path)
        slide_width, slide_height = cur_slide_head.dimensions
        imitate_width = int(slide_width / args.reduction_ratio)
        imitate_height = int(slide_height / args.reduction_ratio)
        imitate_img = np.zeros((imitate_height, imitate_width, 3), dtype=np.uint8)
        # load annotations
        cur_annotation_path = os.path.join(annotation_dir, cur_slide + ".xml")
        lesion_vertices, roi_anno_dict = parse_imagescope_annotations(cur_annotation_path)
        lesion_status = len(lesion_vertices) > 0
        if lesion_status:
            # draw lesion
            lesion_cv_cnt = np.expand_dims((lesion_vertices / args.reduction_ratio).astype(np.int32), axis=1)
            cv2.drawContours(imitate_img, [lesion_cv_cnt, ], 0, (255, 0, 0), 3)
        # draw rois
        # print("--- with {} ROIs".format(len(roi_anno_dict)))        
        for roi_name, roi_cnt in roi_anno_dict.items():
            reduction_cnt = roi_cnt / args.reduction_ratio
            mean_x, mean_y = np.mean(reduction_cnt, axis=0)
            cv2.putText(imitate_img, roi_name, (int(mean_x), int(mean_y)), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 255, 0), 2,  cv2.LINE_AA)
            roi_cv_cnt = np.expand_dims(reduction_cnt.astype(np.int32), axis=1)
            cv2.drawContours(imitate_img, [roi_cv_cnt, ], 0, (0, 0, 255), 2)

        # Save Imitations
        if lesion_status:
            anno_imitation_dir = os.path.join(imitation_dir, "Tumor")
        else:
            anno_imitation_dir = os.path.join(imitation_dir, "Normal")
        if not os.path.exists(anno_imitation_dir):
            os.makedirs(anno_imitation_dir)
        # save imitation image
        anno_imitation_path = os.path.join(anno_imitation_dir, cur_slide + ".png")
        io.imsave(anno_imitation_path, imitate_img)