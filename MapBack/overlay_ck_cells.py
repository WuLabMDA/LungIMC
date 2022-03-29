# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse, pickle
import deepdish as dd
import numpy as np
import tifffile, cv2
from skimage import io


def set_args():
    parser = argparse.ArgumentParser(description = "IMC Cell Filtering")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="Study", choices = ["Study", "Tonsil"])
    parser.add_argument("--merge_dir",              type=str,       default="MergeSeg")
    parser.add_argument("--vis_dir",                type=str,       default="Vis")
    parser.add_argument("--fea_option",             type=str,       default="Transform", choices = ["Transform", "SelfCorrect", "ControlCorrect"])
    parser.add_argument("--seed",                   type=int,       default=1234)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()
    # prepare directory
    merge_dir_root = os.path.join(args.data_root, args.data_type + "Processing", args.merge_dir)
    stain_dir = os.path.join(merge_dir_root, "Stains")
    stain_roi_list = [ele for ele in os.listdir(stain_dir) if os.path.isdir(os.path.join(stain_dir, ele))]
    print("There are {} ROIs.".format(len(stain_roi_list)))
    cnt_dir = os.path.join(merge_dir_root, "Contour")
    for ele in stain_roi_list:
        cur_cnt_path = os.path.join(cnt_dir, ele + ".h5")
        if not os.path.exists(cur_cnt_path):
            print("{} no contours.".format(ele))

    # data_root = os.path.join(args.data_root, args.tissue_type)
    # cellfea_dir = os.path.join(data_root, "CellFea")
    # cnt_dir = os.path.join(cellfea_dir, "CellCNTs")
    # overlay_dir = os.path.join(cellfea_dir, "OverlayCK")
    # if not os.path.exists(overlay_dir):
    #     os.makedirs(overlay_dir)
