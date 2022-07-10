# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse, pytz
import numpy as np
from skimage import io
from datetime import datetime
import warnings


def set_args():
    parser = argparse.ArgumentParser(description = "IMC Cell Feature Extraction")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="Study", choices = ["Study", "Tonsil"])
    parser.add_argument("--merge_dir",              type=str,       default="MergeSeg")
    parser.add_argument("--stain_dir",              type=str,       default="Stains")
    parser.add_argument("--spillover_dir",          type=str,       default="Spillover")
    parser.add_argument("--denoise_dir",            type=str,       default="Denoise")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args() 

    # prepare directory
    root_dir = os.path.join(args.data_root, args.data_type + "Processing")
    roi_spillover_dir = os.path.join(root_dir, args.merge_dir, args.denoise_dir, "ROIs")
    if os.path.exists(roi_spillover_dir):
        shutil.rmtree(roi_spillover_dir)
    os.makedirs(roi_spillover_dir)

    stain_root = os.path.join(root_dir, args.merge_dir, args.stain_dir)
    roi_list = [ele for ele in os.listdir(stain_root) if os.path.isdir(os.path.join(stain_root, ele))]
    for ind, cur_roi in enumerate(roi_list):
        if (ind + 1) % 100 == 0:
            print("Merge {}/{}".format(ind+1, len(roi_list)))
        src_roi_dir = os.path.join(stain_root, cur_roi)
        dst_roi_dir = os.path.join(roi_spillover_dir, cur_roi)
        if not os.path.exists(dst_roi_dir):
            os.makedirs(dst_roi_dir)
        src_stain_path = os.path.join(src_roi_dir, "Ir191_193.tiff")
        dst_stain_path = os.path.join(dst_roi_dir, "Ir191.tiff")
        shutil.copyfile(src_stain_path, dst_stain_path)