# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse, pytz
import numpy as np
from datetime import datetime
import warnings


def set_args():
    parser = argparse.ArgumentParser(description = "IMC Cell Feature Extraction")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="Study", choices = ["Study", "Tonsil"])
    parser.add_argument("--merge_dir",              type=str,       default="MergeSeg")
    parser.add_argument("--denoise_dir",            type=str,       default="Denoise")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args() 

    # prepare directory
    root_dir = os.path.join(args.data_root, args.data_type + "Processing")
    roi_raw_dir = os.path.join(root_dir, args.merge_dir, args.denoise_dir, "RawROIs")
    roi_list = [ele for ele in os.listdir(roi_raw_dir) if os.path.isdir(os.path.join(roi_raw_dir, ele))]
    for ind, cur_roi in enumerate(roi_list):
        if (ind + 1) % 100 == 0:
            print("Merge {}/{}".format(ind+1, len(roi_list)))
        src_roi_dir = os.path.join(roi_raw_dir, cur_roi)
        src_stain_path = os.path.join(src_roi_dir, "Ir191_193.tiff")
        dst_stain_path = os.path.join(src_roi_dir, "Ir191.tiff")
        os.rename(src_stain_path, dst_stain_path)