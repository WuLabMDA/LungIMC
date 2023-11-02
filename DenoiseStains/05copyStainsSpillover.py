# -*- coding: utf-8 -*-

import os, sys
import numpy as np
import argparse, pytz, shutil
from datetime import datetime
from scipy.io import loadmat
import tifffile


def set_args():
    parser = argparse.ArgumentParser(description = "Organize stain 191 mat files")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--mat_dir",                type=str,       default="GroupROI")
    parser.add_argument("--denoise_dir",            type=str,       default="Denoise")
    parser.add_argument("--dst_dir",                type=str,       default="SpilloverROIs")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()
    mat_root = os.path.join(args.data_root, args.data_type, args.mat_dir)
    stain_root = os.path.join(args.data_root, args.data_type, args.denoise_dir, args.dst_dir)
    if not os.path.exists(stain_root):
        os.makedirs(stain_root)
    # Organize
    p_list = sorted([ele for ele in os.listdir(mat_root) if os.path.isdir(os.path.join(mat_root, ele))])
    for ind, cur_p in enumerate(p_list):
        dst_p_name = cur_p[:-4]
        cur_time_str = datetime.now(pytz.timezone('America/Chicago')).strftime("%H:%M:%S")
        print("Organize {:3d}/{:3d} Current time: {}".format(ind+1, len(p_list), cur_time_str))
        cur_p_dir = os.path.join(mat_root, cur_p)
        roi_list = sorted([ele for ele in os.listdir(cur_p_dir) if os.path.isdir(os.path.join(cur_p_dir, ele))])
        for cur_roi in roi_list:
            dst_roi_name = cur_roi[:cur_roi.find("_")]
            roi_dst_dir = os.path.join(stain_root, dst_p_name + dst_roi_name)
            if not os.path.exists(roi_dst_dir):
               os.makedirs(roi_dst_dir)
            src_stain_path = os.path.join(cur_p_dir, cur_roi, "191Ir.mat")
            cur_mat_data = loadmat(src_stain_path)
            cur_marker_img = cur_mat_data["stain_img"].astype(np.uint16)
            dst_stain_path = os.path.join(roi_dst_dir, "Ir191.tiff")
            tifffile.imwrite(dst_stain_path, cur_marker_img)