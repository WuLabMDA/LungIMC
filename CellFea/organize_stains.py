# -*- coding: utf-8 -*-

import os, sys
import argparse
import numpy as np
from scipy.io import loadmat
import pytz
from datetime import datetime
from skimage import io
import tifffile


def set_args():
    parser = argparse.ArgumentParser(description = "Organize stain mat files")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="Study", choices = ["Study", "Tonsil"])
    parser.add_argument("--mat_dir",                type=str,       default="MatROI")
    parser.add_argument("--merge_dir",              type=str,       default="MergeSeg")

    args = parser.parse_args()
    return args



if __name__ == "__main__":
    args = set_args()
    mat_root = os.path.join(args.data_root, args.data_type + "Processing", args.mat_dir)
    stain_root = os.path.join(args.data_root, args.data_type + "Processing", args.merge_dir, "Stains")
    if not os.path.exists(stain_root):
        os.makedirs(stain_root)

    p_list = sorted([ele for ele in os.listdir(mat_root) if os.path.isdir(os.path.join(mat_root, ele))])
    for ind, cur_p in enumerate(p_list):
        print("Organize {}/{}".format(ind+1, len(p_list)))
        print("Current time: " + datetime.now(pytz.timezone('America/Chicago')).strftime("%H:%M:%S"))
        cur_p_dir = os.path.join(mat_root, cur_p)
        roi_list = sorted([ele for ele in os.listdir(cur_p_dir) if os.path.isdir(os.path.join(cur_p_dir, ele))])
        for cur_roi in roi_list:
            dst_roi_name = cur_roi[:cur_roi.find("_")]
            stain_list = [ele for ele in os.listdir(os.path.join(cur_p_dir, cur_roi)) if ele.endswith(".mat")]
            if len(stain_list) != 36:
                print("{}-{}".format(cur_p, dst_roi_name))



        # slide_name = cur_grp[:-5]
        # sub_dirs = sorted([ele for ele in os.listdir(os.path.join(grp_root_dir, cur_grp))])
        # for cur_sub in sub_dirs:
        #     src_img_dir = os.path.join(grp_root_dir, cur_grp, cur_sub)
        #     src_seg_dir = os.path.join(seg_root_dir, cur_grp)
        #     slide_roi_name = "-".join([slide_name, cur_sub])
        #     seg_path = os.path.join(src_seg_dir, cur_sub + "_feature_0.tif")
        #     seg_mask = io.imread(seg_path).astype(np.int32)
        #     # Ignore those ROIs with limited cell numbers
        #     if len(np.unique(seg_mask)) < args.min_cell_num:
        #         continue
        #     dst_roi_dir = os.path.join(segstain_root, slide_roi_name)
        #     if not os.path.exists(dst_roi_dir):
        #         os.makedirs(dst_roi_dir)
        #     # save cell mask
        #     stain_seg_path = os.path.join(dst_roi_dir, slide_roi_name + "_mask.tiff")
        #     tifffile.imwrite(stain_seg_path, seg_mask)
        #     # save marekrs files
        #     marker_list = sorted([os.path.splitext(ele)[0] for ele in os.listdir(src_img_dir) if ele.endswith(".mat")])
        #     for cur_marker in marker_list:
        #         cur_marker_path = os.path.join(src_img_dir, cur_marker + ".mat")
        #         cur_mat_data = loadmat(cur_marker_path)
        #         cur_marker_img = cur_mat_data["stain_img"].astype(np.uint16)
        #         stain_marker_path = os.path.join(dst_roi_dir, cur_marker + ".tiff")
        #         tifffile.imwrite(stain_marker_path, cur_marker_img)
