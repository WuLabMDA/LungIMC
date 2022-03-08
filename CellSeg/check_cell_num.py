# -*- coding: utf-8 -*-

import os, sys
import argparse
import numpy as np


def set_args():
    parser = argparse.ArgumentParser(description = "Filter Segmented Cells")
    parser.add_argument("--data_root",        type=str,       default="/Data")
    parser.add_argument("--data_type",        type=str,       default="Study", choices = ["Study", "Tonsil"])

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    seg_root_dir = os.path.join(args.data_root, args.data_type + "Processing", "SegResults")
    slide_list = sorted([ele for ele in os.listdir(seg_root_dir) if os.path.isdir(os.path.join(seg_root_dir, ele))])

    ttl_roi_num, ttl_cell_num = 0, 0
    for cur_slide in slide_list:
        cur_slide_dir = os.path.join(seg_root_dir, cur_slide)
        roi_list = [ele for ele in os.listdir(cur_slide_dir) if ele.endswith("npy")]
        ttl_roi_num += len(roi_list)
        for cur_roi in roi_list:
            seg_result_path = os.path.join(cur_slide_dir, cur_roi)
            cur_seg = np.load(seg_result_path)
            ttl_cell_num += len(np.unique(cur_seg)) - 1
    print("There are {} ROIs with {} cells.".format(ttl_roi_num, ttl_cell_num))
    avg_cell_roi = ttl_cell_num * 1.0 / ttl_roi_num
    print("On average, there are {:.1f} cells per ROI.".format(avg_cell_roi))
