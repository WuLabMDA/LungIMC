# -*- coding: utf-8 -*-

import os, sys
import argparse


def set_args():
    parser = argparse.ArgumentParser(description = "Filter Segmented Cells")
    parser.add_argument("--data_root",        type=str,       default="/Data")
    parser.add_argument("--data_type",        type=str,       default="Study", choices = ["Study", "Tonsil"])

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    roi_root_dir = os.path.join(args.data_root, args.data_type + "SlidesROIs")
    slide_list = sorted([ele for ele in os.listdir(roi_root_dir) if os.path.isdir(os.path.join(roi_root_dir, ele))])
    ttl_roi_num = 0
    for cur_slide in slide_list:
        cur_roi_dir = os.path.join(roi_root_dir, cur_slide)
        roi_list = sorted([ele for ele in os.listdir(cur_roi_dir) if os.path.isdir(os.path.join(cur_roi_dir, ele))])
        ttl_roi_num += len(roi_list)
    print("{} has {} ROIs.".format(args.data_type, ttl_roi_num))
