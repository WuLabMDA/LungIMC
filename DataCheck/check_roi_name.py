# -*- coding: utf-8 -*-

import os, sys
import argparse


def set_args():
    parser = argparse.ArgumentParser(description = "Filter Segmented Cells")
    parser.add_argument("--data_root",        type=str,       default="/Data")
    parser.add_argument("--data_type",        type=str,       default="Study", choices = ["Study", "Suppl", "Tonsil"])

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    roi_root_dir = os.path.join(args.data_root, args.data_type + "SlidesROIs")
    slide_list = sorted([ele for ele in os.listdir(roi_root_dir) if os.path.isdir(os.path.join(roi_root_dir, ele))])
    for cur_slide in slide_list:
        cur_roi_dir = os.path.join(roi_root_dir, cur_slide)
        roi_list = sorted([ele for ele in os.listdir(cur_roi_dir) if os.path.isdir(os.path.join(cur_roi_dir, ele))])
        for roi_name in roi_list:
            first_loc = roi_name.find("_")
            last_loc = roi_name.rfind("_")
            first_roi_id = roi_name[3:first_loc]
            last_roi_id = roi_name[last_loc+1:]
            if first_roi_id != last_roi_id:
                print("{} {}".format(cur_slide, roi_name))
