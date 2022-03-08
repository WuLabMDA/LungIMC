# -*- coding: utf-8 -*-

import os, sys
import argparse, shutil


def set_args():
    parser = argparse.ArgumentParser(description = "Merge Cell Segmentation")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="Study", choices = ["Study", "Tonsil"])
    parser.add_argument("--vis_dir",                type=str,       default="VisROI")
    parser.add_argument("--merge_dir",              type=str,       default="VisMerge")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    roi_dir_root = os.path.join(args.data_root, args.data_type + "Processing")
    vis_dir_root = os.path.join(roi_dir_root, args.vis_dir)
    # prepare folders for merge files
    merge_dir_root = os.path.join(roi_dir_root, args.merge_dir)
    if os.path.exists(merge_dir_root):
        shutil.rmtree(merge_dir_root)
    merge_raw_dir = os.path.join(merge_dir_root, "Raw")
    os.makedirs(merge_raw_dir)
    merge_seg_dir = os.path.join(merge_dir_root, "Seg")
    os.makedirs(merge_seg_dir)

    # merge
    p_list = sorted([ele for ele in os.listdir(vis_dir_root) if os.path.isdir(os.path.join(vis_dir_root, ele))])
    for ind, cur_p in enumerate(p_list):
        print("Merge {}/{}".format(ind+1, len(p_list)))
        cur_p_name = cur_p[:-4]
        src_p_dir = os.path.join(vis_dir_root, cur_p)
        raw_img_list = sorted([ele for ele in os.listdir(src_p_dir) if ele.endswith("_raw.png")])
        for cur_img_name in raw_img_list:
            roi_name = cur_img_name[:-8]
            src_raw_path = os.path.join(src_p_dir, roi_name + "_raw.png")
            src_seg_path = os.path.join(src_p_dir, roi_name + "_seg.png")
            dst_roi_name = roi_name[:roi_name.find("_")]
            dst_raw_path = os.path.join(merge_raw_dir, cur_p_name + dst_roi_name + ".png")
            dst_seg_path = os.path.join(merge_seg_dir, cur_p_name + dst_roi_name + ".png")
            shutil.copyfile(src_raw_path, dst_raw_path)
            shutil.copyfile(src_seg_path, dst_seg_path)
