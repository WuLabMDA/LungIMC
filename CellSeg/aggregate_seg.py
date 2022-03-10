# -*- coding: utf-8 -*-

import os, sys
import argparse, shutil


def set_args():
    parser = argparse.ArgumentParser(description = "Merge Cell Segmentation")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="Study", choices = ["Study", "Tonsil"])
    parser.add_argument("--seg_dir",                type=str,       default="SegResults")
    parser.add_argument("--merge_dir",              type=str,       default="MergeSeg")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    roi_dir_root = os.path.join(args.data_root, args.data_type + "Processing")
    seg_dir_root = os.path.join(roi_dir_root, args.seg_dir)
    # prepare folders for merge files
    merge_dir_root = os.path.join(roi_dir_root, args.merge_dir)
    merge_mask_dir = os.path.join(merge_dir_root, "Mask")
    if os.path.exists(merge_mask_dir):
        shutil.rmtree(merge_mask_dir)
    os.makedirs(merge_mask_dir)

    # merge
    p_list = sorted([ele for ele in os.listdir(seg_dir_root) if os.path.isdir(os.path.join(seg_dir_root, ele))])
    for ind, cur_p in enumerate(p_list):
        print("Merge {}/{}".format(ind+1, len(p_list)))
        cur_p_name = cur_p[:-4]
        src_p_dir = os.path.join(seg_dir_root, cur_p)
        raw_img_list = sorted([ele for ele in os.listdir(src_p_dir) if ele.endswith(".npy")])
        for cur_img_name in raw_img_list:
            roi_name = cur_img_name[:-4]
            src_seg_path = os.path.join(src_p_dir, roi_name + ".npy")
            dst_roi_name = roi_name[:roi_name.find("_")]
            dst_seg_path = os.path.join(merge_mask_dir, cur_p_name + dst_roi_name + ".npy")
            shutil.copyfile(src_seg_path, dst_seg_path)
