# -*- coding: utf-8 -*-

import os, sys
import argparse, shutil


def set_args():
    parser = argparse.ArgumentParser(description = "Merge Cell Segmentation")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--cellseg_dir",            type=str,       default="CellSeg")        
    parser.add_argument("--result_dir",             type=str,       default="SegCellResults")
    parser.add_argument("--merge_dir",              type=str,       default="MergeCellSeg")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    roi_dir_root = os.path.join(args.data_root, args.data_type, args.cellseg_dir)
    seg_dir_root = os.path.join(roi_dir_root, args.result_dir)
    # prepare folders for merge files
    merge_dir_root = os.path.join(roi_dir_root, args.merge_dir)
    merge_raw_dir = os.path.join(merge_dir_root, "Raw")
    if os.path.exists(merge_raw_dir):
        shutil.rmtree(merge_raw_dir)
    os.makedirs(merge_raw_dir)
    merge_overlay_dir = os.path.join(merge_dir_root, "Overlay")
    if os.path.exists(merge_overlay_dir):
        shutil.rmtree(merge_overlay_dir)
    os.makedirs(merge_overlay_dir)
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
            # src filenames
            roi_name = cur_img_name[:-4]
            src_raw_path = os.path.join(src_p_dir, roi_name + "_raw.png")
            src_overlay_path = os.path.join(src_p_dir, roi_name + "_seg.png")
            src_mask_path = os.path.join(src_p_dir, roi_name + ".npy")
            # dst filenames
            dst_raw_path = os.path.join(merge_raw_dir, roi_name + ".png")
            dst_overlay_path = os.path.join(merge_overlay_dir, roi_name + ".png")
            dst_mask_path = os.path.join(merge_mask_dir, roi_name + ".npy")
            # copy
            shutil.copyfile(src_raw_path, dst_raw_path)
            shutil.copyfile(src_overlay_path, dst_overlay_path)
            shutil.copyfile(src_mask_path, dst_mask_path)
