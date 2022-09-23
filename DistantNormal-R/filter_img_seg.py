# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import pandas as pd


def set_args():
    parser = argparse.ArgumentParser(description = "Filter ROIs")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--raw_steinbock_dir",      type=str,       default="Steinbock")
    parser.add_argument("--dst_steinbock_dir",      type=str,       default="SteinbockDistantNormal")    

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()
    raw_steinbock_dir = os.path.join(args.data_root, args.data_type, args.raw_steinbock_dir)
    raw_img_dir = os.path.join(raw_steinbock_dir, "img")
    raw_seg_dir = os.path.join(raw_steinbock_dir, "masks_deepcell")

    # locate all rois
    dst_steinbock_dir = os.path.join(args.data_root, args.data_type, args.dst_steinbock_dir)
    dst_roi_info_path = os.path.join(dst_steinbock_dir, "StudyROI_Info.xlsx")
    dst_roi_df = pd.read_excel(dst_roi_info_path, index_col=None, sheet_name="Sheet1")
    dst_roi_list = dst_roi_df["ROI_ID"].tolist()

    # copy to destination
    dst_img_dir = os.path.join(dst_steinbock_dir, "img")
    if os.path.exists(dst_img_dir):
        shutil.rmtree(dst_img_dir)
    os.makedirs(dst_img_dir)
    dst_seg_dir = os.path.join(dst_steinbock_dir, "masks_deepcell")
    if os.path.exists(dst_seg_dir):
        shutil.rmtree(dst_seg_dir)
    os.makedirs(dst_seg_dir)

    for ind, cur_roi in enumerate(dst_roi_list):
        print("Copy {}/{}".format(ind+1, len(dst_roi_list)))
        raw_img_path = os.path.join(raw_img_dir, cur_roi + ".tiff")
        raw_seg_path = os.path.join(raw_seg_dir, cur_roi + ".tiff")
        shutil.copy(raw_img_path, dst_img_dir)
        shutil.copy(raw_seg_path, dst_seg_dir)
