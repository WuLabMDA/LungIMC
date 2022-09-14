# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import pandas as pd


def set_args():
    parser = argparse.ArgumentParser(description = "Filter ROIs")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--raw_steinbock_dir",      type=str,       default="Steinbock")
    parser.add_argument("--dst_steinbock_dir",      type=str,       default="SteinbockNormal")    

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    # load raw slide information
    raw_steinbock_dir = os.path.join(args.data_root, args.data_type, args.raw_steinbock_dir)
    raw_slide_info_path = os.path.join(raw_steinbock_dir, "StudySlide_Info.xlsx")
    raw_slide_df = pd.read_excel(raw_slide_info_path, index_col=None, sheet_name="Sheet1")

    dst_steinbock_dir = os.path.join(args.data_root, args.data_type, args.dst_steinbock_dir)
    dst_roi_info_path = os.path.join(dst_steinbock_dir, "StudyROI_Info.xlsx")
    dst_roi_df = pd.read_excel(dst_roi_info_path, index_col=None, sheet_name="Sheet1")
    dst_roi_list = dst_roi_df["ROI_ID"].tolist()
    normal_slide_set = [ele[:ele.find("-ROI")] for ele in dst_roi_list]
    normal_slide_list = [ele for ele in set(normal_slide_set)]

    dst_slide_df = raw_slide_df[raw_slide_df["Slide_ID"].isin(normal_slide_list)]
    dst_slide_info_path = os.path.join(dst_steinbock_dir, "StudySlide_Info.xlsx")
    dst_slide_df.to_excel(dst_slide_info_path, index=False, header=True) 
