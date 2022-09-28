# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import pandas as pd


def set_args():
    parser = argparse.ArgumentParser(description = "Filter ROIs")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--raw_steinbock_dir",      type=str,       default="SteinbockAll")
    parser.add_argument("--dst_steinbock_dir",      type=str,       default="SteinbockDistantNormal")    

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    raw_steinbock_dir = os.path.join(args.data_root, args.data_type, args.raw_steinbock_dir)
    raw_roi_info_path = os.path.join(raw_steinbock_dir, "StudyROI_Info.xlsx")
    raw_roi_df = pd.read_excel(raw_roi_info_path, index_col=None, sheet_name="Sheet1")
    # filter Normal & DistantNormal ROIs
    normal_roi_df = raw_roi_df[raw_roi_df["ROI_Location"] == "DistantNormal"]
    normal_japan_roi_df = normal_roi_df[normal_roi_df["ROI_ID"].str.startswith("H")]
    dst_steinbock_dir = os.path.join(args.data_root, args.data_type, args.dst_steinbock_dir)
    dst_roi_info_path = os.path.join(dst_steinbock_dir, "StudyROI_Info.xlsx")
    normal_japan_roi_df.to_excel(dst_roi_info_path, index=False, header=True)