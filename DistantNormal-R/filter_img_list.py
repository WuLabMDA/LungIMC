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

    # load raw slide information
    raw_steinbock_dir = os.path.join(args.data_root, args.data_type, args.raw_steinbock_dir)
    raw_image_path = os.path.join(raw_steinbock_dir, "images.csv")
    raw_image_df = pd.read_csv(raw_image_path)

    dst_steinbock_dir = os.path.join(args.data_root, args.data_type, args.dst_steinbock_dir)
    dst_roi_info_path = os.path.join(dst_steinbock_dir, "StudyROI_Info.xlsx")
    dst_roi_df = pd.read_excel(dst_roi_info_path, index_col=None, sheet_name="Sheet1")
    dst_roi_list = [ele+".tiff" for ele in dst_roi_df["ROI_ID"].tolist()]

    dst_image_df = raw_image_df[raw_image_df["image"].isin(dst_roi_list)]
    dst_image_path = os.path.join(dst_steinbock_dir, "images.csv")
    dst_image_df.to_csv(dst_image_path, index=False, header=True) 