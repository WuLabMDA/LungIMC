# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse, json, pickle
import deepdish as dd
import numpy as np
from skimage import io
from scipy.io import loadmat
import tifffile, cv2
import warnings
warnings.simplefilter("ignore", UserWarning)


def set_args():
    parser = argparse.ArgumentParser(description = "Cell Segmentation Control Extraction")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="Study", choices = ["Study", "Tonsil"])
    parser.add_argument("--merge_dir",              type=str,       default="MergeSeg")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    roi_dir_root = os.path.join(args.data_root, args.data_type + "Processing")
    merge_dir_root = os.path.join(roi_dir_root, args.merge_dir)
    mask_dir_root = os.path.join(merge_dir_root, "Mask")
    tiff_dir_root = os.path.join(merge_dir_root, "SegTIF16")
    if os.path.exists(tiff_dir_root):
        shutil.rmtree(tiff_dir_root)
    os.makedirs(tiff_dir_root)

    # deal with all rois
    roi_list = sorted([os.path.splitext(ele)[0] for ele in os.listdir(mask_dir_root) if ele.endswith(".npy")])
    # cell_num_list = []
    for ind, cur_roi in enumerate(roi_list):
        if (ind + 1) % 100 == 0:
            print("Merge {}/{}".format(ind+1, len(roi_list)))
        roi_cnt_dict = {}
        mask_path = os.path.join(mask_dir_root, cur_roi + ".npy")
        mask = np.load(mask_path).astype(np.uint16)
        seg_tif_path = os.path.join(tiff_dir_root, cur_roi + ".tiff") 
        io.imsave(seg_tif_path, mask, plugin="tifffile")