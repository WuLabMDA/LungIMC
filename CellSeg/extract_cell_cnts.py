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
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--cellseg_dir",            type=str,       default="CellSeg")  
    parser.add_argument("--merge_dir",              type=str,       default="MergeCellSeg")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    roi_dir_root = os.path.join(args.data_root, args.data_type, args.cellseg_dir)
    merge_dir_root = os.path.join(roi_dir_root, args.merge_dir)
    mask_dir_root = os.path.join(merge_dir_root, "Mask")
    cnt_dir_root = os.path.join(merge_dir_root, "Contour")
    if os.path.exists(cnt_dir_root):
        shutil.rmtree(cnt_dir_root)
    os.makedirs(cnt_dir_root)

    # deal with all rois
    roi_list = sorted([os.path.splitext(ele)[0] for ele in os.listdir(mask_dir_root) if ele.endswith(".npy")])
    for ind, cur_roi in enumerate(roi_list):
        print("Work on {}/{}".format(ind+1, len(roi_list)))
        roi_cnt_dict = {}
        mask_path = os.path.join(mask_dir_root, cur_roi + ".npy")
        mask = np.load(mask_path)
        mask_h, mask_w = mask.shape
        cell_num = np.max(mask)
        roi_cnt_dict["height"] = mask_h
        roi_cnt_dict["width"] = mask_w
        roi_cnt_dict["cell_num"] = cell_num
        cell_cnt_dict = {}
        for inst_id in range(1, 1+cell_num):
            # locate the cell mask region
            inst_map = np.array(mask == inst_id, np.uint8)
            contours, hierarchy = cv2.findContours(inst_map, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
            cell_cnt = contours[0]
            cell_cnt_dict[inst_id] = cell_cnt
        roi_cnt_dict["cell_cnt"] = cell_cnt_dict
        # save to json
        roi_cnt_path = os.path.join(cnt_dir_root, cur_roi + ".pkl")
        pickle_out = open(roi_cnt_path, "wb")
        pickle.dump(roi_cnt_dict, pickle_out)
        pickle_out.close()