# -*- coding: utf-8 -*-

import os, sys
import numpy as np
from skimage import io
import shutil, argparse, pytz
import pickle, cv2, warnings


def set_args():
    parser = argparse.ArgumentParser(description = "IMC Cell Feature Extraction")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--cellseg_dir",            type=str,       default="CellSeg")        
    parser.add_argument("--merge_dir",              type=str,       default="MergeCellSeg")
    parser.add_argument("--cnt_dir",                type=str,       default="Contour")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()
    # prepare directory
    root_dir = os.path.join(args.data_root, args.data_type)
    # cell contour directory
    cnt_dir_root = os.path.join(root_dir, args.cellseg_dir, args.merge_dir, args.cnt_dir)
    roi_lst = [os.path.splitext(ele)[0] for ele in os.listdir(cnt_dir_root) if ele.endswith(".pkl")]

    # traverse one-on-one 
    for ind, cur_roi in enumerate(roi_lst):
        print("Cell over in {}/{}".format(ind+1, len(roi_lst)))
        roi_cnt_path = os.path.join(cnt_dir_root, cur_roi + ".pkl")
        roi_cnt_dict = None
        with open(roi_cnt_path, "rb") as file:
            roi_cnt_dict = pickle.load(file)
        mask_h = roi_cnt_dict["height"]
        mask_w = roi_cnt_dict["width"]
        mask_img = np.ones((mask_h, mask_w, 3), dtype=np.uint8) * 255
        cell_cnt_dict = roi_cnt_dict["cell_cnt"]
        cnt_keys = [ele for ele in cell_cnt_dict.keys()]
        cur_cnts = []
        for cur_key in cnt_keys:
            cur_cnt = cell_cnt_dict[cur_key]
            cur_cnts.append(cur_cnt)
        cv2.drawContours(mask_img, cur_cnts, -1, (0, 0, 0), 1)
        ##
        roi_cell_path = os.path.join(cnt_dir_root, cur_roi + ".png")
        io.imsave(roi_cell_path, mask_img)


        
        
