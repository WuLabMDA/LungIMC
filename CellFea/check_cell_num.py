# -*- coding: utf-8 -*-

import os, sys
import argparse
import numpy as np


def set_args():
    parser = argparse.ArgumentParser(description = "Filter Segmented Cells")
    parser.add_argument("--data_root",        type=str,       default="/Data")
    parser.add_argument("--data_type",        type=str,       default="LungROIProcessing")
    parser.add_argument("--cellinfo_dir",     type=str,       default="CellInfo")        
    parser.add_argument("--fea_dir",          type=str,       default="CellFeas")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()
    cell_fea_dir = os.path.join(args.data_root, args.data_type, args.cellinfo_dir, args.fea_dir)
    roi_list = [ele for ele in os.listdir(cell_fea_dir) if ele.endswith("npy")]
    ttl_roi_num = len(roi_list)
    ttl_cell_num = 0
    for cur_roi in roi_list:
        cell_fea_path = os.path.join(cell_fea_dir, cur_roi)
        cell_feas = np.load(cell_fea_path)
        ttl_cell_num += cell_feas.shape[0]
    print("There are {} ROIs with {} cells.".format(ttl_roi_num, ttl_cell_num))
    avg_cell_roi = ttl_cell_num * 1.0 / ttl_roi_num
    print("On average, there are {:.1f} cells per ROI.".format(avg_cell_roi))
