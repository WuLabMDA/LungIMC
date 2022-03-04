# -*- coding: utf-8 -*-

import os, sys
import argparse, shutil
import numpy as np


def set_args():
    parser = argparse.ArgumentParser(description = "IMC Cell Segmentation")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="Study", choices = ["Study", "Tonsil"])
    parser.add_argument("--roi_dir",                type=str,       default="SegROI")
    parser.add_argument("--result_dir",             type=str,       default="SegResults")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()
    seg_root_dir = os.path.join(args.data_root, args.data_type + "Processing")
    segroi_dir = os.path.join(seg_root_dir, args.roi_dir)
    segresult_dir = os.path.join(seg_root_dir, args.result_dir)
    if os.path.exists(segresult_dir):
        shutil.rmtree(segresult_dir)
    os.makedirs(segresult_dir)

    # stains
    nucs = ["191Ir", ]
    mems = ["NaKATPase", ]
    cell_stains = [nucs[0], mems[0]]
    # traverse patients one-by-one
    patient_list = [ele for ele in os.listdir(segroi_dir)]
    for ind, p_id in enumerate(patient_list):
        print("Segment {}/{}".format(ind+1, len(patient_list)))
