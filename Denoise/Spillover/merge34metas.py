# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse, pytz
import numpy as np
from skimage import io
from scipy.io import loadmat
from datetime import datetime
import warnings
import xtiff


def set_args():
    parser = argparse.ArgumentParser(description = "IMC Cell Feature Extraction")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="Study", choices = ["Study", "Tonsil"])
    parser.add_argument("--merge_dir",              type=str,       default="MergeSeg")
    parser.add_argument("--stain_dir",              type=str,       default="Stains")
    parser.add_argument("--spillover_dir",          type=str,       default="Spillover")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args() 

    # antibody list for spillover correction
    antibody_list = ['CD45RO', 'aSMA', 'ICOS', 'HLA_DR', 'CD68', 'MPO', 'TIGIT', 'CD11c', 'CD73', 'PD_L1', 
                    'CD163', 'CD45', 'CD11b', 'CD14', 'FoxP3', 'TIM3', 'LAG3', 'CD31', 'IDO_1', 'Ki67', 
                    'VISTA', 'B2M', 'PD_1', 'CD8a', 'CD33', 'B7_H3', 'GranzymeB', 'CD94', 'CD19', 'CD3e', 
                    'CD4', 'CK', 'CTLA_4', 'NaKATPase']    

    # prepare directory
    root_dir = os.path.join(args.data_root, args.data_type + "Processing")
    raw_spillover_dir = os.path.join(root_dir, args.merge_dir, args.spillover_dir, "Raw")
    if os.path.exists(raw_spillover_dir):
        shutil.rmtree(raw_spillover_dir)
    os.makedirs(raw_spillover_dir)

    stain_root = os.path.join(root_dir, args.merge_dir, args.stain_dir)
    roi_list = [ele for ele in os.listdir(stain_root) if os.path.isdir(os.path.join(stain_root, ele))]
    for ind, cur_roi in enumerate(roi_list):
        if (ind + 1) % 100 == 0:
            print("Merge {}/{}".format(ind+1, len(roi_list)))
        cur_roi_dir = os.path.join(stain_root, cur_roi)
        roi_data = []
        roi_data.append([io.imread(os.path.join(cur_roi_dir, ele + ".tiff")) for ele in antibody_list])
        roi_data = np.stack(roi_data, axis=0)
        roi_data = roi_data.astype(np.uint16)
        tiff_path = os.path.join(raw_spillover_dir, cur_roi + ".tiff")
        xtiff.to_tiff(roi_data, tiff_path)
