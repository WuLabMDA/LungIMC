# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse, pytz
import numpy as np
from skimage import io
from datetime import datetime
import warnings


def set_args():
    parser = argparse.ArgumentParser(description = "IMC Cell Feature Extraction")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--denoise_dir",            type=str,       default="Denoise")
    parser.add_argument("--spillover_dir",          type=str,       default="SpilloverCorrection")

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
    root_dir = os.path.join(args.data_root, args.data_type)
    correct_spillover_dir = os.path.join(root_dir, args.denoise_dir, args.spillover_dir, "Correct")
    roi_spillover_dir = os.path.join(root_dir, args.denoise_dir, "SpilloverROIs")
    roi_list = [os.path.splitext(ele)[0] for ele in os.listdir(correct_spillover_dir) if ele.endswith(".tiff")]
    for ind, cur_roi in enumerate(roi_list):
        if (ind + 1) % 10 == 0:
            print("Split {}/{}".format(ind+1, len(roi_list)))
        cur_roi_path = os.path.join(correct_spillover_dir, cur_roi + ".tiff")
        correct_img = io.imread(cur_roi_path, plugin="tifffile")
        # save images to ROI
        dst_roi_dir = os.path.join(roi_spillover_dir, cur_roi)
        if not os.path.exists(dst_roi_dir):
            os.makedirs(dst_roi_dir)        
        for ii, antibody in enumerate(antibody_list):
            cur_antibody = correct_img[ii, :, :]
            cur_antibody_path = os.path.join(dst_roi_dir, antibody + ".tiff")
            io.imsave(cur_antibody_path, cur_antibody, plugin="tifffile")