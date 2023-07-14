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
    parser = argparse.ArgumentParser(description = "Merge stains for Cell Segmentation")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--denoise_dir",            type=str,       default="Denoise")
    parser.add_argument("--cellseg_dir",            type=str,       default="CellSeg")
    parser.add_argument("--segstain_dir",           type=str,       default="SegStainROIs")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()  

    # prepare directory
    root_dir = os.path.join(args.data_root, args.data_type)
    raw_roi_dir = os.path.join(root_dir, args.denoise_dir, "DenoisedROIs")
    segstain_roi_dir = os.path.join(root_dir, args.cellseg_dir, args.segstain_dir)
    if os.path.exists(segstain_roi_dir):
        shutil.rmtree(segstain_roi_dir)
    os.makedirs(segstain_roi_dir)

    roi_list = [ele for ele in os.listdir(raw_roi_dir) if os.path.isdir(os.path.join(raw_roi_dir, ele))]
    for ind, cur_roi in enumerate(roi_list):
        if (ind + 1) % 100 == 0:
            print("Merge {}/{}".format(ind+1, len(roi_list)))
        cur_roi_dir = os.path.join(raw_roi_dir, cur_roi)
        dst_roi_dir = os.path.join(segstain_roi_dir, cur_roi, cur_roi)
        os.makedirs(dst_roi_dir)
        # save membrane image
        mem_tiff_path = os.path.join(cur_roi_dir, "NaKATPase.tiff")
        mem_img = io.imread(mem_tiff_path, plugin="tifffile")
        mem_max = np.percentile(mem_img, 99)
        mem_min = np.percentile(mem_img, 1)
        if mem_max - mem_min < 1:
            mem_img = np.zeros_like(mem_img)
        else:
            mem_img[mem_img > mem_max] = mem_max
            mem_img[mem_img < mem_min] = mem_min
            mem_img = (1.0 * (mem_img - mem_min) / (mem_max - mem_min) * 255.0)
        mem_img = mem_img.astype(np.uint8)
        mem_tif_path = os.path.join(dst_roi_dir, "NaKATPase.tif")
        io.imsave(mem_tif_path, mem_img)
        # save nuclear image
        nuc_tiff_path = os.path.join(cur_roi_dir, "Ir191.tiff")
        nuc_img = io.imread(nuc_tiff_path, plugin="tifffile")
        nuc_max = np.percentile(nuc_img, 99)
        nuc_min = np.percentile(nuc_img, 1)
        if nuc_max - nuc_min < 1:
            nuc_img = np.zeros_like(nuc_img)
        else:
            nuc_img[nuc_img > nuc_max] = nuc_max
            nuc_img[nuc_img < nuc_min] = nuc_min
            nuc_img = (1.0 * (nuc_img - nuc_min) / (nuc_max - nuc_min) * 255.0)
        nuc_img = nuc_img.astype(np.uint8)
        nuc_tif_path = os.path.join(dst_roi_dir, "Ir191.tif")
        io.imsave(nuc_tif_path, nuc_img)