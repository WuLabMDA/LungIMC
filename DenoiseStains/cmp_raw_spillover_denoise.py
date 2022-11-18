# -*- coding: utf-8 -*-

import os, sys
import numpy as np
import argparse, pytz, shutil
from skimage import io
from datetime import datetime


def set_args():
    parser = argparse.ArgumentParser(description = "Comparison of Raw/SpilloverCorrection/Denoising")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--denoise_dir",            type=str,       default="Denoise")
    parser.add_argument('--divide_ratio',           type=float,     default=20.0)

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = set_args()

    denoise_root_dir = os.path.join(args.data_root, args.data_type, args.denoise_dir)
    raw_dir = os.path.join(denoise_root_dir, "RawROIs")
    spillover_dir = os.path.join(denoise_root_dir, "SpilloverROIs")
    denoise_dir = os.path.join(denoise_root_dir, "DenoisedROIs")
    cmp_dir = os.path.join(denoise_root_dir, "CmpRawSpilloverDenoise")

    roi_lst = sorted([ele for ele in os.listdir(raw_dir) if os.path.isdir(os.path.join(raw_dir, ele))])
    for ind, cur_roi in enumerate(roi_lst):
        if (ind + 1) % 10 == 0:
            print("Processing {}/{}".format(ind+1, len(roi_lst)))
            
        raw_roi_dir = os.path.join(raw_dir, cur_roi)
        spillover_roi_dir = os.path.join(spillover_dir, cur_roi)
        denoise_roi_dir = os.path.join(denoise_dir, cur_roi)
        cmp_roi_dir = os.path.join(cmp_dir, cur_roi)
        if not os.path.exists(cmp_roi_dir):
            os.makedirs(cmp_roi_dir)

        img_lst = [os.path.splitext(ele)[0] for ele in os.listdir(raw_roi_dir) if ele.endswith(".tiff")]
        for cur_img in img_lst:
            raw_img_path = os.path.join(raw_roi_dir, cur_img + ".tiff")
            spillover_img_path = os.path.join(spillover_roi_dir, cur_img + ".tiff")
            denoise_img_path = os.path.join(denoise_roi_dir, cur_img + ".tiff")
            # transform raw
            raw_img = io.imread(raw_img_path, plugin='tifffile').astype(np.float32)
            high_thresh = np.max(raw_img) * 1.0 / args.divide_ratio
            raw_img[raw_img > high_thresh] = high_thresh
            raw_img = (255.0 * (raw_img / high_thresh)).astype(np.uint8)
            # transform spillover
            spillover_img = io.imread(spillover_img_path, plugin='tifffile').astype(np.float32)
            spillover_img[spillover_img > high_thresh] = high_thresh
            spillover_img = (255.0 * (spillover_img / high_thresh)).astype(np.uint8)
            # transform denoise image
            denoise_img = io.imread(denoise_img_path, plugin='tifffile').astype(np.float32)
            denoise_img[denoise_img > high_thresh] = high_thresh
            denoise_img = (255.0 * (denoise_img / high_thresh)).astype(np.uint8)
            # combine three images together
            cmp_img = np.zeros((raw_img.shape[0], raw_img.shape[1] * 3), dtype=np.uint8)
            cmp_img[:,:raw_img.shape[1]] = raw_img
            cmp_img[:,raw_img.shape[1]:2*raw_img.shape[1]] = spillover_img
            cmp_img[:,2*raw_img.shape[1]:] = denoise_img
            # save combined image
            cmp_img_path = os.path.join(cmp_roi_dir, cur_img + ".png")
            io.imsave(cmp_img_path, cmp_img)


