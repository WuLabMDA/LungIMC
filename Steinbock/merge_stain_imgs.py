# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse, pytz
import pandas as pd
import numpy as np
from skimage import io
from scipy.io import loadmat
from datetime import datetime
import warnings
import tifffile


def set_args():
    parser = argparse.ArgumentParser(description = "Merge stains for spillover correction")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--denoise_dir",            type=str,       default="Denoise")
    parser.add_argument("--steinbock_dir",          type=str,       default="Steinbock")    

    args = parser.parse_args()
    return args



if __name__ == "__main__":
    args = set_args()

    steinbock_dir = os.path.join(args.data_root, args.data_type, args.steinbock_dir)
    if not os.path.exists(steinbock_dir):
        os.makedirs(steinbock_dir)
    stain_tif_dir = os.path.join(args.data_root, args.data_type, args.steinbock_dir, "img")
    if os.path.exists(stain_tif_dir):
        shutil.rmtree(stain_tif_dir)
    os.makedirs(stain_tif_dir)

    # load stain panel
    stain_panel_path = os.path.join(args.data_root, args.data_type, args.steinbock_dir, "panel.csv")
    stain_panel = pd.read_csv(stain_panel_path)
    antibody_list = stain_panel["name"].tolist()

    image_lst, width_lst, height_lst, channel_lst = [], [], [], []
    stain_roi_dir = os.path.join(args.data_root, args.data_type, args.denoise_dir, "DenoisedROIs")
    roi_list = [ele for ele in os.listdir(stain_roi_dir) if os.path.isdir(os.path.join(stain_roi_dir, ele))]
    for ind, cur_roi in enumerate(roi_list):
        if (ind + 1) % 10 == 0:
            print("Merge {}/{}".format(ind+1, len(roi_list)))
        cur_roi_dir = os.path.join(stain_roi_dir, cur_roi)
        roi_data = []
        roi_data.append([io.imread(os.path.join(cur_roi_dir, ele + ".tiff")) for ele in antibody_list])
        roi_data = np.stack(roi_data, axis=0).astype(np.float32)
        # save image
        tif_img_name = cur_roi + ".tiff"
        tif_img_path = os.path.join(stain_tif_dir, tif_img_name)
        tifffile.imwrite(tif_img_path,
            data=roi_data[np.newaxis, np.newaxis, :, :, :, np.newaxis],
            imagej=roi_data.dtype in (np.uint8, np.uint16, np.float32))
        # add meta info
        image_lst.append(tif_img_name)
        width_lst.append(roi_data.shape[2])
        height_lst.append(roi_data.shape[1])
        channel_lst.append(roi_data.shape[0])
    
    # save meta info
    stain_images_path = os.path.join(args.data_root, args.data_type, args.steinbock_dir, "images.csv")
    stain_iamges_df = pd.DataFrame(list(zip(image_lst, width_lst, height_lst, channel_lst, width_lst, height_lst)),
        columns =["image", "width_px", "height_px", "num_channels", "acquisition_width_um", "acquisition_height_um"])
    stain_iamges_df.to_csv(stain_images_path, index=False)