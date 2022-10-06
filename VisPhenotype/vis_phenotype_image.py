# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import pandas as pd
import numpy as np
from skimage import io
import tifffile
import cv2

from seg_utils import bounding_box


def set_args():
    parser = argparse.ArgumentParser(description = "Check Cell Phenotype on Raw Image")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--denoise_dir",            type=str,       default="Denoise")
    parser.add_argument("--steinbock_dir",          type=str,       default="SteinbockAll")
    parser.add_argument("--phenotype_dir",          type=str,       default="PhenotypeAll")
    parser.add_argument("--phenotype_id",           type=int,       default=3)


    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = set_args()

    antibody_list = ["CD8a", "CD3e", "CK"]
    denoised_img_dir = os.path.join(args.data_root, args.data_type, args.denoise_dir, "DenoisedROIs")

    steinbock_dir = os.path.join(args.data_root, args.data_type, args.steinbock_dir)
    # segmentation dir
    cell_seg_dir = os.path.join(steinbock_dir, "masks_deepcell")

    # load cell id & phenotype information
    cell_id_phenotype_path = os.path.join(steinbock_dir, "cell_id_phenotype.csv")
    cell_id_phenotype_df = pd.read_csv(cell_id_phenotype_path, index_col=None)
    cell_ids = cell_id_phenotype_df["CellID"].tolist()
    cell_phenotypes = cell_id_phenotype_df["CellPhenotype"].tolist()
    interested_ids = [cell_ids[ind] for ind in np.arange(len(cell_ids)) if cell_phenotypes[ind] == args.phenotype_id]
    # collect ROIs with ids
    roi_id_dict = {}
    for cur_id in interested_ids:
        cur_roi = cur_id[:cur_id.find("_")]
        cell_id = int(cur_id[cur_id.find("_")+1:])
        if cur_roi in roi_id_dict:
            roi_id_dict[cur_roi].append(cell_id)
        else:
            roi_id_dict[cur_roi] = [cell_id, ]

    cell_phenotype_dir = os.path.join(args.data_root, args.data_type, args.phenotype_dir, "Phenotype"+str(args.phenotype_id))
    if os.path.exists(cell_phenotype_dir):
        shutil.rmtree(cell_phenotype_dir)
    os.makedirs(cell_phenotype_dir)

    # superimpose cells ontop antibodies
    for ind, (roi_name, cell_list) in enumerate(roi_id_dict.items()):
        print("Superimpose on ROI {} {}/{} ".format(roi_name, ind+1, len(roi_id_dict)))
        if ind > 20:
            continue
        img_list = []
        for antibody in antibody_list:
            cur_antibody_path = os.path.join(denoised_img_dir, roi_name, antibody + ".tiff")
            cur_antibody = io.imread(cur_antibody_path, plugin="tifffile").astype(np.float32)
            high_thresh = np.quantile(cur_antibody, 0.96)
            if high_thresh < 0.01:
                cur_antibody = np.zeros_like(cur_antibody)
            else:
                cur_antibody[cur_antibody > high_thresh] = high_thresh
                cur_antibody[cur_antibody < 0.0] = 0.0
                cur_antibody = cur_antibody * 1.0 / high_thresh
            img_list.append(cur_antibody)
        img_arr = (np.dstack(img_list) * 255).astype(np.uint8)
        # overlay cell onto image
        cell_seg_path = os.path.join(cell_seg_dir, roi_name + ".tiff")
        cell_seg = io.imread(cell_seg_path, plugin="tifffile").astype(np.int32)
        for cell_id in cell_list:
            inst_map = np.array(cell_seg==cell_id, np.uint8)
            y1, y2, x1, x2  = bounding_box(inst_map)
            y1 = y1 - 2 if y1 - 2 >= 0 else y1
            x1 = x1 - 2 if x1 - 2 >= 0 else x1
            x2 = x2 + 2 if x2 + 2 <= cell_seg.shape[1] - 1 else x2
            y2 = y2 + 2 if y2 + 2 <= cell_seg.shape[0] - 1 else y2
            inst_cell_crop = inst_map[y1:y2, x1:x2]
            contours, hierarchy = cv2.findContours(inst_cell_crop, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
            contours = sorted(contours, key=lambda x: cv2.contourArea(x), reverse=True)
            img_crop = img_arr[y1:y2, x1:x2]
            cv2.drawContours(img_crop, contours=contours, contourIdx=0, color=(180, 180, 180), thickness=1)
            img_arr[y1:y2, x1:x2] = img_crop
        cell_overlay_path = os.path.join(cell_phenotype_dir, roi_name + ".png")
        io.imsave(cell_overlay_path, img_arr)
        

