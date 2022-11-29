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
    parser.add_argument("--steinbock_dir",          type=str,       default="Steinbock")
    parser.add_argument("--phenotype_dir",          type=str,       default="CellPhenotype")
    parser.add_argument('--divide_ratio',           type=float,     default=20.0)

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = set_args()

    antibody_list = ["CD8a", "FoxP3", "CK"]

    cell_phenotypes = ["B-cell", "CD3T-cell", "CD4T-cell", "CD8T-cell", "Dendritic", 
        "Endothelial", "Epithelial", "Macrophage", "MDSC", "Monocytes", 
        "Neutrophils", "NK", "SmoothMuscle-Stromal", "T-reg", "Unknown"]
    phenotype_colors = [(128, 0, 0), (170, 110, 40), (128, 128, 0), (0, 128, 128), (0, 0, 128),
        (245, 130, 48), (255, 225, 25), (210, 245, 60), (70, 240, 240), (145, 30, 180),
        (250, 190, 212), (255, 215, 180), (255, 250, 200), (170, 255, 195), (220, 190, 255)]

    # steinbock dir
    steinbock_dir = os.path.join(args.data_root, args.data_type, args.steinbock_dir)
    # image dir 
    roi_img_dir = os.path.join(args.data_root, args.data_type, args.denoise_dir, "DenoisedROIs")
    # segmentation dir
    cell_seg_dir = os.path.join(steinbock_dir, "masks_deepcell")

    # cell phenotype dir
    phenotype_dir = os.path.join(args.data_root, args.data_type, args.phenotype_dir)
    # load cell id & phenotype information
    cell_phenotype_path = os.path.join(phenotype_dir, "ReferenceIDS.xlsx")
    cell_phenotype_df = pd.read_excel(cell_phenotype_path)
    cell_ids = cell_phenotype_df["IDS"].tolist()
    cell_types = cell_phenotype_df["celltypes"].tolist()
    cell_colors = []
    for cur_phenotype in cell_types:
        for ind, cell_type in enumerate(cell_phenotypes):
            if cur_phenotype.startswith(cell_type):
                cell_colors.append(phenotype_colors[ind])
    # print("Number of cells: {}".format(len(cell_ids)))
    # print("Number of colors: {}".format(len(cell_colors)))

    # accumulate cell ids
    roi_id_dict, roi_color_dict = {}, {}
    for cur_cell, cur_color in zip(cell_ids, cell_colors):
        roi_name = cur_cell[:cur_cell.find("_")]
        roi_id = int(cur_cell[cur_cell.find("_")+1:])
        if roi_name not in roi_id_dict.keys():
            roi_id_dict[roi_name] = [roi_id, ]
            roi_color_dict[roi_name] = [cur_color, ]
        else:
            roi_id_dict[roi_name].append(roi_id)
            roi_color_dict[roi_name].append(cur_color)

    cell_phenotype_dir = os.path.join(args.data_root, args.data_type, args.phenotype_dir, "EntireCD8a")
    if os.path.exists(cell_phenotype_dir):
        shutil.rmtree(cell_phenotype_dir)
    os.makedirs(cell_phenotype_dir)    

    # superimpose cells ontop antibodies
    for ind, (roi_name, cell_list) in enumerate(roi_id_dict.items()):
        print("Superimpose on ROI {} {}/{} ".format(roi_name, ind+1, len(roi_id_dict)))
        img_list = []
        for antibody in antibody_list:
            cur_antibody_path = os.path.join(roi_img_dir, roi_name, antibody + ".tiff")
            cur_antibody = io.imread(cur_antibody_path, plugin="tifffile").astype(np.float32)
            high_thresh = np.max(cur_antibody) * 1.0 / args.divide_ratio     
            if high_thresh < 0.01:
                cur_antibody = np.zeros_like(cur_antibody)
            else:
                cur_antibody[cur_antibody > high_thresh] = high_thresh
                cur_antibody[cur_antibody < 0.0] = 0.0
                cur_antibody = cur_antibody * 1.0 / high_thresh
            img_list.append(cur_antibody)
        img_arr = (np.dstack(img_list) * 255).astype(np.uint8)

        color_list = roi_color_dict[roi_name]
        assert len(cell_list) == len(color_list), "On ROI {}, cell and color number not matched"
        # overlay cell onto image
        cell_seg_path = os.path.join(cell_seg_dir, roi_name + ".tiff")
        cell_seg = io.imread(cell_seg_path, plugin="tifffile").astype(np.int32)
        for cind, cell_id in enumerate(cell_list):
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
            cv2.drawContours(img_crop, contours=contours, contourIdx=0, color=color_list[cind], thickness=1)
            img_arr[y1:y2, x1:x2] = img_crop
        cell_overlay_path = os.path.join(cell_phenotype_dir, roi_name + ".png")
        io.imsave(cell_overlay_path, img_arr)
        

