# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import pandas as pd
import numpy as np
from skimage import io
import tifffile
import cv2
import matplotlib.pyplot as plt
import seaborn as sns

from seg_utils import bounding_box


def set_args():
    parser = argparse.ArgumentParser(description = "Check Cell Phenotype on Raw Image")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanSampling35", choices=["HumanWholeIMC", "HumanSampling35"])
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--denoise_dir",            type=str,       default="Denoise")    
    parser.add_argument("--steinbock_dir",          type=str,       default="Steinbock")
    parser.add_argument("--phenotype_dir",          type=str,       default="CellPhenotype")
    parser.add_argument('--divide_ratio',           type=float,     default=10.0)

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = set_args()

    antibody_list = ['B2M', 'B7_H3', 'CD11b', 'CD11c', 'CD14', 'CD163', 'CD19', 
        'CD31', 'CD33', 'CD3e', 'CD4', 'CD45', 'CD45RO', 'CD68', 
        'CD73', 'CD8a', 'CD94', 'CK', 'CTLA_4', 'FoxP3', 'GranzymeB', 
        'HLA_DR', 'ICOS', 'IDO_1', 'Ir191', 'Ki67', 'LAG3', 'MPO', 
        'NaKATPase', 'PD_1', 'PD_L1', 'TIGIT', 'TIM3', 'VISTA', 'aSMA']

    # steinbock dir
    steinbock_dir = os.path.join(args.data_root, args.data_set, args.data_type, args.steinbock_dir)

    # image dir 
    roi_img_dir = os.path.join(args.data_root, args.data_set, args.data_type, args.denoise_dir, "DenoisedROIs")
    # segmentation dir
    cell_seg_dir = os.path.join(steinbock_dir, "masks_deepcell")
    seg_roi_lst = [os.path.splitext(ele)[0] for ele in os.listdir(cell_seg_dir) if ele.endswith(".tiff")]

    # cell phenotype dir
    phenotype_dir = os.path.join(args.data_root, args.data_set, args.data_type, args.phenotype_dir)
    # load cell id & phenotype information
    cell_phenotype_path = os.path.join(phenotype_dir, "ReferenceIDS41ROIs.xlsx")
    cell_phenotype_df = pd.read_excel(cell_phenotype_path)
    cell_ids = cell_phenotype_df["cell ids"].tolist()
    cell_phenotypes = cell_phenotype_df["celltypes"].tolist()

    # accumulate cell ids
    roi_id_dict = {}
    for cur_cell, cur_phenotype in zip(cell_ids, cell_phenotypes):
        if cur_phenotype != "Unknown":
            continue
        roi_name = cur_cell[:cur_cell.find("_")]
        roi_id = int(cur_cell[cur_cell.find("_")+1:])
        if roi_name not in roi_id_dict.keys():
            roi_id_dict[roi_name] = [roi_id, ]
        else:
            roi_id_dict[roi_name].append(roi_id)

    cell_intensity_dict = {ele: [] for ele in antibody_list}
    # collect all information
    for ind, (roi_name, cell_list) in enumerate(roi_id_dict.items()):
        if roi_name not in seg_roi_lst:
            continue      
        
        cur_roi_img_dir = os.path.join(roi_img_dir, roi_name)
        img_list = []
        for antibody in antibody_list:
            cur_antibody_path = os.path.join(cur_roi_img_dir, antibody + ".tiff")
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

        # load cell segmentation
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
            cell_mask = np.zeros((inst_cell_crop.shape[0], inst_cell_crop.shape[1]), dtype=np.uint8)
            cv2.drawContours(cell_mask, contours=contours, contourIdx=0, color=1, thickness=-1)
            for aa, antibody in enumerate(antibody_list):
                cell_intensity_dict[antibody].append(cv2.mean(img_crop[:,:,aa], mask=cell_mask))
    
    # print mean values
    for antibody in cell_intensity_dict.keys():
        print("{}: {}".format(antibody, np.mean(cell_intensity_dict[antibody])))
