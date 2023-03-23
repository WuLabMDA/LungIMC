# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz
import pandas as pd
import numpy as np
from datetime import datetime
from skimage import io
import cv2


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--cellphenotype_dir",      type=str,       default="CellPhenotyping")    
    parser.add_argument("--steinbock_dir",          type=str,       default="Steinbock")
    parser.add_argument("--result_dir",             type=str,       default="Results")
    parser.add_argument("--plot_format",            type=str,       default=".png", choices=[".png", ".pdf"])        

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    # load image width/height information
    steinbock_dir = os.path.join(args.data_root, args.data_set, args.data_type, args.steinbock_dir)
    cellmask_dir = os.path.join(steinbock_dir, "masks_deepcell")
    roi_info_path = os.path.join(steinbock_dir, "images.csv")
    roi_df = pd.read_csv(roi_info_path)
    img_lst = [os.path.splitext(ele)[0] for ele in roi_df["image"].tolist()]
    img_widths = roi_df["width_px"].tolist()
    img_heights = roi_df["height_px"].tolist()

    # load cell type / cn type
    phenotype_dir = os.path.join(args.data_root, args.data_set, args.cellphenotype_dir)  
    cell_type_path = os.path.join(phenotype_dir, "cell_ct_cn_morphs.csv")
    cell_df = pd.read_csv(cell_type_path)
    cell_ids = cell_df["cell_id"].tolist()
    cell_types = cell_df["cell_type"].tolist()
    cell_cns = cell_df["cell_cn"].tolist()

    cell_vis_dir = os.path.join(args.data_root, args.data_set, args.result_dir, "VisCN-Seg")
    if os.path.exists(cell_vis_dir):
        shutil.rmtree(cell_vis_dir)
    os.makedirs(cell_vis_dir)    

    cell_color_dict = {"1": (141, 211, 199), "2": (228, 26, 28), "3": (190, 186, 218), "4": (251, 128, 114), "5": (128, 177, 211),
                       "6": (253, 180, 98), "7": (51, 160, 44), "8": (252, 205, 229), "9": (166, 86, 40), "10": (106, 61, 154)}

    for img_ind in np.arange(len(img_lst)):
        print("CN Seg on {:4d}/{}".format(img_ind+1, len(img_lst)))
        img_name = img_lst[img_ind]
        img_height = img_heights[img_ind]
        img_width = img_widths[img_ind]
        # load cell seg
        roi_cell_seg_path = os.path.join(cellmask_dir, img_name + ".tiff")
        roi_cell_seg = io.imread(roi_cell_seg_path)
        # initizate image
        phenotype_img = np.zeros((img_height, img_width, 3), dtype=np.uint8)
        cell_inds = [ind for ind in np.arange(len(cell_ids)) if cell_ids[ind].startswith(img_name)]
        roi_cell_cns = [cell_cns[ind] for ind in cell_inds]
        roi_cell_ids = [int(cell_ids[ind][len(img_name)+1:]) for ind in cell_inds]
        roi_loc_path = os.path.join(steinbock_dir, "regionprops", img_name + ".csv")
        roi_cen_df = pd.read_csv(roi_loc_path)
        for cell_ind in np.arange(len(cell_inds)):
            cell_cn = roi_cell_cns[cell_ind]
            cell_id = roi_cell_ids[cell_ind]
            roi_cell_mask = (roi_cell_seg == cell_id).astype(np.uint8)
            cell_color = cell_color_dict[str(cell_cn)]
            phenotype_img[:,:,0] += roi_cell_mask * cell_color[0]
            phenotype_img[:,:,1] += roi_cell_mask * cell_color[1]
            phenotype_img[:,:,2] += roi_cell_mask * cell_color[2]
        cell_phenotype_path = os.path.join(cell_vis_dir, img_name + args.plot_format)
        io.imsave(cell_phenotype_path, phenotype_img)







