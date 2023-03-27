# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, itertools
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
    parser.add_argument("--steinbock_dir",          type=str,       default="Steinbock")
    parser.add_argument("--feature_dir",            type=str,       default="FeatureAnalysis") 

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set) 
    feature_root_dir = os.path.join(dataset_dir, args.feature_dir)
    celltype_dir = os.path.join(dataset_dir, "CellPhenotyping", "Cell_CT_CN_Morphs")   

    steinbock_dir = os.path.join(dataset_dir, args.data_type, args.steinbock_dir)
    roi_info_path = os.path.join(steinbock_dir, "images.csv")
    roi_info_df = pd.read_csv(roi_info_path)
    img_lst = [os.path.splitext(ele)[0] for ele in roi_info_df["image"].tolist()]
    img_heights = [int(ele) for ele in roi_info_df["acquisition_height_um"].tolist()]
    img_widths = [int(ele) for ele in roi_info_df["acquisition_width_um"].tolist()]

    celltype_lst = ["Epithelial-Cell", "Endothelial-Cell", "Fibroblast", "CD4-T-Cell", "CD8-T-Cell", 
                    "T-Reg-Cell", "B-Cell", "Macrophage", "Monocyte", "Dendritic-Cell", 
                    "Neutrophil", "MDSC", "NK-Cell", "Proliferating-Cell", "Undefined"]
    cellcn_lst = [str(ind+1) for ind in np.arange(10)]

    # Construct ct-cn combintation list
    fea_names = ["ROI_ID", ]    
    ct_cn_combs = []
    fea_ind = 0
    fea_dict = {}
    for cell_type in celltype_lst:
        for cn_type in cellcn_lst:
            ct_cn_comb = cell_type + "-" + cn_type
            ct_cn_combs.append(ct_cn_comb)
            fea_dict[ct_cn_comb] = fea_ind
            fea_ind += 1
    fea_names.extend(ct_cn_combs)
    fea_df = pd.DataFrame(columns=fea_names)
 
    # Extract features images by images
    roi_num = len(img_lst)
    for ind in np.arange(roi_num):
        roi_name = img_lst[ind]
        print("Extract {:4d}/{:4d} on roi: {}".format(ind+1, roi_num, roi_name))
        roi_celltype_path = os.path.join(celltype_dir, roi_name + ".csv")
        roi_celltype_df = pd.read_csv(roi_celltype_path)
        roi_cell_num = len(roi_celltype_df)
        cell_type_lst = roi_celltype_df["cell_type"].tolist()
        cell_cn_lst = roi_celltype_df["cell_cn"].tolist()   
        roi_area = img_heights[ind] * img_widths[ind]

        fea_lst = [roi_name, ]
        ct_cn_feas = [0.0, ] * len(ct_cn_combs)
        for ct, cn in zip(cell_type_lst, cell_cn_lst):
            ct_cn_feas[fea_dict[ct + "-" + str(cn)]] += 1
        norm_ct_cn_feas = [ele * 1.0 / roi_area for ele in ct_cn_feas]
        fea_lst.extend(norm_ct_cn_feas)
        fea_df.loc[len(fea_df)] = fea_lst

    # Save the features
    cell_fea_path = os.path.join(feature_root_dir, "JointDistCTCNFeas.csv")
    fea_df.to_csv(cell_fea_path, index=False)