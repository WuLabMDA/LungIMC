# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, itertools
import pandas as pd
import numpy as np
from scipy.stats import entropy
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
    # image list
    steinbock_dir = os.path.join(dataset_dir, args.data_type, args.steinbock_dir)
    roi_info_path = os.path.join(steinbock_dir, "images.csv")
    roi_info_df = pd.read_csv(roi_info_path)
    img_lst = [os.path.splitext(ele)[0] for ele in roi_info_df["image"].tolist()]
    
    fea_names = ["ROI_ID", "CN-Richness", "CN-ShannonEntropy", "CN-SimpsonIndex"]
    fea_df = pd.DataFrame(columns=fea_names)    
    cellcn_lst = [ind+1 for ind in np.arange(10)]

    # extract features images by images
    roi_num = len(img_lst)
    for ind in np.arange(roi_num):
        roi_name = img_lst[ind]
        print("Extract {:4d}/{:4d} on roi: {}".format(ind+1, roi_num, roi_name))
        roi_celltype_path = os.path.join(celltype_dir, roi_name + ".csv")
        roi_celltype_df = pd.read_csv(roi_celltype_path)
        roi_cell_num = len(roi_celltype_df)
        cell_cn_lst = roi_celltype_df["cell_cn"].tolist()
        roi_cn_ratios = [(cellcn_lst.count(cell_cn) * 1.0 / roi_cell_num) for cell_cn in cell_cn_lst]
        import pdb; pdb.set_trace()
        # richness
        cn_richness = sum([cn_ratio != 0.0 for cn_ratio in roi_cn_ratios])
        # Shannon index
        shannon_entropy = entropy(roi_cn_ratios)
        # Simpson index
        simpson_index = sum([np.power(cn_ratio, 2) for cn_ratio in roi_cn_ratios])
        # collect features
        fea_lst = [roi_name, cn_richness, shannon_entropy, simpson_index]
        fea_df.loc[len(fea_df)] = fea_lst

    cell_fea_path = os.path.join(feature_root_dir, "CN_DiversityFeas.csv")
    fea_df.to_csv(cell_fea_path, index=False)