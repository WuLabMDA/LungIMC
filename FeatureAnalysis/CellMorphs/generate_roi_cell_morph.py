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
    parser.add_argument("--feature_dir",            type=str,       default="FeatureAnalysis")    

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    # load cell_cn_types
    phenotype_dir = os.path.join(dataset_dir, args.cellphenotype_dir)  
    cell_type_path = os.path.join(phenotype_dir, "cell_type_cn_morphs.csv")
    celltype_df = pd.read_csv(cell_type_path)
    print("Total cell number is: {}".format(len(celltype_df)))

    # load image list
    roi_info_path = os.path.join(dataset_dir, args.data_type, args.steinbock_dir, "images.csv")
    roi_df = pd.read_csv(roi_info_path)
    img_lst = [os.path.splitext(ele)[0] for ele in roi_df["image"].tolist()]

    # prepare folder
    feature_root_dir = os.path.join(dataset_dir, args.feature_dir)
    celltype_dir = os.path.join(feature_root_dir, "CellMorphs")
    if os.path.exists(celltype_dir):
        shutil.rmtree(celltype_dir)
    os.makedirs(celltype_dir)

    # split
    split_ttl_num = 0
    for img_ind in np.arange(len(img_lst)):
        roi_img_name = img_lst[img_ind]
        roi_df = celltype_df[celltype_df["cell_id"].str.startswith(roi_img_name)]
        split_ttl_num += len(roi_df)
        cell_ids = [ele[len(roi_img_name)+1:] for ele in roi_df["cell_id"].tolist()]
        roi_df.insert(1, "seg_id", cell_ids)
        roi_celltype_path = os.path.join(celltype_dir, roi_img_name + ".csv")
        roi_df.to_csv(roi_celltype_path, index=False)
    print("Split cell num: {}".format(split_ttl_num))