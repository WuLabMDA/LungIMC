# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz
import pandas as pd
import numpy as np


def set_args():
    parser = argparse.ArgumentParser(description = "Load data into squidpy image container")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--meta_dir",               type=str,       default="Metadata")
    parser.add_argument("--steinbock_dir",          type=str,       default="Steinbock")
    parser.add_argument("--enrichment_dir",         type=str,       default="Enrichment")
    parser.add_argument("--path_stage",             type=str,       default="Normal", choices=["Normal", "AAH", "AIS", "MIA", "ADC"])      

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = set_args()

    # directory setting
    roi_root_dir = os.path.join(args.data_root, args.data_type)
    metadata_dir = os.path.join(args.data_root, args.meta_dir)    
    enrichment_dir = os.path.join(roi_root_dir, args.enrichment_dir, args.path_stage)
    if os.path.exists(enrichment_dir):
        shutil.rmtree(enrichment_dir)
    os.makedirs(enrichment_dir)

    # load roi stage information
    roi_info_path = os.path.join(metadata_dir, "StudyROI_Info.xlsx")
    study_roi_df = pd.read_excel(roi_info_path)
    study_roi_lst, roi_diag_lst = study_roi_df["ROI_ID"].tolist(), study_roi_df["ROI_Diag"].tolist()
    enrichment_roi_lst = [roi for roi, diag in zip(study_roi_lst, roi_diag_lst) if diag == args.path_stage]

    enrichment_roi_names = [ele + ".tiff" for ele in enrichment_roi_lst]
    steinbock_meta_img_path = os.path.join(roi_root_dir, args.steinbock_dir, "images.csv")
    steinbock_img_df = pd.read_csv(steinbock_meta_img_path)
    stage_img_df = steinbock_img_df[steinbock_img_df["image"].isin(enrichment_roi_names)]
    stage_meta_img_path = os.path.join(enrichment_dir, "images.csv")
    stage_img_df.to_csv(stage_meta_img_path, index=False)
 
    # copy images & segmentations
    enrichment_img_dir = os.path.join(enrichment_dir, "img")
    os.makedirs(enrichment_img_dir)
    enrichment_seg_dir = os.path.join(enrichment_dir, "masks_deepcell")
    os.makedirs(enrichment_seg_dir)
    steinbock_img_dir = os.path.join(roi_root_dir, args.steinbock_dir, "img")
    steinbock_seg_dir = os.path.join(roi_root_dir, args.steinbock_dir, "masks_deepcell")
    for ele in enrichment_roi_lst:
        src_img_path = os.path.join(steinbock_img_dir, ele + ".tiff")
        shutil.copy(src_img_path, enrichment_img_dir)
        src_seg_path = os.path.join(steinbock_seg_dir, ele + ".tiff")
        shutil.copy(src_seg_path, enrichment_seg_dir)




