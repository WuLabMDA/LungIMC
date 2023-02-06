# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import scanpy as sc
import anndata
import squidpy as sq


def set_args():
    parser = argparse.ArgumentParser(description = "Load data into squidpy image container")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--meta_dir",               type=str,       default="Metadata")
    parser.add_argument("--steinbock_dir",          type=str,       default="Steinbock")
    parser.add_argument("--result_dir",             type=str,       default="Results")
    parser.add_argument("--path_stage",             type=str,       default="Normal", choices=["Normal", "AAH", "AIS", "MIA", "ADC"])      

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    # directory setting
    steinbock_root_dir = os.path.join(args.data_root, args.data_type)
    metadata_dir = os.path.join(args.data_root, args.meta_dir)
    enrich_result_dir = os.path.join(args.data_root, args.result_dir, "Enrichment")
    if not os.path.exists(enrich_result_dir):
        os.makedirs(enrich_result_dir)

    # load roi stage information
    roi_info_path = os.path.join(metadata_dir, "StudyROI_Info.xlsx")
    study_roi_df = pd.read_excel(roi_info_path)
    study_roi_lst, roi_diag_lst = study_roi_df["ROI_ID"].tolist(), study_roi_df["ROI_Diag"].tolist()
    enrichment_roi_lst = [roi for roi, diag in zip(study_roi_lst, roi_diag_lst) if diag == args.path_stage]
    print("{} has {} ROIs.".format(args.path_stage, len(enrichment_roi_lst)))
    import pdb; pdb.set_trace()
