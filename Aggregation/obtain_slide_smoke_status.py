# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, pickle
from datetime import datetime
import pandas as pd
import numpy as np


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])                        
    parser.add_argument("--aggregation_dir",        type=str,       default="Aggregation")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    patient_info_path = os.path.join(dataset_dir, "NatureFigures", "Fig03", "Patient_Info_PackYear.xlsx")
    lesion_info_path = os.path.join(dataset_dir, "Metadata", "Lesion_Info.xlsx")

    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)
    if not os.path.exists(slide_agg_dir):
        os.makedirs(slide_agg_dir)     

    # split patient into 3 categories ["Never", "Light", "Heavy"]
    patient_df = pd.read_excel(patient_info_path)
    patient_packs = [ele for ele in patient_df["Pack*year"] if ele > 0.01]
    
    smoke_level_lst = []
    patient_pack_thresh = 20.0
    for ele in patient_df["Pack*year"]:
        if ele < 0.01:
            smoke_level_lst.append("Never")
        elif ele < patient_pack_thresh:
            smoke_level_lst.append("Light")
        else:
            smoke_level_lst.append("Heavy")
    smoke_levels = ["Never", "Light", "Heavy"]
    for level in smoke_levels:
        print("{} has {} patients.".format(level, smoke_level_lst.count(level)))
    patient_df["SmokeLevel"] = smoke_level_lst
    patient_smoke_dict = {}
    patient_ids = [ele for ele in patient_df["PatientID"]]
    for pid, smoke in zip(patient_ids, smoke_level_lst):
        patient_smoke_dict[pid] = smoke

    # load lesions
    lesion_df = pd.read_excel(lesion_info_path)
    lesion_df = lesion_df[lesion_df["Slide_Diag"] != "Normal"]
    lesion_df = lesion_df[~(lesion_df.Slide_Diag == "AAH") | ~(lesion_df.Slide_ID == "2571-1D")]
    slide_smoke_lst = []
    lesion_pid_lst = [ele for ele in lesion_df["Patient_ID"]]
    for pid in lesion_pid_lst:
        slide_smoke_lst.append(patient_smoke_dict[pid])
    lesion_df["SmokeLevel"] = slide_smoke_lst
    # save to local
    slide_smoke_path = os.path.join(slide_agg_dir, "lesion_smoke_info.xlsx")
    lesion_df.to_excel(slide_smoke_path, index = False)