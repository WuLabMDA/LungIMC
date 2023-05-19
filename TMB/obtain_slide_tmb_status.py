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
    parser.add_argument("--tmb_dir",                type=str,       default="TMB")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    lesion_info_path = os.path.join(dataset_dir, "Metadata", "Lesion_Info_Aggregation.csv")

    slide_tmb_dir = os.path.join(dataset_dir, args.tmb_dir)
    if not os.path.exists(slide_tmb_dir):
        os.makedirs(slide_tmb_dir)     

    # load lesions
    lesion_df = pd.read_csv(lesion_info_path)
    import pdb; pdb.set_trace()

    # # lesion_df = lesion_df[lesion_df["Slide_Diag"] != "Normal"]
    # lesion_df = lesion_df[~(lesion_df.Slide_Diag == "AAH") | ~(lesion_df.Slide_ID == "2571-1D")]
    # slide_smoke_lst = []
    # lesion_pid_lst = [ele for ele in lesion_df["Patient_ID"]]
    # for pid in lesion_pid_lst:
    #     slide_smoke_lst.append(patient_smoke_dict[pid])
    # lesion_df["SmokeLevel"] = slide_smoke_lst


    # # save to local
    # slide_tmb_path = os.path.join(slide_tmb_dir, "lesion_tmb_info.xlsx")
    # tmb_df.to_excel(slide_tmb_path, index = False)