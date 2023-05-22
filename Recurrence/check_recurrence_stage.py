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
    parser.add_argument("--recur_dir",              type=str,       default="Recurrence")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    lesion_info_path = os.path.join(dataset_dir, "Metadata", "Lesion_Info_Aggregation.csv")
    # load lesions
    lesion_df = pd.read_csv(lesion_info_path)

    recurrence_json_path = os.path.join(dataset_dir, args.recur_dir, "LesionRecurrence.json")
    lesion_recur_dict = None
    with open(recurrence_json_path) as fp:
        lesion_recur_dict = json.load(fp)
    lesion_lst = [ele for ele in lesion_recur_dict.keys()]
    lesion_df = lesion_df[lesion_df["Slide_ID"].isin(lesion_lst)]
    lesion_diags = set(lesion_df["Slide_Diag"].tolist())
    print("Unique diagnosis are: {}".format(lesion_diags))