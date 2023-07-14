# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz
from datetime import datetime
import pandas as pd
import numpy as np
from scipy import stats


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
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)
    slide_smoke_path = os.path.join(slide_agg_dir, "lesion_smoke_info.xlsx")
    lesion_df = pd.read_excel(slide_smoke_path)
    # filtering TMB/NA
    lesion_df = lesion_df[~pd.isnull(lesion_df["TMB"])]
    # lesion_df.loc[lesion_df["SmokeLevel"] == "Light", "SmokeLevel"] = "Never" # combine Light to Never
    lesion_df = lesion_df[lesion_df["SmokeLevel"].isin(["Never", "Heavy"])]
    print("Within {} lesion: ".format(lesion_df.shape[0]))
    tmb_map_dict = {"Low": 0, "High": 1}
    smoke_map_dict = {"Never": 0, "Heavy": 1}
    tmb_lst = [tmb_map_dict[ele] for ele in lesion_df["TMB"]]
    smoke_lst = [smoke_map_dict[ele] for ele in lesion_df["SmokeLevel"]]
    p_corr, _ = stats.pearsonr(tmb_lst, smoke_lst)
    print("Pearson correlation: {:.3f}".format(p_corr))
    s_corr, _ = stats.spearmanr(tmb_lst, smoke_lst)
    print("Spearman correlation: {:.3f}".format(s_corr))

    