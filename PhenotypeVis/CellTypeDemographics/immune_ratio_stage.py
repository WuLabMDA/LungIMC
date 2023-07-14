# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns


def set_args():
    parser = argparse.ArgumentParser(description = "Cell type distribution based on sex")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--steinbock_dir",          type=str,       default="Steinbock")
    parser.add_argument("--meta_dir",               type=str,       default="Metadata")
    parser.add_argument("--result_dir",             type=str,       default="Results")
    parser.add_argument("--plot_format",            type=str,       default=".png", choices=[".png", ".pdf"])        

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    # directory setting
    roi_root_dir = os.path.join(args.data_root, args.data_type)
    metadata_dir = os.path.join(args.data_root, args.meta_dir)  
    stat_result_dir = os.path.join(args.data_root, args.result_dir, "Demographics")
    if not os.path.exists(stat_result_dir):
        os.makedirs(stat_result_dir)      

    # load patient information
    roi_info_path = os.path.join(metadata_dir, "ROI_Info.xlsx")
    roi_df = pd.read_excel(roi_info_path)

    # load patient gender information
    roi_lst = roi_df["ROI_ID"].tolist()
    diag_lst = roi_df["ROI_Diag"].tolist()
    roi_diag_dict = {roi:diag for roi, diag in zip(roi_lst, diag_lst)}


    # load cell phenotype information
    cell_phenotype_path = os.path.join(roi_root_dir, args.steinbock_dir, "cell_phenotypes.json")
    cell_phenotype_dict = None
    with open(cell_phenotype_path) as fp:
        cell_phenotype_dict = json.load(fp)
    cell_ids = [ele for ele in cell_phenotype_dict.keys()]
    cell_phenotypes = [ele for ele in cell_phenotype_dict.values()]
    cell_rois = [ele[:ele.rfind("_")] for ele in cell_ids]    

    major_type_dict = {"Epithelial-Cell": "NonImmune", "Endothelial-Cell": "NonImmune", "Fibroblast": "NonImmune", "Other-Immune": "Immune",
                        "NK-Cell": "Immune", "B-Cell": "Immune", "CD4-T-Cell": "Immune", "CD8-T-Cell": "Immune", "T-Reg-Cell": "Immune",
                        "Dendritic-Cell": "Immune", "Neutrophil": "Immune", "Monocyte": "Immune", "Macrophage": "Immune", "MDSC": "Immune"}
    roi_cell_dict = {}
    for roi, phenotype in zip(cell_rois, cell_phenotypes):
        if roi not in roi_cell_dict:
            roi_cell_dict[roi] = [major_type_dict[phenotype]]
        else:
            roi_cell_dict[roi].append(major_type_dict[phenotype])
    immune_ratio_lst = []
    for roi in roi_lst:
        cell_lst = roi_cell_dict[roi]
        immune_ratio_lst.append(np.sum([ele == "Immune" for ele in cell_lst]) * 1.0 / len(cell_lst))

    
    immune_df = pd.DataFrame(list(zip(roi_lst, diag_lst, immune_ratio_lst)), columns=["ROI", "Stage", "ImmuneRatio"])
    sns.swarmplot(data=immune_df, x="Stage", y="ImmuneRatio", order=["Normal", "AAH", "AIS", "MIA", "ADC"], size=2)

    # save plot
    plot_name = "immune_ratio_stage"
    plot_path = os.path.join(stat_result_dir, plot_name + args.plot_format)
    plt.savefig(plot_path, transparent=False, dpi=300)    