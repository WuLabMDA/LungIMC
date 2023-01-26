# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt


def set_args():
    parser = argparse.ArgumentParser(description = "Assess the cell numbers of 3 major types")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
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
    stat_result_dir = os.path.join(args.data_root, args.result_dir, "Stats")
    if not os.path.exists(stat_result_dir):
        os.makedirs(stat_result_dir)
    
    # load roi stage information
    roi_info_path = os.path.join(metadata_dir, "StudyROI_Info.xlsx")
    study_roi_df = pd.read_excel(roi_info_path)
    study_roi_lst = study_roi_df["ROI_ID"].tolist()
    roi_diag_lst = study_roi_df["ROI_Diag"].tolist()
    roi_stage_dict = {roi:diag for roi, diag in zip(study_roi_lst, roi_diag_lst)}

    # load cell phenotype information
    cell_phenotype_path = os.path.join(roi_root_dir, "CellPhenotypes.json")
    cell_phenotype_dict = None
    with open(cell_phenotype_path) as fp:
        cell_phenotype_dict = json.load(fp)
    cell_ids = [ele for ele in cell_phenotype_dict.keys()]
    cell_phenotypes = [ele for ele in cell_phenotype_dict.values()]
    cell_rois = [ele[:ele.rfind("_")] for ele in cell_ids]
    
    # obtain all cell stages
    cell_stages = [roi_stage_dict[ele] for ele in cell_rois]
    # merge major cell types
    major_type_dict = {"Epithelial-Cell": "Epithelial", "Endothelial-Cell": "Stromal", "Fibroblast": "Stromal", "Other-Immune": "Immune",
                       "NK-Cell": "Immune", "B-Cell": "Immune", "CD4-T-Cell": "Immune", "CD8-T-Cell": "Immune", "T-Reg-Cell": "Immune",
                       "Dendritic-Cell": "Immune", "Neutrophil": "Immune", "Monocyte": "Immune", "Macrophage": "Immune", "MDSC": "Immune"}
    cell_types = [major_type_dict[ele] for ele in cell_phenotypes]

    # collect information
    stage_lst = ["AAH", "AIS", "MIA", "ADC"]
    ratio_dict = {}
    stromal_lst, immune_lst, epithelial_lst = [], [], []
    for cur_stage in stage_lst:
        cur_cell_types = [cell_type for cell_type, cell_stage in zip(cell_types, cell_stages) if cell_stage == cur_stage]
        # print("{} has {} cells".format(cur_stage, len(cur_cell_types)))
        stromal_num = len([ele for ele in cur_cell_types if  ele == "Stromal"]) 
        stromal_lst.append(stromal_num)
        immune_num = len([ele for ele in cur_cell_types if  ele == "Immune"]) 
        immune_lst.append(immune_num)
        epithelial_num = len([ele for ele in cur_cell_types if ele == "Epithelial"])
        epithelial_lst.append(epithelial_num)
    type_stage_df = pd.DataFrame(list(zip(stage_lst, epithelial_lst, stromal_lst, immune_lst)), 
                                          columns =["Stage", "Epithelial", "Stromal", "Immune"])
    
    # plot a Stacked Bar Chart using matplotlib
    type_stage_df.plot(x = "Stage", kind = "bar", stacked = True, title = "Major Cell Type Number")
    plot_name = "major_3type_stage_num"
    plot_path = os.path.join(stat_result_dir, plot_name + args.plot_format)
    plt.savefig(plot_path, transparent=False, dpi=300)     