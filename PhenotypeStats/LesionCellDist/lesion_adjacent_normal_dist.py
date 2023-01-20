# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt


def set_args():
    parser = argparse.ArgumentParser(description = "Assess cell fraction")
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
    # Filter Tumor ROIs
    study_roi_df = study_roi_df[study_roi_df["ROI_Location"] == "Tumor"]
    roi_id_lst = study_roi_df["ROI_ID"].tolist()
    roi_diag_lst = study_roi_df["ROI_Diag"].tolist()
    # collect stage lesion list
    stage_lst = ["AAH", "AIS", "MIA", "ADC"]
    stage_lesion_dict = {}
    for cur_stage in stage_lst:
        lesion_lst = []
        for roi_id, roi_diag in zip(roi_id_lst, roi_diag_lst):
            if roi_diag == cur_stage:
                lesion_name = roi_id[:roi_id.rfind("-ROI")]
                if lesion_name not in lesion_lst:
                    lesion_lst.append(lesion_name)
        stage_lesion_dict[cur_stage] = lesion_lst
    lesion_merge_lst = stage_lesion_dict["AAH"] + stage_lesion_dict["AIS"] + stage_lesion_dict["MIA"] + stage_lesion_dict["ADC"]

   # load adjacent normal rois
    normal_info_path = os.path.join(metadata_dir, "Normal", "ROIs_AdjacentNormal.xlsx")
    normal_info_df = pd.read_excel(normal_info_path)
    slide_id_lst, roi_id_lst = normal_info_df["Slide_ID"].tolist(), normal_info_df["ROI_ID"].tolist()
    slide_roi_map, roi_slide_map = {}, {}
    for slide_id, roi_id in zip(slide_id_lst, roi_id_lst):
        if slide_id not in slide_roi_map:
            slide_roi_map[slide_id] = [roi_id, ]
        else:
            slide_roi_map[slide_id].append(roi_id)
        if roi_id not in roi_slide_map:
            roi_slide_map[roi_id] = [slide_id, ]
        else:
            roi_slide_map[roi_id].append(slide_id)

    # load cell phenotype information
    cell_phenotype_path = os.path.join(roi_root_dir, "CellPhenotypes.json")
    cell_phenotype_dict = None
    with open(cell_phenotype_path) as fp:
        cell_phenotype_dict = json.load(fp)
    lesion_cell_dict = {}
    cell_roi_lst = []
    for cell_id, cell_type in cell_phenotype_dict.items():
        roi_id = cell_id[:cell_id.rfind("_")]
        if roi_id not in roi_id_lst:
            continue
        cell_roi_lst.append(roi_id)
        roi_lesion_lst = roi_slide_map[roi_id]
        for lesion_id in roi_lesion_lst:
            if lesion_id not in lesion_cell_dict:
                lesion_cell_dict[lesion_id] = [cell_type, ]
            else:
                lesion_cell_dict[lesion_id].append(cell_type)
    
    # load cell phenotype information
    cell_phenotype_path = os.path.join(roi_root_dir, "CellPhenotypes.json")
    cell_phenotype_dict = None
    with open(cell_phenotype_path) as fp:
        cell_phenotype_dict = json.load(fp)
    lesion_cell_dict = {}
    cell_roi_lst = []
    for cell_id, cell_type in cell_phenotype_dict.items():
        roi_id = cell_id[:cell_id.rfind("_")]
        if roi_id not in roi_id_lst:
            continue
        cell_roi_lst.append(roi_id)
        roi_lesion_lst = roi_slide_map[roi_id]
        for lesion_id in roi_lesion_lst:
            if lesion_id not in lesion_cell_dict:
                lesion_cell_dict[lesion_id] = [cell_type, ]
            else:
                lesion_cell_dict[lesion_id].append(cell_type)
    # Organize cell 
    lesion_id_lst = []
    epithelial_lst, endothelial_lst, fibroblast_lst, cd4t_lst, cd8t_lst, treg_lst, bcell_lst = [], [], [], [], [], [], []
    macrophage_lst, monocyte_lst, dendrtic_lst, neutrophil_lst, mdsc_lst, nk_lst, other_immnue_lst = [], [], [], [], [], [], []
    for lesion_id in lesion_merge_lst:
        lesion_id_lst.append(lesion_id)
        if lesion_id not in lesion_cell_dict:
            epithelial_lst.append(0.0)
            endothelial_lst.append(0.0)
            fibroblast_lst.append(0.0)
            cd4t_lst.append(0.0)
            cd8t_lst.append(0.0)
            treg_lst.append(0.0)
            bcell_lst.append(0.0)
            macrophage_lst.append(0.0)
            monocyte_lst.append(0.0)
            dendrtic_lst.append(0.0)
            neutrophil_lst.append(0.0)
            mdsc_lst.append(0.0)
            nk_lst.append(0.0)
            other_immnue_lst.append(0.0)
        else:
            cell_lst = lesion_cell_dict[lesion_id]
            ttl_cell_num = len(cell_lst)
            epithelial_lst.append(len([ele for ele in cell_lst if  ele == "Epithelial-Cell"]) * 1.0 / ttl_cell_num)
            endothelial_lst.append(len([ele for ele in cell_lst if  ele == "Endothelial-Cell"]) * 1.0 / ttl_cell_num)
            fibroblast_lst.append(len([ele for ele in cell_lst if  ele == "Fibroblast"]) * 1.0 / ttl_cell_num)
            cd4t_lst.append(len([ele for ele in cell_lst if  ele == "CD4-T-Cell"]) * 1.0 / ttl_cell_num)
            cd8t_lst.append(len([ele for ele in cell_lst if  ele == "CD8-T-Cell"]) * 1.0 / ttl_cell_num)
            treg_lst.append(len([ele for ele in cell_lst if  ele == "T-Reg-Cell"]) * 1.0 / ttl_cell_num)
            bcell_lst.append(len([ele for ele in cell_lst if  ele == "B-Cell"]) * 1.0 / ttl_cell_num)
            macrophage_lst.append(len([ele for ele in cell_lst if  ele == "Macrophage"]) * 1.0 / ttl_cell_num)
            monocyte_lst.append(len([ele for ele in cell_lst if  ele == "Monocyte"]) * 1.0 / ttl_cell_num)
            dendrtic_lst.append(len([ele for ele in cell_lst if  ele == "Dendritic-Cell"]) * 1.0 / ttl_cell_num)
            neutrophil_lst.append(len([ele for ele in cell_lst if  ele == "Neutrophil"]) * 1.0 / ttl_cell_num)
            mdsc_lst.append(len([ele for ele in cell_lst if  ele == "MDSC"]) * 1.0 / ttl_cell_num)
            nk_lst.append(len([ele for ele in cell_lst if  ele == "NK-Cell"]) * 1.0 / ttl_cell_num)
            other_immnue_lst.append(len([ele for ele in cell_lst if  ele == "Other-Immune"]) * 1.0 / ttl_cell_num)

    df_zip_lst = list(zip(lesion_id_lst, epithelial_lst, endothelial_lst, fibroblast_lst, cd4t_lst, cd8t_lst, treg_lst, bcell_lst,
                          macrophage_lst, monocyte_lst, dendrtic_lst, neutrophil_lst, mdsc_lst, nk_lst, other_immnue_lst))
    df_col_lst = ["Lesion", "Epithelial-Cell", "Endothelial-Cell", "Fibroblast", "CD4-T-Cell", "CD8-T-Cell", "T-Reg-Cell", "B-Cell", 
                  "Macrophage", "Monocyte", "Dendritic-Cell", "Neutrophil", "MDSC", "NK-Cell", "Other-Immune"]
    lesion_ratio_df = pd.DataFrame(df_zip_lst, columns = df_col_lst)
    # lesion_ratio_df = lesion_ratio_df.sort_values(by = ["Epithelial-Cell"], ascending=False)
    # plot a Stacked Bar Chart using matplotlib
    
    lesion_ratio_df.plot(x="Lesion", kind="bar", stacked=True, title="Cell Type Fraction", figsize=(25, 12))
    plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
    plot_name = "AdjacentNormal_lesion_fraction"
    plot_path = os.path.join(stat_result_dir, plot_name + args.plot_format)
    plt.savefig(plot_path, transparent=False, dpi=300)