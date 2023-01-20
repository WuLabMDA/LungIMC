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

    lesion_merge_lst = ['2493-6M', '2513-1G', 'H18-0130-6', '2571-4B', 'H09-0054-3', 'H11-0383-6', 'H12-0235-10', 'H14-0599-8', 'H11-0383-2', 'H12-0330-4', 
                        'H18-0130-2', '2538-11O', 'H12-0330-1', '2318-1I', 'H17-0285-8', 'H16-0685-1', '2608-I9', 'H16-0223-4', 'H12-0235-9', '2241-5C', 
                        '2608-I5', '2166-1E', '2538-11M', '2571-1D', 'H13-0307-6', 'H11-0330-3', 'H17-0458-11', 'H16-0588-13', 'H17-0221-7', 'H08-0377-6', 
                        '2017-1B', '2017-1C', 'H09-0259-3', 'H12-0305-4', '2166-1F', 'H09-0259-1', 'H11-0330-4', '2323-5G', '2420-1H', '2405-1D', 
                        '2325-1F', '2513-1H', 'H14-0608-3', 'H17-0221-1', 'H16-0206-3', 'H18-0130-7', 'H18-0360-7', 'H18-0052-2', 'H18-0518-10', 'H17-0179-9', 
                        'H18-0454-5', '2571-1D', 'H10-0526-3', 'H14-0290-14', 'H13-0307-1', 'H10-0271-3', 'H18-0201-9', 'H15-0129-1', 'H18-0712-3', 'H18-0656-2', 
                        'H18-0255-6', 'H18-0331-10', 'H17-0285-1', 'H10-0702-2', '2513-1C', 'H13-0583-8', '2571-1C', 'H18-0255-9', 'H16-0223-6', 'H18-0130-5', 
                        'H18-0165-2', 'H15-0519-4', 'H10-0526-4', 'H18-0271-8', 'H17-0221-5', '2618-E2', 'H10-0553-3', 'H18-0455-4', 'H13-0583-3', 'H18-0201-15', 
                        '2241-1C', '2248-4E', 'H16-0588-10', 'H18-0271-3', 'H18-0331-11', 'H14-0599-3', '2420-1D', '2318-1A', 'H16-0206-12', 'H16-0223-9', 
                        '2532-9C', 'H17-0179-7', '2538-11P', '2621-D5', '2248-1E', 'H12-0330-2', 'H12-0235-6', '2582-1B', 'H17-0458-5', '2166-1B', 
                        'H18-0205-2', '2621-D3', 'H13-0307-7', '2538-11B', '2493-1F', '2323-1D', 'H13-0583-5', '2017-1G', 'H14-0290-12', 'H18-0254-2', 
                        'H18-0518-7', '2405-2E', '2485-1F', '697s-1H']
    
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
    plot_name = "AdjacentNormalCellFraction"
    plot_path = os.path.join(stat_result_dir, plot_name + args.plot_format)
    plt.savefig(plot_path, transparent=False, dpi=300)