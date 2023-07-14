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


    # load roi information
    roi_info_path = os.path.join(metadata_dir, "ROI_Info.xlsx")
    roi_df = pd.read_excel(roi_info_path)
    roi_lst = roi_df["ROI_ID"].tolist()            

    # load cell phenotype information
    cell_phenotype_path = os.path.join(roi_root_dir, args.steinbock_dir, "cell_phenotypes.json")
    cell_phenotype_dict = None
    with open(cell_phenotype_path) as fp:
        cell_phenotype_dict = json.load(fp)
    cell_ids = [ele for ele in cell_phenotype_dict.keys()]
    cell_phenotypes = [ele for ele in cell_phenotype_dict.values()]
    cell_rois = [ele[:ele.rfind("_")] for ele in cell_ids]    

    major_type_lst = ["Epithelial-Cell", "Endothelial-Cell", "Fibroblast", "Other-Immune",
                        "NK-Cell", "B-Cell", "CD4-T-Cell", "CD8-T-Cell", "T-Reg-Cell",
                        "Dendritic-Cell", "Neutrophil", "Monocyte", "Macrophage", "MDSC"]
    roi_cell_dict = {}
    for roi, phenotype in zip(cell_rois, cell_phenotypes):
        if roi not in roi_cell_dict:
            roi_cell_dict[roi] = [phenotype]
        else:
            roi_cell_dict[roi].append(phenotype)

    epithelial_ratio_lst, endothelial_raio_lst, fibrobalst_ratio_lst = [], [], []
    otherimmune_ratio_lst, nkcell_ratio_lst, bcell_ratio_lst = [], [], []
    cd4tcell_ratio_lst, cd8tcell_ratio_lst, tregcell_raito_lst = [], [], []
    dendritic_ratio_lst, neutrophil_ratio_lst, monocyte_raito_lst = [], [], []
    macrophage_ratio_lst, mdsc_ratio_lst = [], []

    for roi in roi_lst:
        cell_lst = roi_cell_dict[roi]
        epithelial_ratio_lst.append(np.sum([ele == "Epithelial-Cell" for ele in cell_lst]) * 1.0 / len(cell_lst))
        endothelial_raio_lst.append(np.sum([ele == "Endothelial-Cell" for ele in cell_lst]) * 1.0 / len(cell_lst))
        fibrobalst_ratio_lst.append(np.sum([ele == "Fibroblast" for ele in cell_lst]) * 1.0 / len(cell_lst))
        otherimmune_ratio_lst.append(np.sum([ele == "Other-Immune" for ele in cell_lst]) * 1.0 / len(cell_lst))
        nkcell_ratio_lst.append(np.sum([ele == "NK-Cell" for ele in cell_lst]) * 1.0 / len(cell_lst))
        bcell_ratio_lst.append(np.sum([ele == "B-Cell" for ele in cell_lst]) * 1.0 / len(cell_lst))
        cd4tcell_ratio_lst.append(np.sum([ele == "CD4-T-Cell" for ele in cell_lst]) * 1.0 / len(cell_lst))
        cd8tcell_ratio_lst.append(np.sum([ele == "CD8-T-Cell" for ele in cell_lst]) * 1.0 / len(cell_lst))
        tregcell_raito_lst.append(np.sum([ele == "T-Reg-Cell" for ele in cell_lst]) * 1.0 / len(cell_lst))
        dendritic_ratio_lst.append(np.sum([ele == "Dendritic-Cell" for ele in cell_lst]) * 1.0 / len(cell_lst))
        neutrophil_ratio_lst.append(np.sum([ele == "Neutrophil" for ele in cell_lst]) * 1.0 / len(cell_lst))
        monocyte_raito_lst.append(np.sum([ele == "Monocyte" for ele in cell_lst]) * 1.0 / len(cell_lst))
        macrophage_ratio_lst.append(np.sum([ele == "Macrophage" for ele in cell_lst]) * 1.0 / len(cell_lst))
        mdsc_ratio_lst.append(np.sum([ele == "MDSC" for ele in cell_lst]) * 1.0 / len(cell_lst))
    mean_ratios = [np.mean(epithelial_ratio_lst), np.mean(endothelial_raio_lst), np.mean(fibrobalst_ratio_lst), np.mean(otherimmune_ratio_lst),
        np.mean(nkcell_ratio_lst), np.mean(bcell_ratio_lst), np.mean(cd4tcell_ratio_lst), np.mean(cd8tcell_ratio_lst), np.mean(tregcell_raito_lst),
        np.mean(dendritic_ratio_lst), np.mean(neutrophil_ratio_lst), np.mean(monocyte_raito_lst), np.mean(macrophage_ratio_lst), np.mean(mdsc_ratio_lst)]
    mean_order = np.argsort(mean_ratios)[::-1]
    order_cell_types = [major_type_lst[order] for order in mean_order]

    cell_type_lst = []
    cell_ratio_lst = []
    for roi in roi_lst:
        cell_lst = roi_cell_dict[roi]
        for cell_type in major_type_lst:
            cell_type_lst.append(cell_type)
            cell_ratio_lst.append(np.sum([ele == cell_type for ele in cell_lst]) * 1.0 / len(cell_lst))

    cell_type_ratio_df = pd.DataFrame(list(zip(cell_type_lst, cell_ratio_lst)), 
        columns=["CellType", "CellRatio"])
    #define figure size
    sns.set(rc={"figure.figsize":(10, 12)})
    sns.set(font_scale=0.5)
    sns.catplot(data=cell_type_ratio_df, x="CellType", y="CellRatio", order=order_cell_types)
    plt.xticks(rotation=24)
    # save plot
    plot_name = "cell_type_ratios"
    plot_path = os.path.join(stat_result_dir, plot_name + args.plot_format)
    plt.savefig(plot_path, transparent=False, dpi=300)    