# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, itertools
import pandas as pd
import numpy as np
from datetime import datetime
from tensorly.decomposition import tucker


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--steinbock_dir",          type=str,       default="Steinbock")
    parser.add_argument("--feature_dir",            type=str,       default="FeatureAnalysis") 

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set) 
    feature_root_dir = os.path.join(dataset_dir, args.feature_dir)
    celltype_dir = os.path.join(dataset_dir, "CellPhenotyping", "Cell_CT_CN_Morphs")   

    # load features
    test_stage = "AAH"
    merge_roi_fea_path = os.path.join(feature_root_dir, "ROI_Fea_Merge.csv")
    merge_roi_fea_df = pd.read_csv(merge_roi_fea_path)
    stage_fea_df = merge_roi_fea_df[merge_roi_fea_df["ROI_Stage"] == test_stage]

    # filter ct-cn features
    celltype_lst = ["Epithelial-Cell", "Endothelial-Cell", "Fibroblast", "CD4-T-Cell", "CD8-T-Cell", 
                    "T-Reg-Cell", "B-Cell", "Macrophage", "Monocyte", "Dendritic-Cell", 
                    "Neutrophil", "MDSC", "NK-Cell", "Proliferating-Cell", "Undefined"]
    cellcn_lst = [str(ind+1) for ind in np.arange(10)]
    ct_cn_combs = []
    for cell_type in celltype_lst:
        for cn_type in cellcn_lst:
            ct_cn_comb = cell_type + "-CN" + cn_type
            ct_cn_combs.append(ct_cn_comb)  
    ct_cn_fea_df = stage_fea_df[ct_cn_combs]
    
    # construct tensor
    stage_roi_num = len(ct_cn_fea_df)
    ct_cn_fea_arr = ct_cn_fea_df.to_numpy()
    roi_fea_tensor = np.zeros((stage_roi_num, len(celltype_lst), len(cellcn_lst)), dtype=np.float32)
    for ind in range(stage_roi_num):
        ct_cn_fea = ct_cn_fea_arr[ind, :]
        roi_fea_tensor[ind, :, :] = ct_cn_fea.reshape(len(celltype_lst), len(cellcn_lst))

    # conduct tensor decomposition
    core, factors = tucker(roi_fea_tensor, rank=(2, 6, 6))
    # tensor CT
    ct_factor_path = os.path.join(feature_root_dir, "tensor_ct_factor.csv")
    ct_factor_df = pd.DataFrame(np.transpose(factors[1]), columns=celltype_lst)
    ct_factor_df.to_csv(ct_factor_path)
    # tensor CN
    cn_factor_path = os.path.join(feature_root_dir, "tensor_cn_factor.csv")
    cn_factor_df = pd.DataFrame(np.transpose(factors[2]), columns=cellcn_lst)
    cn_factor_df.to_csv(cn_factor_path)

    
    core1_lst = core[0, :, :].diagonal()
    core2_lst = core[1, :, :].diagonal()   
    # print("Core 1 diagnoal:")
    # print("Core 2 diagnoal:")


    ct_cn_core_path = os.path.join(feature_root_dir, "tensor_core_diagnoals.csv")
    ct_cn_core_df = pd.DataFrame(list(zip(core1_lst, core2_lst)), columns=["Core1", "Core2"])
    ct_cn_core_df.to_csv(ct_cn_core_path, index=False)
