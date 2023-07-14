# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz
import pandas as pd
import numpy as np
from datetime import datetime
from sklearn.impute import SimpleImputer


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])
    parser.add_argument("--cellphenotype_dir",      type=str,       default="CellPhenotyping")    
    parser.add_argument("--feature_dir",            type=str,       default="FeatureAnalysis")    

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    feature_root_dir = os.path.join(dataset_dir, args.feature_dir)
    # load interaction sigval
    phenotype_dir = os.path.join(dataset_dir, args.cellphenotype_dir)  
    interaction_path = os.path.join(phenotype_dir, "DelaunayInteraction", "DelaunayCN8Interaction-Sigval.csv")
    interaction_df = pd.read_csv(interaction_path)
    # print("Total number of pairs is: {}".format(len(interaction_df)))

    roi_names = interaction_df["roi_names"].tolist()
    unique_roi_names = list(set(roi_names))
    num_rois = len(unique_roi_names)
    celltype_lst = ["1", "2", "3", "4", "5", "6", "7", "8"]
    ct_cn_dict = {"1":"Undefined-CN1", "2":"Epithelial1-CN2", "3":"Proliferating-CN3", "4":"Epithelial2-CN4",
                  "5":"Endothelial-CN5", "6":"Fibroblast-CN6", "7":"Macrophage-CN7", "8":"PanImmune-CN8"}    
    interaction_lst = []
    for from_type in celltype_lst:
        for to_type in celltype_lst:
            interaction_lst.append(ct_cn_dict[from_type] + "-" + ct_cn_dict[to_type])
    
    # Initiate data frame
    roi_df = pd.DataFrame(index = range(num_rois), columns = range(len(interaction_lst)))
    roi_df.columns = interaction_lst
    for ind in interaction_df.index:
        roi_name = interaction_df["roi_names"][ind]
        from_phenotype = interaction_df["from_phenotypes"][ind]
        to_phenotype = interaction_df["to_phenotypes"][ind]
        interaction_name = ct_cn_dict[str(from_phenotype)] + "-" + ct_cn_dict[str(to_phenotype)]
        row_ind = unique_roi_names.index(roi_name)
        roi_df[interaction_name][row_ind] = interaction_df["sig_vals"][ind]
    roi_arr = roi_df.to_numpy()
    # Create cell interaction features
    interaction_fea_df = pd.DataFrame(data = roi_arr, columns = interaction_lst)
    interaction_fea_df.insert(0, "ROI_ID", unique_roi_names, True)
    # Save features
    cell_fea_path = os.path.join(feature_root_dir, "DelaunayCN8InteractionFeas.csv")
    interaction_fea_df.to_csv(cell_fea_path, index=False)