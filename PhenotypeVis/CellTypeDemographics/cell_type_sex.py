# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt


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
    patient_info_path = os.path.join(metadata_dir, "Patient_Info.xlsx")
    patient_df = pd.read_excel(patient_info_path)

    # load patient gender information
    patient_lst = patient_df["PatientID"].tolist()
    gender_lst = patient_df["Gender"].tolist()
    pat_sex_dict = {patient:gender for patient, gender in zip(patient_lst, gender_lst)}


    # load cell phenotype information
    cell_phenotype_path = os.path.join(roi_root_dir, args.steinbock_dir, "cell_phenotypes.json")
    cell_phenotype_dict = None
    with open(cell_phenotype_path) as fp:
        cell_phenotype_dict = json.load(fp)
    cell_ids = [ele for ele in cell_phenotype_dict.keys()]
    cell_phenotypes = [ele for ele in cell_phenotype_dict.values()]
    cell_rois = [ele[:ele.rfind("_")] for ele in cell_ids]    

   major_type_dict = {"Epithelial-Cell": "Epithelial", "Endothelial-Cell": "Stromal", "Fibroblast": "Stromal", "Other-Immune": "Immune",
                       "NK-Cell": "Immune", "B-Cell": "Immune", "CD4-T-Cell": "Immune", "CD8-T-Cell": "Immune", "T-Reg-Cell": "Immune",
                       "Dendritic-Cell": "Immune", "Neutrophil": "Immune", "Monocyte": "Immune", "Macrophage": "Immune", "MDSC": "Immune"}
