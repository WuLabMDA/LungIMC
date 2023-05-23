# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz
from datetime import datetime
import pandas as pd
import numpy as np


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])
    parser.add_argument("--metadata_dir",           type=str,       default="Metadata")   
    parser.add_argument("--feature_dir",            type=str,       default="FeatureAnalysis")    
    parser.add_argument("--recur_dir",              type=str,       default="Recurrence")


    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    metadata_dir = os.path.join(dataset_dir, args.metadata_dir) 
    feature_root_dir = os.path.join(dataset_dir, args.feature_dir)

    slide_recur_dir = os.path.join(dataset_dir, args.recur_dir)
    recurrence_json_path = os.path.join(slide_recur_dir, "LesionRecurrence.json")
    lesion_recur_dict = None
    with open(recurrence_json_path) as fp:
        lesion_recur_dict = json.load(fp)
    lesion_lst = [ele for ele in lesion_recur_dict.keys()]

    # CT Proportion/Density
    ct_proportion_density_fea_path = os.path.join(feature_root_dir, "CT_ProportionDensityFeas.csv")
    ct_proportion_density_df = pd.read_csv(ct_proportion_density_fea_path)
    # CT Morphology
    ct_morph_fea_path = os.path.join(feature_root_dir, "CT_MorphFeas.csv")
    ct_morph_df = pd.read_csv(ct_morph_fea_path)
    ct_morph_df = ct_morph_df.set_index("ROI_ID")
    ct_morph_df  = ct_morph_df.reindex(index=ct_proportion_density_df["ROI_ID"])
    ct_morph_df = ct_morph_df.reset_index()
    ct_morph_df = ct_morph_df.iloc[:, 1:]
    # CT States
    ct_state_fea_path = os.path.join(feature_root_dir, "CT_StateFeas.csv")
    ct_state_df = pd.read_csv(ct_state_fea_path)
    ct_state_df = ct_state_df.set_index("ROI_ID")
    ct_state_df  = ct_state_df.reindex(index=ct_proportion_density_df["ROI_ID"])
    ct_state_df = ct_state_df.reset_index()
    ct_state_df = ct_state_df.iloc[:, 1:]    
    # Interaction Delaunay
    interaction_delaunay_fea_path = os.path.join(feature_root_dir, "InteractionDelaunay50Feas.csv")
    interaction_delaunay_df = pd.read_csv(interaction_delaunay_fea_path)
    interaction_delaunay_df = interaction_delaunay_df.set_index("ROI_ID")
    interaction_delaunay_df  = interaction_delaunay_df.reindex(index=ct_proportion_density_df["ROI_ID"])
    interaction_delaunay_df = interaction_delaunay_df.reset_index()
    ct_interaction_df = interaction_delaunay_df.iloc[:, 1:]

    # CN Proportion/Density
    cn_proportion_density_fea_path = os.path.join(feature_root_dir, "CN_ProportionDensityFeas.csv")
    cn_proportion_density_df = pd.read_csv(cn_proportion_density_fea_path)
    cn_proportion_density_df = cn_proportion_density_df.set_index("ROI_ID")
    cn_proportion_density_df  = cn_proportion_density_df.reindex(index=ct_proportion_density_df["ROI_ID"])
    cn_proportion_density_df = cn_proportion_density_df.reset_index()
    cn_proportion_density_df = cn_proportion_density_df.iloc[:, 1:]          
    # CN Morphology
    cn_morph_fea_path = os.path.join(feature_root_dir, "CN_MorphFeas.csv")
    cn_morph_df = pd.read_csv(cn_morph_fea_path)
    cn_morph_df = cn_morph_df.set_index("ROI_ID")
    cn_morph_df  = cn_morph_df.reindex(index=ct_proportion_density_df["ROI_ID"])
    cn_morph_df = cn_morph_df.reset_index()
    cn_morph_df = cn_morph_df.iloc[:, 1:]
    # CN States
    cn_state_fea_path = os.path.join(feature_root_dir, "CN_StateFeas.csv")
    cn_state_df = pd.read_csv(cn_state_fea_path)
    cn_state_df = cn_state_df.set_index("ROI_ID")
    cn_state_df = cn_state_df.reindex(index=ct_proportion_density_df["ROI_ID"])
    cn_state_df = cn_state_df.reset_index()
    cn_state_df = cn_state_df.iloc[:, 1:]      
    # CN Interaction 
    cn_interaction_fea_path = os.path.join(feature_root_dir, "DelaunayCN8InteractionFeas.csv")
    cn_interaction_delaunay_df = pd.read_csv(cn_interaction_fea_path)
    cn_interaction_delaunay_df = cn_interaction_delaunay_df.set_index("ROI_ID")
    cn_interaction_delaunay_df  = cn_interaction_delaunay_df.reindex(index=ct_proportion_density_df["ROI_ID"])
    cn_interaction_delaunay_df = cn_interaction_delaunay_df.reset_index()
    cn_interaction_df = cn_interaction_delaunay_df.iloc[:, 1:]

    # Aggregate features
    roi_fea_df = pd.concat([ct_proportion_density_df, ct_morph_df, ct_state_df, ct_interaction_df,
                            cn_proportion_density_df, cn_morph_df, cn_state_df, cn_interaction_df], axis=1)
    
    # Insert ROI stage information
    roi_info_path = os.path.join(metadata_dir, "ROI_Info.xlsx")
    roi_info_df = pd.read_excel(roi_info_path)
    roi_id_lst = [ele for ele in roi_info_df["ROI_ID"].tolist()]
    roi_loc_lst = [ele for ele in roi_info_df["ROI_Location"].tolist()]
    roi_diag_lst = [ele for ele in roi_info_df["ROI_Diag"].tolist()]

    # obtain roi stages
    roi_stage_lst = []
    for loc, diag in zip(roi_loc_lst, roi_diag_lst):
        if loc == "Tumor":
            roi_stage_lst.append(diag)
        else:
            roi_stage_lst.append("Normal")

    # construct ROI data frame
    roi_stage_df = pd.DataFrame(list(zip(roi_id_lst, roi_stage_lst)), columns=["ROI_ID", "ROI_Stage"])
    roi_stage_df = roi_stage_df.set_index("ROI_ID")
    roi_stage_df  = roi_stage_df.reindex(index=ct_proportion_density_df["ROI_ID"])
    roi_stage_df = roi_stage_df.reset_index()
    roi_stage_lst = [ele for ele in roi_stage_df["ROI_Stage"].tolist()]

    # insert stage column and sort
    roi_fea_df.insert(loc=1, column="ROI_Stage", value=roi_stage_lst)
    # # print ROI number counts
    # print("In total, there are {} ROIs.".format(roi_fea_df.shape[0]))

    # Feature processing
    fea_lst = list(roi_fea_df.columns[2:])
    for cur_fea in fea_lst:
        # z-scoring feature
        roi_fea_df[cur_fea] = (roi_fea_df[cur_fea] - roi_fea_df[cur_fea].mean()) / roi_fea_df[cur_fea].std(ddof=0)
        # clipping (-3.0, 3.0)
        roi_fea_df[cur_fea] = roi_fea_df[cur_fea].clip(-3.0, 3.0)

    # filtering lesion recurrence
    roi_fea_df = roi_fea_df[roi_fea_df["ROI_Stage"] != "Normal"]
    roi_recurrence_inds = []
    for ind, roi_name in enumerate(roi_fea_df["ROI_ID"]):
        lesion_name = roi_name[:-7]
        if lesion_name in lesion_lst:
            roi_recurrence_inds.append(ind)
    roi_fea_df = roi_fea_df.iloc[roi_recurrence_inds, :]

    print("There are {} ROIs with Recurrence status.".format(roi_fea_df.shape[0]))
    roi_recur_lst = [] 
    for ind, roi_name in enumerate(roi_fea_df["ROI_ID"]):
        lesion_name = roi_name[:-7]
        roi_recur_lst.append(lesion_recur_dict[lesion_name])
    roi_fea_df.insert(loc = 2, column = "Recur", value = roi_recur_lst)

    unique_lesions = set([ele[:-7] for ele in roi_fea_df["ROI_ID"].tolist()])
    print("There are {} unique lesions with recur information".format(len(unique_lesions)))
    # Merge features
    agg_roi_fea_path = os.path.join(slide_recur_dir, "lesion_roi_feas.csv")
    roi_fea_df.to_csv(agg_roi_fea_path, index=False)