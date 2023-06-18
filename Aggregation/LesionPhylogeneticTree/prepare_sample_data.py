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
    parser.add_argument("--aggregation_dir",        type=str,       default="Aggregation")    


    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    metadata_dir = os.path.join(dataset_dir, args.metadata_dir) 
    feature_root_dir = os.path.join(dataset_dir, args.feature_dir)
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)
    if not os.path.exists(slide_agg_dir):
        os.makedirs(slide_agg_dir)

    fea_num_lst = []
    # CT Proportion/Density
    ct_proportion_density_fea_path = os.path.join(feature_root_dir, "CT_ProportionDensityFeas.csv")
    ct_proportion_density_df = pd.read_csv(ct_proportion_density_fea_path)
    fea_num_lst.append(ct_proportion_density_df.shape[1] - 1)
    # CT States
    ct_state_fea_path = os.path.join(feature_root_dir, "CT_StateFeas.csv")
    ct_state_df = pd.read_csv(ct_state_fea_path)
    fea_num_lst.append(ct_state_df.shape[1] - 1)
    ct_state_df = ct_state_df.set_index("ROI_ID")
    ct_state_df  = ct_state_df.reindex(index=ct_proportion_density_df["ROI_ID"])
    ct_state_df = ct_state_df.reset_index()
    ct_state_df = ct_state_df.iloc[:, 1:]
    # CT Morphology
    ct_morph_fea_path = os.path.join(feature_root_dir, "CT_MorphFeas.csv")
    ct_morph_df = pd.read_csv(ct_morph_fea_path)
    fea_num_lst.append(ct_morph_df.shape[1] - 1)
    ct_morph_df = ct_morph_df.set_index("ROI_ID")
    ct_morph_df  = ct_morph_df.reindex(index=ct_proportion_density_df["ROI_ID"])
    ct_morph_df = ct_morph_df.reset_index()
    ct_morph_df = ct_morph_df.iloc[:, 1:]
    # CT Diveristy
    ct_diveristy_fea_path = os.path.join(feature_root_dir, "CT_DiversityFeas.csv")
    ct_diversity_df = pd.read_csv(ct_diveristy_fea_path)
    fea_num_lst.append(ct_diversity_df.shape[1] - 1)
    ct_diversity_df = ct_diversity_df.set_index("ROI_ID")
    ct_diversity_df  = ct_diversity_df.reindex(index=ct_proportion_density_df["ROI_ID"])
    ct_diversity_df = ct_diversity_df.reset_index()
    ct_diversity_df = ct_diversity_df.iloc[:, 1:]  

    # CN Proportion/Density
    cn_proportion_density_fea_path = os.path.join(feature_root_dir, "CN_ProportionDensityFeas.csv")
    cn_proportion_density_df = pd.read_csv(cn_proportion_density_fea_path)
    fea_num_lst.append(cn_proportion_density_df.shape[1] - 1)
    cn_proportion_density_df = cn_proportion_density_df.set_index("ROI_ID")
    cn_proportion_density_df  = cn_proportion_density_df.reindex(index=ct_proportion_density_df["ROI_ID"])
    cn_proportion_density_df = cn_proportion_density_df.reset_index()
    cn_proportion_density_df = cn_proportion_density_df.iloc[:, 1:]        
    # CN States
    cn_state_fea_path = os.path.join(feature_root_dir, "CN_StateFeas.csv")
    cn_state_df = pd.read_csv(cn_state_fea_path)
    fea_num_lst.append(cn_state_df.shape[1] - 1)
    cn_state_df = cn_state_df.set_index("ROI_ID")
    cn_state_df = cn_state_df.reindex(index=ct_proportion_density_df["ROI_ID"])
    cn_state_df = cn_state_df.reset_index()
    cn_state_df = cn_state_df.iloc[:, 1:]    
    # CN Morphology
    cn_morph_fea_path = os.path.join(feature_root_dir, "CN_MorphFeas.csv")
    cn_morph_df = pd.read_csv(cn_morph_fea_path)
    fea_num_lst.append(cn_morph_df.shape[1] - 1)
    cn_morph_df = cn_morph_df.set_index("ROI_ID")
    cn_morph_df  = cn_morph_df.reindex(index=ct_proportion_density_df["ROI_ID"])
    cn_morph_df = cn_morph_df.reset_index()
    cn_morph_df = cn_morph_df.iloc[:, 1:]    
    # CN Diveristy
    cn_diveristy_fea_path = os.path.join(feature_root_dir, "CN_DiversityFeas.csv")
    cn_diversity_df = pd.read_csv(cn_diveristy_fea_path)
    fea_num_lst.append(cn_diversity_df.shape[1] - 1)
    cn_diversity_df = cn_diversity_df.set_index("ROI_ID")
    cn_diversity_df  = cn_diversity_df.reindex(index=ct_proportion_density_df["ROI_ID"])
    cn_diversity_df = cn_diversity_df.reset_index()
    cn_diversity_df = cn_diversity_df.iloc[:, 1:]  

    # Joint CT-CN Dist 
    joint_dist_ct_cn_fea_path = os.path.join(feature_root_dir, "JointDistCTCNFeas.csv")
    joint_dist_ct_cn_df = pd.read_csv(joint_dist_ct_cn_fea_path)
    fea_num_lst.append(joint_dist_ct_cn_df.shape[1] - 1)
    joint_dist_ct_cn_df = joint_dist_ct_cn_df.set_index("ROI_ID")
    joint_dist_ct_cn_df  = joint_dist_ct_cn_df.reindex(index=ct_proportion_density_df["ROI_ID"])
    joint_dist_ct_cn_df = joint_dist_ct_cn_df.reset_index()
    joint_dist_ct_cn_df = joint_dist_ct_cn_df.iloc[:, 1:]    

    # Interaction Delaunay
    interaction_delaunay_fea_path = os.path.join(feature_root_dir, "InteractionDelaunay50Feas.csv")
    interaction_delaunay_df = pd.read_csv(interaction_delaunay_fea_path)
    fea_num_lst.append(interaction_delaunay_df.shape[1] - 1)
    interaction_delaunay_df = interaction_delaunay_df.set_index("ROI_ID")
    interaction_delaunay_df  = interaction_delaunay_df.reindex(index=ct_proportion_density_df["ROI_ID"])
    interaction_delaunay_df = interaction_delaunay_df.reset_index()
    interaction_delaunay_df = interaction_delaunay_df.iloc[:, 1:]

    # Aggregate features
    roi_fea_df = pd.concat([ct_proportion_density_df, ct_morph_df, ct_state_df, cn_proportion_density_df, cn_morph_df,
                            interaction_delaunay_df], axis=1)
    
    
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
            roi_stage_lst.append(loc)
    # construct ROI data frame
    roi_stage_df = pd.DataFrame(list(zip(roi_id_lst, roi_stage_lst)), columns=["ROI_ID", "ROI_Stage"])
    roi_stage_df = roi_stage_df.set_index("ROI_ID")
    roi_stage_df  = roi_stage_df.reindex(index=ct_proportion_density_df["ROI_ID"])
    roi_stage_df = roi_stage_df.reset_index()
    roi_stage_lst = [ele for ele in roi_stage_df["ROI_Stage"].tolist()]

    # filter patient
    roi_fea_df.insert(loc=1, column="ROI_Stage", value=roi_stage_lst)
    # normalize all features between 0.0 to 1.0
    roi_fea_columns = [ele for ele in roi_fea_df.columns.tolist()]
    roi_fea_names = roi_fea_columns[2:]
    for fea in roi_fea_names:
        roi_fea_df[fea] = (roi_fea_df[fea] - roi_fea_df[fea].mean()) / roi_fea_df[fea].std(ddof=0)
        roi_fea_df[fea] = (roi_fea_df[fea].clip(-2.0, 2.0) + 2.0) / 4.0

    # filter AdjacentNormal & DistantNormal
    roi_fea_df = roi_fea_df[roi_fea_df["ROI_Stage"] != "AdjacentNormal"]

    # filter patient
    # interested_pid = "2405"
    # interested_pid = "H18-0271"
    interested_pid = "H16-0223"
    lesion_names = [ele[:-7] for ele in roi_fea_df["ROI_ID"].tolist()]
    patient_names = [ele[:ele.rindex("-")] for ele in lesion_names]
    patient_inds = []
    for ind, cur_name in enumerate(patient_names):
        if cur_name == interested_pid:
            patient_inds.append(ind)
    pat_fea_df = roi_fea_df.iloc[patient_inds]
    
    # save features
    phylo_dir = os.path.join(slide_agg_dir, "Phylo")
    if not os.path.exists(phylo_dir):
        os.makedirs(phylo_dir)
    pat_roi_fea_path = os.path.join(phylo_dir, "{}_roi_feas.csv".format(interested_pid))
    pat_fea_df.to_csv(pat_roi_fea_path, index=False)