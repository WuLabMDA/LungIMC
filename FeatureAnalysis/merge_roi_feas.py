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
    parser.add_argument("--feature_dir",            type=str,       default="FeatureAnalysis")    

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set) 
    feature_root_dir = os.path.join(dataset_dir, args.feature_dir)

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

    # Joint Dist CT-CN
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

    roi_fea_df = pd.concat([ct_proportion_density_df, ct_state_df, ct_morph_df, ct_diversity_df,
                            cn_proportion_density_df, cn_state_df, cn_morph_df, cn_diversity_df,
                            joint_dist_ct_cn_df, interaction_delaunay_df], axis=1)
    # Merge features
    merge_roi_fea_path = os.path.join(feature_root_dir, "ROI_Fea_Merge.csv")
    roi_fea_df.to_csv(merge_roi_fea_path, index=False)