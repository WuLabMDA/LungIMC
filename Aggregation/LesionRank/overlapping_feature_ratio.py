# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, pickle, math
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])                        
    parser.add_argument("--aggregation_dir",        type=str,       default="Aggregation")
    parser.add_argument("--rank_dir",               type=str,       default="Rank")
    parser.add_argument("--dominant",               type=str,       default="Heavy", choices=["Heavy", "Never"])
    parser.add_argument("--path_stage",             type=str,       default="AAH", choices=["AAH", "AIS", "MIA", "ADC"])
    parser.add_argument("--plot_format",            type=str,       default=".png", choices=[".png", ".pdf"])    

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)
    # load ROI information
    roi_rank_dir = os.path.join(slide_agg_dir, args.rank_dir, "ROI")
    roi_rank_volcano_path = os.path.join(roi_rank_dir, "{}_ROI_VolcanoFeatures.csv".format(args.path_stage))
    roi_volcano_df = pd.read_csv(roi_rank_volcano_path)
    roi_volcano_df = roi_volcano_df[roi_volcano_df["Pvalue"] < 0.05]
    if args.dominant == "Heavy":
        roi_volcano_df = roi_volcano_df[roi_volcano_df["Log2FC"] > 0.0]
    else:
        roi_volcano_df = roi_volcano_df[roi_volcano_df["Log2FC"] < 0.0]
    roi_features = [ele for ele in roi_volcano_df["Feature"].tolist()]
    # print("There are {} ROI features.".format(len(roi_features)))

    # load lesion information
    fea_num_lst = np.arange(len(roi_features)) + 1
    lesion_rank_dir = os.path.join(slide_agg_dir, args.rank_dir, "Lesion")

    # average
    avg_lesion_rank_name = "{}_Lesion_VolcanoFeatures_avg.csv".format(args.path_stage)
    avg_lesion_rank_volcano_path = os.path.join(lesion_rank_dir, avg_lesion_rank_name)
    avg_lesion_volcano_df = pd.read_csv(avg_lesion_rank_volcano_path)
    if args.dominant == "Heavy":
        avg_lesion_volcano_df = avg_lesion_volcano_df[avg_lesion_volcano_df["Log2FC"] > 0.0]
    else:
        avg_lesion_volcano_df = avg_lesion_volcano_df[avg_lesion_volcano_df["Log2FC"] < 0.0]
    avg_lesion_volcano_df.sort_values(by="Pvalue")
    avg_lesion_features = [ele for ele in avg_lesion_volcano_df["Feature"].tolist()]
    avg_overlap_ratios = []
    for fea_num in fea_num_lst:
        lesion_feas = avg_lesion_features[:fea_num]
        in_num = np.sum([ele in roi_features for ele in lesion_feas])
        avg_overlap_ratios.append(in_num * 1.0 / fea_num)
    avg_overlap_area = np.sum(avg_overlap_ratios) / len(avg_overlap_ratios)

    # min
    min_lesion_rank_name = "{}_Lesion_VolcanoFeatures_min.csv".format(args.path_stage)
    min_lesion_rank_volcano_path = os.path.join(lesion_rank_dir, min_lesion_rank_name)
    min_lesion_volcano_df = pd.read_csv(min_lesion_rank_volcano_path)
    if args.dominant == "Heavy":
        min_lesion_volcano_df = min_lesion_volcano_df[min_lesion_volcano_df["Log2FC"] > 0.0]
    else:
        min_lesion_volcano_df = min_lesion_volcano_df[min_lesion_volcano_df["Log2FC"] < 0.0]
    min_lesion_volcano_df.sort_values(by="Pvalue")
    min_lesion_features = [ele for ele in min_lesion_volcano_df["Feature"].tolist()]
    min_overlap_ratios = []
    for fea_num in fea_num_lst:
        lesion_feas = min_lesion_features[:fea_num]
        in_num = np.sum([ele in roi_features for ele in lesion_feas])
        min_overlap_ratios.append(in_num * 1.0 / fea_num)
    min_overlap_area = np.sum(min_overlap_ratios) / len(min_overlap_ratios)    

    # max
    max_lesion_rank_name = "{}_Lesion_VolcanoFeatures_max.csv".format(args.path_stage)
    max_lesion_rank_volcano_path = os.path.join(lesion_rank_dir, max_lesion_rank_name)
    max_lesion_volcano_df = pd.read_csv(max_lesion_rank_volcano_path)
    if args.dominant == "Heavy":
        max_lesion_volcano_df = max_lesion_volcano_df[max_lesion_volcano_df["Log2FC"] > 0.0]
    else:
        max_lesion_volcano_df = max_lesion_volcano_df[max_lesion_volcano_df["Log2FC"] < 0.0]
    max_lesion_volcano_df.sort_values(by="Pvalue")
    max_lesion_features = [ele for ele in max_lesion_volcano_df["Feature"].tolist()]
    max_overlap_ratios = []
    for fea_num in fea_num_lst:
        lesion_feas = max_lesion_features[:fea_num]
        in_num = np.sum([ele in roi_features for ele in lesion_feas])
        max_overlap_ratios.append(in_num * 1.0 / fea_num)
    max_overlap_area = np.sum(max_overlap_ratios) / len(max_overlap_ratios)  

    # bmin
    bmin_lesion_rank_name = "{}_Lesion_VolcanoFeatures_bmin.csv".format(args.path_stage)
    bmin_lesion_rank_volcano_path = os.path.join(lesion_rank_dir, bmin_lesion_rank_name)
    bmin_lesion_volcano_df = pd.read_csv(bmin_lesion_rank_volcano_path)
    if args.dominant == "Heavy":
        bmin_lesion_volcano_df = bmin_lesion_volcano_df[bmin_lesion_volcano_df["Log2FC"] > 0.0]
    else:
        bmin_lesion_volcano_df = bmin_lesion_volcano_df[bmin_lesion_volcano_df["Log2FC"] < 0.0]
    bmin_lesion_volcano_df.sort_values(by="Pvalue")
    bmin_lesion_features = [ele for ele in bmin_lesion_volcano_df["Feature"].tolist()]
    bmin_overlap_ratios = []
    for fea_num in fea_num_lst:
        lesion_feas = bmin_lesion_features[:fea_num]
        in_num = np.sum([ele in roi_features for ele in lesion_feas])
        bmin_overlap_ratios.append(in_num * 1.0 / fea_num)
    bmin_overlap_area = np.sum(bmin_overlap_ratios) / len(bmin_overlap_ratios)   

    # bmax
    bmax_lesion_rank_name = "{}_Lesion_VolcanoFeatures_bmax.csv".format(args.path_stage)
    bmax_lesion_rank_volcano_path = os.path.join(lesion_rank_dir, bmax_lesion_rank_name)
    bmax_lesion_volcano_df = pd.read_csv(bmax_lesion_rank_volcano_path)
    if args.dominant == "Heavy":
        bmax_lesion_volcano_df = bmax_lesion_volcano_df[bmax_lesion_volcano_df["Log2FC"] > 0.0]
    else:
        bmax_lesion_volcano_df = bmax_lesion_volcano_df[bmax_lesion_volcano_df["Log2FC"] < 0.0]
    bmax_lesion_volcano_df.sort_values(by="Pvalue")
    bmax_lesion_features = [ele for ele in bmax_lesion_volcano_df["Feature"].tolist()]
    bmax_overlap_ratios = []
    for fea_num in fea_num_lst:
        lesion_feas = bmax_lesion_features[:fea_num]
        in_num = np.sum([ele in roi_features for ele in lesion_feas])
        bmax_overlap_ratios.append(in_num * 1.0 / fea_num)
    bmax_overlap_area = np.sum(bmax_overlap_ratios) / len(bmax_overlap_ratios)     

    # plot 
    fig, ax = plt.subplots()
    ax.plot(fea_num_lst, avg_overlap_ratios, '--', linewidth=2, label="Average-{:.3f}".format(avg_overlap_area))
    ax.plot(fea_num_lst, max_overlap_ratios, '--', linewidth=2, label="Max-{:.3f}".format(max_overlap_area))
    ax.plot(fea_num_lst, min_overlap_ratios, '--', linewidth=2, label="Min-{:.3f}".format(min_overlap_area))
    ax.plot(fea_num_lst, bmax_overlap_ratios, '--', linewidth=2, label="B-Max-{:.3f}".format(bmax_overlap_area))
    ax.plot(fea_num_lst, bmin_overlap_ratios, '--', linewidth=2, label="B-Min-{:.3f}".format(bmin_overlap_area))    
    ax.legend()

    plot_name = "{}-Lesion-Feature-{}-Dominant".format(args.path_stage, args.dominant)
    plot_path = os.path.join(slide_agg_dir, args.rank_dir, plot_name + args.plot_format)
    plt.savefig(plot_path, transparent=False, dpi=300)    


