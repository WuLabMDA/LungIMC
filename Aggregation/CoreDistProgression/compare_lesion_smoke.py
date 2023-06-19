# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz, pickle
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])  
    parser.add_argument("--feature_dir",            type=str,       default="FeatureAnalysis")  
    parser.add_argument("--aggregation_dir",        type=str,       default="Aggregation")
    parser.add_argument("--dist_dir",               type=str,       default="Distance")
    parser.add_argument("--dist_prog_dir",          type=str,       default="DistProgression")
    parser.add_argument("--path_stage",             type=str,       default="AAH", choices=["AAH", "AIS", "MIA", "ADC"])

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)
    slide_distprog_dir = os.path.join(slide_agg_dir, args.dist_prog_dir)
    vis_prog_dir = os.path.join(slide_distprog_dir, "VisProgression")
    vis_smoke_prog_dir = os.path.join(vis_prog_dir, "{}-Smoke".format(args.path_stage))
    if not os.path.exists(vis_smoke_prog_dir):
        os.makedirs(vis_smoke_prog_dir)

    # save lesion information
    core_distprog_dict = None
    lesion_core_distprog_path = os.path.join(slide_distprog_dir, "lesion_roi_info.json")
    with open(lesion_core_distprog_path) as fp:
        core_distprog_dict = json.load(fp)
    print("There are {} lesions in total.".format(len(core_distprog_dict)))
    
    stage_distprog_dict = {}
    for cur_lesion in core_distprog_dict.keys():
        if core_distprog_dict[cur_lesion]["LesionStage"] == args.path_stage:
            stage_distprog_dict[cur_lesion] = core_distprog_dict[cur_lesion]
    print("{} has {} lesions.".format(args.path_stage, len(stage_distprog_dict)))

    # interested feature list
    study_fea_lst = ["CD8-T-Cell-Proportion", "CD8-T-Cell-Density", "Macrophage-Proportion", "Macrophage-Density", 
                     "B-Cell-Proportion", "B-Cell-Density", "Endothelial-Cell-Proportion", "Endothelial-Cell-Density",
                     "CD8-T-Cell-Epithelial-Cell", "Macrophage-Epithelial-Cell", "Endothelial-Cell-Fibroblast", "Neutrophil-T-Reg-Cell"]
    for cur_fea in study_fea_lst:
        fig, ax = plt.subplots()
        for cur_lesion in stage_distprog_dict.keys():
            lesion_dict = stage_distprog_dict[cur_lesion]
            dist_lst = lesion_dict["CoreDist"]
            # print(dist_lst)
            fea_lst = lesion_dict[cur_fea]
            if lesion_dict["SmokeStatus"] == "Never":
                ax.plot(dist_lst, fea_lst, color = "green", label = "Never Smoker")
            elif lesion_dict["SmokeStatus"] == "Heavy":
                ax.plot(dist_lst, fea_lst, color = "red", label = "Heavy Smoker")
        handles, labels = plt.gca().get_legend_handles_labels()
        newLabels, newHandles = [], []
        for handle, label in zip(handles, labels):
            if label not in newLabels:
                newLabels.append(label)
                newHandles.append(handle)
        label_lst = ["Never Smoker", "Heavy Smoker"]
        handle_lst = [newHandles[newLabels.index(ele)] for ele in label_lst]
        plt.legend(handle_lst, label_lst, loc = "upper right")
        # plt.ylim(0, 1.0)
        plt.xlabel("Distance to Lesion Core (" + u"\u03bcm)")
        plt.title("{} Never vs. Heavy Smoker: {}".format(args.path_stage, cur_fea)) 
        # save plotting
        plot_path = os.path.join(vis_smoke_prog_dir, cur_fea + ".png")
        plt.savefig(plot_path, transparent=False, dpi=300)
