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

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    dataset_dir = os.path.join(args.data_root, args.data_set)
    slide_agg_dir = os.path.join(dataset_dir, args.aggregation_dir)
    slide_distprog_dir = os.path.join(slide_agg_dir, args.dist_prog_dir)
    vis_prog_dir = os.path.join(slide_distprog_dir, "VisProgression")
    vis_stage_prog_dir = os.path.join(vis_prog_dir, "StageProgression")
    if not os.path.exists(vis_stage_prog_dir):
        os.makedirs(vis_stage_prog_dir)

    # save lesion information
    core_distprog_dict = None
    lesion_core_distprog_path = os.path.join(slide_distprog_dir, "lesion_roi_info.json")
    with open(lesion_core_distprog_path) as fp:
        core_distprog_dict = json.load(fp)
    print("There are {} lesions in total.".format(len(core_distprog_dict)))

    # interested feature list
    study_fea_lst = ["CD8-T-Cell-Proportion", "CD8-T-Cell-Density", "Macrophage-Proportion", "Macrophage-Density", 
                     "B-Cell-Proportion", "B-Cell-Density", "Endothelial-Cell-Proportion", "Endothelial-Cell-Density",
                     "CD8-T-Cell-Epithelial-Cell", "Macrophage-Epithelial-Cell", "Endothelial-Cell-Fibroblast", "Neutrophil-T-Reg-Cell"]
    for cur_fea in study_fea_lst:
        fig, ax = plt.subplots(nrows=1, ncols=1)
        for cur_lesion in core_distprog_dict.keys():
            lesion_dict = core_distprog_dict[cur_lesion]
            dist_lst = lesion_dict["CoreDist"]
            # print(dist_lst)
            fea_lst = lesion_dict[cur_fea]
            if lesion_dict["LesionStage"] == "AAH":
                ax.plot(dist_lst, fea_lst, color = "green", label = "AAH")
            elif lesion_dict["LesionStage"] == "AIS":
                ax.plot(dist_lst, fea_lst, color = "cyan", label = "AIS")
            elif lesion_dict["LesionStage"] == "MIA":
                ax.plot(dist_lst, fea_lst, color = "violet", label = "MIA")
            elif lesion_dict["LesionStage"] == "ADC":
                ax.plot(dist_lst, fea_lst, color = "tomato", label = "ADC") 
            else:
                print("Unknow stage: {}".format(lesion_dict["LesionStage"]))
        handles, labels = plt.gca().get_legend_handles_labels()
        newLabels, newHandles = [], []
        for handle, label in zip(handles, labels):
            if label not in newLabels:
                newLabels.append(label)
                newHandles.append(handle)
        label_lst = ["AAH", "AIS", "MIA", "ADC"]
        handle_lst = [newHandles[newLabels.index(ele)] for ele in label_lst]
        plt.legend(handle_lst, label_lst, loc = "upper right")
        plt.xlim(0, 25000)
        plt.ylim(0, 1.0)
        plt.xlabel("Distance to Lesion Core (" + u"\u03bcm)")
        plt.title("Lesion Evolution: {}".format(cur_fea))        
        # save plotting
        plot_path = os.path.join(vis_stage_prog_dir, cur_fea + ".png")
        plt.savefig(plot_path, transparent=False, dpi=300)