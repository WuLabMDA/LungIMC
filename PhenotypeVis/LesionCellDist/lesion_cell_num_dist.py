# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt


def set_args():
    parser = argparse.ArgumentParser(description = "Assess cell fraction")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--meta_dir",               type=str,       default="Metadata")
    parser.add_argument("--cellphenotype_dir",      type=str,       default="CellPhenotyping")    
    parser.add_argument("--result_dir",             type=str,       default="Results")
    parser.add_argument("--plot_format",            type=str,       default=".png", choices=[".png", ".pdf"])        

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    # directory setting
    roi_root_dir = os.path.join(args.data_root, args.data_type)
    metadata_dir = os.path.join(args.data_root, args.meta_dir)
    cellphenotype_dir = os.path.join(args.data_root, args.cellphenotype_dir)
    stat_result_dir = os.path.join(args.data_root, args.result_dir, "DistStats")
    if not os.path.exists(stat_result_dir):
        os.makedirs(stat_result_dir)    
        
    lesion_merge_lst = ['2493-6M', '2513-1G', 'H18-0130-6', '2571-4B', 'H09-0054-3', 'H11-0383-6', 'H12-0235-10', 'H14-0599-8', 'H11-0383-2', 'H12-0330-4', 
                        'H18-0130-2', '2538-11O', 'H12-0330-1', '2318-1I', 'H17-0285-8', 'H16-0685-1', '2608-I9', 'H16-0223-4', 'H12-0235-9', '2241-5C', 
                        '2608-I5', '2166-1E', '2538-11M', '2571-1D', 'H13-0307-6', 'H11-0330-3', 'H17-0458-11', 'H16-0588-13', 'H17-0221-7', 'H08-0377-6', 
                        '2017-1B', '2017-1C', 'H09-0259-3', 'H12-0305-4', '2166-1F', 'H09-0259-1', 'H11-0330-4', '2323-5G', '2420-1H', '2405-1D', 
                        '2325-1F', '2513-1H', 'H14-0608-3', 'H17-0221-1', 'H16-0206-3', 'H18-0130-7', 'H18-0360-7', 'H18-0052-2', 'H18-0518-10', 'H17-0179-9', 
                        'H18-0454-5', '2571-1D', 'H10-0526-3', 'H14-0290-14', 'H13-0307-1', 'H10-0271-3', 'H18-0201-9', 'H15-0129-1', 'H18-0712-3', 'H18-0656-2', 
                        'H18-0255-6', 'H18-0331-10', 'H17-0285-1', 'H10-0702-2', '2513-1C', 'H13-0583-8', '2571-1C', 'H18-0255-9', 'H16-0223-6', 'H18-0130-5', 
                        'H18-0165-2', 'H15-0519-4', 'H10-0526-4', 'H18-0271-8', 'H17-0221-5', '2618-E2', 'H10-0553-3', 'H18-0455-4', 'H13-0583-3', 'H18-0201-15', 
                        '2241-1C', '2248-4E', 'H16-0588-10', 'H18-0271-3', 'H18-0331-11', 'H14-0599-3', '2420-1D', '2318-1A', 'H16-0206-12', 'H16-0223-9', 
                        '2532-9C', 'H17-0179-7', '2538-11P', '2621-D5', '2248-1E', 'H12-0330-2', 'H12-0235-6', '2582-1B', 'H17-0458-5', '2166-1B', 
                        'H18-0205-2', '2621-D3', 'H13-0307-7', '2538-11B', '2493-1F', '2323-1D', 'H13-0583-5', '2017-1G', 'H14-0290-12', 'H18-0254-2', 
                        'H18-0518-7', '2405-2E', '2485-1F', '697s-1H']
    lesion_cell_num_dict = {}
    lesion_roi_dict = {}
    for leion_name in lesion_merge_lst:
        lesion_cell_num_dict[leion_name] = 0
        lesion_roi_dict[leion_name] = set()
    
    # load cell phenotype information
    cell_phenotype_path = os.path.join(cellphenotype_dir, "cell_ct_cn_morphs.csv")
    cell_phenotype_df = pd.read_csv(cell_phenotype_path)
    cell_id_lst = cell_phenotype_df["cell_id"]
    for cell_id in cell_id_lst:
        lesion_id = cell_id[:cell_id.rfind("-ROI")]
        lesion_roi = cell_id[:cell_id.rfind("_")]
        if lesion_id in lesion_merge_lst:
            lesion_cell_num_dict[lesion_id] += 1
            lesion_roi_dict[lesion_id].add(lesion_roi)
    
    average_num_lst = []
    for leion_name in lesion_merge_lst:
        lesion_cell_num = lesion_cell_num_dict[leion_name]
        lesion_roi_num = len(lesion_roi_dict[leion_name])
        average_num = int(np.floor(0.5 + lesion_cell_num * 1.0 / lesion_roi_num))
        average_num_lst.append(average_num)

    lesion_cell_df = pd.DataFrame({"LesionName": lesion_merge_lst, "CellNum": average_num_lst})
    lesion_cell_df.plot.bar(x="LesionName", y="CellNum", width=0.8, figsize=(25, 12))
    plot_name = "LesionTotalCellNum"
    plot_path = os.path.join(stat_result_dir, plot_name + args.plot_format)
    plt.savefig(plot_path, transparent=False, dpi=300)