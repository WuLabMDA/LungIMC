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
    parser.add_argument("--result_dir",             type=str,       default="Results")
    parser.add_argument("--plot_format",            type=str,       default=".png", choices=[".png", ".pdf"])        

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    # directory setting
    roi_root_dir = os.path.join(args.data_root, args.data_type)
    metadata_dir = os.path.join(args.data_root, args.meta_dir)
    stat_result_dir = os.path.join(args.data_root, args.result_dir, "Stats")
    if not os.path.exists(stat_result_dir):
        os.makedirs(stat_result_dir)    

    # load roi stage information
    roi_info_path = os.path.join(metadata_dir, "StudyROI_Info.xlsx")
    study_roi_df = pd.read_excel(roi_info_path)
    # Filter Tumor ROIs
    study_roi_df = study_roi_df[study_roi_df["ROI_Location"] == "Tumor"]
    roi_id_lst = study_roi_df["ROI_ID"].tolist()
    roi_diag_lst = study_roi_df["ROI_Diag"].tolist()
    # collect stage lesion list
    stage_lst = ["AAH", "AIS", "MIA", "ADC"]
    stage_lesion_dict = {}
    for cur_stage in stage_lst:
        lesion_lst = []
        for roi_id, roi_diag in zip(roi_id_lst, roi_diag_lst):
            if roi_diag == cur_stage:
                lesion_name = roi_id[:roi_id.rfind("-ROI")]
                if lesion_name not in lesion_lst:
                    lesion_lst.append(lesion_name)
        stage_lesion_dict[cur_stage] = lesion_lst
    # # print lesion statistics
    # for key in stage_lesion_dict.keys():
    #     print("{} has {} lesions.".format(key, len(stage_lesion_dict[key])))
    
    # filter ADC lesions
    adc_lesion_lst = stage_lesion_dict["ADC"]
    adc_lesion_dict = {}
    for roi_id, roi_diag in zip(roi_id_lst, roi_diag_lst):
        lesion_name = roi_id[:roi_id.rfind("-ROI")]
        if lesion_name not in adc_lesion_lst:
            continue
        if roi_diag != "ADC":
            continue
        if lesion_name not in adc_lesion_dict:
            adc_lesion_dict[lesion_name] = [roi_id, ]
        else:
            adc_lesion_dict[lesion_name].append(roi_id)
    import pdb; pdb.set_trace()
    