# -*- coding: utf-8 -*-

import os, sys
import argparse, shutil
from skimage import io
import numpy as np

from deepcell.utils.plot_utils import create_rgb_image


def set_args():
    parser = argparse.ArgumentParser(description = "IMC Cell Segmentation")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--nuc_stain",              type=str,       default="191Ir")
    parser.add_argument("--mem_stain",              type=str,       default="NaKATPase")
    parser.add_argument("--data_type",              type=str,       default="Study", choices = ["Study", "Tonsil"])
    parser.add_argument("--roi_dir",                type=str,       default="SegROI")
    parser.add_argument("--vis_dir",                type=str,       default="VisROI")
    parser.add_argument("--result_dir",             type=str,       default="SegResults")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()
    seg_root_dir = os.path.join(args.data_root, args.data_type + "Processing")
    segroi_dir = os.path.join(seg_root_dir, args.roi_dir)
    visroi_dir = os.path.join(seg_root_dir, args.vis_dir)
    if os.path.exists(visroi_dir):
        shutil.rmtree(visroi_dir)
    os.makedirs(visroi_dir)
    segresult_dir = os.path.join(seg_root_dir, args.result_dir)
    if os.path.exists(segresult_dir):
        shutil.rmtree(segresult_dir)
    os.makedirs(segresult_dir)

    # traverse patients one-by-one
    patient_list = sorted([ele for ele in os.listdir(segroi_dir)])
    for ind, p_id in enumerate(patient_list[:3]):
        print("Segment {:3d}/{:3d} ID {}".format(ind+1, len(patient_list), p_id))
        cur_p_dir = os.path.join(segroi_dir, p_id)
        input_roi_arr_list = []
        cur_roi_list = sorted([ele for ele in os.listdir(cur_p_dir) if os.path.isdir(os.path.join(cur_p_dir, ele))])
        for cur_roi_name in cur_roi_list:
            cur_roi_dir = os.path.join(cur_p_dir, cur_roi_name)
            stain_list = [os.path.join(cur_roi_dir, ele + ".tif") for ele in [args.nuc_stain, args.mem_stain]]
            input_roi_arr_list.append(np.stack([io.imread(ele) for ele in stain_list], axis=2))
        input_roi_arrs = np.stack(input_roi_arr_list, axis=0).astype(np.float64)
        # Vis original image
        rgb_images = create_rgb_image(input_roi_arrs, channel_colors=['green', 'blue'])
        cur_roi_vis_dir = os.path.join(visroi_dir, p_id)
        os.makedirs(cur_roi_vis_dir)
        for ind, cur_roi_name in enumerate(cur_roi_list):
            cur_roi_vis_path = os.path.join(cur_roi_vis_dir, cur_roi_name + "_raw.png")
            io.imsave(cur_roi_vis_path, rgb_images[ind, ...])
