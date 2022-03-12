# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse, pytz
import numpy as np
from skimage import io
from scipy.io import loadmat
from datetime import datetime
import tifffile, cv2
import warnings


def set_args():
    parser = argparse.ArgumentParser(description = "IMC Cell Feature Extraction")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="Study", choices = ["Study", "Tonsil"])
    parser.add_argument("--merge_dir",              type=str,       default="MergeSeg")
    parser.add_argument("--mask_dir",               type=str,       default="Mask")
    parser.add_argument("--stain_dir",              type=str,       default="Stains")
    parser.add_argument("--fea_dir",                type=str,       default="CellFeas")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()
    # prepare directory
    root_dir = os.path.join(args.data_root, args.data_type + "Processing")
    mask_root = os.path.join(root_dir, args.merge_dir, args.mask_dir)
    stain_root = os.path.join(root_dir, args.merge_dir, args.stain_dir)
    cellfea_dir = os.path.join(root_dir, args.fea_dir)
    if os.path.exists(cellfea_dir):
        shutil.rmtree(cellfea_dir)
    os.makedirs(cellfea_dir)
    # Stain list
    stain_list = ['B2M', 'B7_H3', 'CD11b', 'CD11c', 'CD14', 'CD163', 'CD19', 'CD31', 'CD33', 'CD3e',
                  'CD4', 'CD45', 'CD45RO', 'CD68', 'CD73', 'CD8a', 'CD94', 'CK', 'CTLA_4', 'FoxP3',
                  'GranzymeB', 'HLA_DR', 'ICOS', 'IDO_1', 'Ir191_193', 'Ki67', 'LAG3', 'MPO', 'NaKATPase', 'PD_1',
                  'PD_L1', 'TIGIT', 'TIM3', 'VISTA', 'aSMA']
    # deal with all rois
    roi_list = sorted([os.path.splitext(ele)[0] for ele in os.listdir(mask_root) if ele.endswith(".npy")])
    for ind, cur_roi in enumerate(roi_list):
        cur_time_str = datetime.now(pytz.timezone('America/Chicago')).strftime("%H:%M:%S")
        print("Work on {}/{} Current time: {}".format(ind+1, len(roi_list), cur_time_str))
        cur_roi_mask_path = os.path.join(mask_root, cur_roi + ".npy")
        seg_mask = np.load(cur_roi_mask_path)
        inst_list = sorted(np.unique(seg_mask).tolist())
        inst_list.remove(0) # remove background
        # collect stain images
        stain_imgs = np.zeros((seg_mask.shape[0], seg_mask.shape[1], len(stain_list)), np.float32)
        cur_roi_stain_dir = os.path.join(stain_root, cur_roi)
        for sind, cur_stain in enumerate(stain_list):
            cur_stain_path = os.path.join(cur_roi_stain_dir, cur_stain + ".tiff")
            stain_imgs[:, :, sind] = tifffile.imread(cur_stain_path).astype(np.float32)
        # extract cell features
        roi_cell_feas = []
        for inst_id in inst_list:
            scell_fea = []
            # locate the cell mask region
            inst_map = np.array(seg_mask == inst_id, np.uint8)
            contours, hierarchy = cv2.findContours(inst_map, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
            cell_cnt = contours[0]
            # extract cell centers
            cnt_m = cv2.moments(cell_cnt)
            if cnt_m["m00"] != 0:
                cx = int(cnt_m["m10"] / cnt_m["m00"])
                cy = int(cnt_m["m01"] / cnt_m["m00"])
            else:
                cx, cy = 0, 0
            scell_fea.extend([cx, cy])
            cell_stain_pixels = stain_imgs[inst_map==1]
            cell_stain_means = np.mean(cell_stain_pixels, axis=0)
            scell_fea.extend(cell_stain_means.tolist())
            roi_cell_feas.append(scell_fea)
        # save roi cell features
        np.save(os.path.join(cellfea_dir, cur_roi + ".npy"), np.asarray(roi_cell_feas))
