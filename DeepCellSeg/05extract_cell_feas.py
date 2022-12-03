# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse, pytz
import pickle
import numpy as np
from skimage import io
from scipy.io import loadmat
from datetime import datetime
import tifffile, cv2
import warnings


def set_args():
    parser = argparse.ArgumentParser(description = "IMC Cell Feature Extraction")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--cellseg_dir",            type=str,       default="CellSeg")        
    parser.add_argument("--merge_dir",              type=str,       default="MergeCellSeg")
    parser.add_argument("--denoise_dir",            type=str,       default="Denoise")
    parser.add_argument("--steinbock_dir",          type=str,       default="Steinbock")
    parser.add_argument("--fea_dir",                type=str,       default="CellFeas")
    parser.add_argument("--cnt_dir",                type=str,       default="Contour")
    parser.add_argument("--tif_dir",                type=str,       default="masks_deepcell")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()
    # prepare directory
    root_dir = os.path.join(args.data_root, args.data_type)
    mask_root = os.path.join(root_dir, args.cellseg_dir, args.merge_dir, "Mask")
    stain_root = os.path.join(root_dir, args.denoise_dir, "DenoisedROIs")
    # cell feature directory
    cellfea_dir = os.path.join(root_dir, args.cellseg_dir, args.merge_dir, args.fea_dir)
    if os.path.exists(cellfea_dir):
        shutil.rmtree(cellfea_dir)
    os.makedirs(cellfea_dir)
    # cell contour directory
    cnt_dir_root = os.path.join(root_dir, args.cellseg_dir, args.merge_dir, args.cnt_dir)
    if os.path.exists(cnt_dir_root):
        shutil.rmtree(cnt_dir_root)
    os.makedirs(cnt_dir_root)
    # cell tif mask directory
    tifmask_root = os.path.join(root_dir, args.steinbock_dir, args.tif_dir)
    if os.path.exists(tifmask_root):
        shutil.rmtree(tifmask_root)
    os.makedirs(tifmask_root)    

    # Stain list
    stain_list = ['B2M', 'B7_H3', 'CD11b', 'CD11c', 'CD14', 'CD163', 'CD19', 'CD31', 'CD33', 'CD3e',
                  'CD4', 'CD45', 'CD45RO', 'CD68', 'CD73', 'CD8a', 'CD94', 'CK', 'CTLA_4', 'FoxP3',
                  'GranzymeB', 'HLA_DR', 'ICOS', 'IDO_1', 'Ir191', 'Ki67', 'LAG3', 'MPO', 'NaKATPase', 'PD_1',
                  'PD_L1', 'TIGIT', 'TIM3', 'VISTA', 'aSMA']

    # deal with all rois
    roi_list = sorted([os.path.splitext(ele)[0] for ele in os.listdir(mask_root) if ele.endswith(".npy")])
    for ind, cur_roi in enumerate(roi_list):
        cur_time_str = datetime.now(pytz.timezone('America/Chicago')).strftime("%H:%M:%S")
        print("Work on {}/{} Current time: {}".format(ind+1, len(roi_list), cur_time_str))
        cur_roi_mask_path = os.path.join(mask_root, cur_roi + ".npy")
        seg_mask = np.load(cur_roi_mask_path)
        mask_h, mask_w = seg_mask.shape
        inst_list = sorted(np.unique(seg_mask).tolist())
        inst_list.remove(0) # remove background
        # collect stain images
        stain_imgs = np.zeros((seg_mask.shape[0], seg_mask.shape[1], len(stain_list)), np.float32)
        cur_roi_stain_dir = os.path.join(stain_root, cur_roi)
        for sind, cur_stain in enumerate(stain_list):
            cur_stain_path = os.path.join(cur_roi_stain_dir, cur_stain + ".tiff")
            stain_imgs[:, :, sind] = io.imread(cur_stain_path, plugin="tifffile").astype(np.float32)
        # extract cell features
        cell_idx = 0
        roi_cell_feas = []
        roi_cnt_dict = {}
        filter_mask = np.zeros((mask_h, mask_w), dtype=np.uint16)
        cell_cnt_dict = {}
        for inst_id in inst_list:
            scell_fea = []
            # locate the cell mask region
            inst_map = np.array(seg_mask == inst_id, np.uint8)
            contours, hierarchy = cv2.findContours(inst_map, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
            cell_cnt = contours[0]
            # extract cell centers
            cnt_m = cv2.moments(cell_cnt)
            if cnt_m["m00"] == 0:
                continue
            cx = int(cnt_m["m10"] / cnt_m["m00"])
            cy = int(cnt_m["m01"] / cnt_m["m00"])            
            cnt_area = cv2.contourArea(cell_cnt)
            if cnt_area < 10.0 or cnt_area > 160.0:                
                continue
            # cell feature
            scell_fea.extend([cx, cy, cnt_area])
            cell_stain_pixels = stain_imgs[inst_map==1]
            cell_stain_means = np.mean(cell_stain_pixels, axis=0)
            scell_fea.extend(cell_stain_means.tolist())
            roi_cell_feas.append(scell_fea)
            # cell contour
            cell_idx += 1
            cell_cnt_dict[cell_idx] = cell_cnt
            # cell mask
            cv2.drawContours(filter_mask, contours=[cell_cnt, ], contourIdx=0, color=cell_idx, thickness=-1)
        # save roi cell features
        roi_cell_feas = np.asarray(roi_cell_feas)
        np.save(os.path.join(cellfea_dir, cur_roi + ".npy"), roi_cell_feas)
        # save to pkl
        roi_cnt_dict["cell_cnt"] = cell_cnt_dict
        roi_cnt_dict["height"] = mask_h
        roi_cnt_dict["width"] = mask_w
        roi_cnt_dict["cell_num"] = len(cell_cnt_dict)
        roi_cnt_path = os.path.join(cnt_dir_root, cur_roi + ".pkl")
        pickle_out = open(roi_cnt_path, "wb")
        pickle.dump(roi_cnt_dict, pickle_out)
        pickle_out.close()
        # save to tiff
        tif_path = os.path.join(tifmask_root, cur_roi + ".tiff")
        io.imsave(tif_path, filter_mask, plugin="tifffile")