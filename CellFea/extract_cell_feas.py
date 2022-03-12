# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import numpy as np
from skimage import io
from scipy.io import loadmat
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
        print("Work on {}/{}".format(ind+1, len(roi_list)))
        cur_roi_mask_path = os.path.join(mask_root, cur_roi + ".npy")
        if not os.path.exists(cur_roi_mask_path):
            print("No mask")
        seg_mask = np.squeeze(tifffile.imread(cur_roi_mask_path))
        inst_list = sorted(list(np.unique(seg_mask)))
        inst_list.remove(0) # remove background
        # locate stain images
        stain_imgs = np.zeros((seg_mask.shape[0], seg_mask.shape[1], len(stain_list), np.float32)
        cur_roi_stain_dir = os.path.join(stain_root, cur_roi)
        for sind, cur_stain in enumerate(stain_list):
            cur_stain_path = os.path.join(cur_roi_stain_dir, cur_stain + ".tiff")
            if not os.path.exists(cur_stain_path):
                print("No stain {}".format(cur_stain))
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
            cx = int(cnt_m["m10"] / cnt_m["m00"])
            cy = int(cnt_m["m01"] / cnt_m["m00"])
            scell_fea.extend([cx, cy])
            cell_stain_pixels = stain_imgs[inst_map==1]
            cell_stain_means = np.mean(cell_stain_pixels, axis=0)
            scell_fea.extend(cell_stain_means.tolist())
            roi_cell_feas.append(scell_fea)
        # save roi cell features
        np.save(os.path.join(cellfea_dir, cur_roi + ".npy"), np.asarray(roi_cell_feas))

        #     inst_map = np.array(seg_mask == inst_id, np.uint8)
        #     y1, y2, x1, x2  = fea_utils.bounding_box(inst_map)
        #     y1 = y1 - 2 if y1 - 2 >= 0 else y1
        #     x1 = x1 - 2 if x1 - 2 >= 0 else x1
        #     x2 = x2 + 2 if x2 + 2 <= seg_mask.shape[1] - 1 else x2
        #     y2 = y2 + 2 if y2 + 2 <= seg_mask.shape[0] - 1 else y2
        #     inst_map_crop = inst_map[y1:y2, x1:x2]
        #     contours, hierarchy = cv2.findContours(inst_map_crop, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        #     contours = sorted(contours, key=lambda x: cv2.contourArea(x), reverse=True)
        #     cell_cnt = contours[0]
        #     # obtain central coordinates and area
        #     cnt_area = float(cv2.contourArea(cell_cnt) * 1.0)
        #     if cnt_area < args.min_cell_area:
        #         continue
        #     cnt_m = cv2.moments(cell_cnt)
        #     cx = int(cnt_m["m10"] / cnt_m["m00"])
        #     cy = int(cnt_m["m01"] / cnt_m["m00"])
        #     # combine features
        #     cnt_x, cnt_y, cnt_w, cnt_h = cv2.boundingRect(cell_cnt)
        #     rect_area = float(cnt_w * cnt_h * 1.0)
        #     cnt_extent = cnt_area / rect_area
        #     cnt_aspect_ratio = float(cnt_w) / float(cnt_h)
        #     scell_fea = [cx, cy]
        #     ## stain features
        #     cell_mask = np.zeros((y2-y1, x2-x1), np.uint8)
        #     cv2.drawContours(cell_mask, contours=contours, contourIdx=0, color=1, thickness=-1)
        #     cell_stains = stain_imgs[y1:y2, x1:x2]
        #     cell_stain_pixels = cell_stains[cell_mask==1]
        #     cell_stain_means = np.mean(cell_stain_pixels, axis=0)
        #     scell_fea.extend(cell_stain_means.tolist())
        #     roi_cell_feas.append(scell_fea)
        #
        # # save roi cell features
        # np.save(os.path.join(cellfea_dir, cur_roi + ".npy"), np.asarray(roi_cell_feas))
