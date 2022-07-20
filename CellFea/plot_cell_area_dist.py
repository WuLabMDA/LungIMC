# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse, pytz
import numpy as np
from skimage import io
from scipy.io import loadmat
from datetime import datetime
import tifffile, cv2
import warnings

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def set_args():
    parser = argparse.ArgumentParser(description = "IMC Cell Feature Extraction")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--fig_dir",                type=str,       default="LungPlots")
    parser.add_argument("--fea_dir",                type=str,       default="CellFeas")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()
    # prepare directory
    root_dir = os.path.join(args.data_root, args.data_type)
    cellfea_dir = os.path.join(root_dir, args.fea_dir)

    # deal with all rois
    roi_list = sorted([os.path.splitext(ele)[0] for ele in os.listdir(cellfea_dir) if ele.endswith(".npy")])
    print("There are {} ROIs in {}".format(len(roi_list), args.data_type))

    cell_area_list = []
    for ind, cur_roi in enumerate(roi_list):
        cur_roi_fea_path = os.path.join(cellfea_dir, cur_roi + ".npy")
        cur_roi_fea = np.load(cur_roi_fea_path)
        cell_area_list.extend(cur_roi_fea[:, 2].tolist()) # x, y, area
    cell_area_arr = np.asarray(cell_area_list)
    cell_min = np.percentile(cell_area_arr, q=1)
    cell_max = np.percentile(cell_area_arr, q=99)

    # the histogram of the data
    n, bins, patches = plt.hist(cell_area_list, 200, density=True, facecolor='g', alpha=0.75)
    plt.xlabel('Area')
    plt.ylabel('Num')
    plt.title('Histogram of Cell Area')
    plt.xlim(0, 200)
    # plt.ylim(0, 0.03)
    plt.grid(True)    

    plt_dir = os.path.join(args.fig_dir, args.fea_dir)
    if not os.path.exists(plt_dir):
        os.makedirs(plt_dir)
    export_pdf = PdfPages(os.path.join(plt_dir, "cell_area_dist.pdf"))
    export_pdf.savefig(dpi=300)
    export_pdf.close()