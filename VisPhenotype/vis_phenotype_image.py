# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import pandas as pd
import numpy as np
from skimage import io
import tifffile


def set_args():
    parser = argparse.ArgumentParser(description = "Check Cell Phenotype on Raw Image")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--denoise_dir",            type=str,       default="Denoise")
    parser.add_argument("--steinbock_dir",          type=str,       default="SteinbockAll")
    parser.add_argument("--phenotype_dir",          type=str,       default="PhenotypeAll")
    parser.add_argument("--phenotype_id",           type=int,       default=3)


    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = set_args()

    antibody_list = ["CD8a", "CD3e", "CK"]
    denoised_img_dir = os.path.join(args.data_root, args.data_type, args.denoise_dir, "DenoisedROIs")

    steinbock_dir = os.path.join(args.data_root, args.data_type, args.steinbock_dir)
    # segmentation dir
    cell_seg_dir = os.path.join(steinbock_dir, "masks_deepcell")

    # load cell id & phenotype information
    cell_id_phenotype_path = os.path.join(steinbock_dir, "cell_id_phenotype.csv")
    cell_id_phenotype_df = pd.read_csv("./cell_id_phenotype.csv", index_col=None)
    cell_ids = cell_id_phenotype_df["CellID"].tolist()
    cell_phenotypes = cell_id_phenotype_df["CellPhenotype"].tolist()
    interested_ids = [cell_ids[ind] for ind in np.arange(len(cell_ids)) if cell_phenotypes[ind] == args.phenotype_id]
    import pdb; pdb.set_trace()


