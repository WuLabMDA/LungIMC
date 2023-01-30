# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse
import json, pytz
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import scanpy as sc
import anndata
import squidpy as sq


def set_args():
    parser = argparse.ArgumentParser(description = "Load data into squidpy image container")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--steinbock_dir",          type=str,       default="Steinbock")    
    parser.add_argument("--object_name",            type=str,       default="objects.h5ad")    

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    # directory setting
    roi_root_dir = os.path.join(args.data_root, args.data_type)

    # load cell phenotype information
    cell_phenotype_path = os.path.join(roi_root_dir, "CellPhenotypes.json")
    cell_phenotype_dict = None
    with open(cell_phenotype_path) as fp:
        cell_phenotype_dict = json.load(fp)
    cell_ids = [ele for ele in cell_phenotype_dict.keys()]
    # collect roi list
    roi_lst = []
    for cell_id in cell_ids:
        roi_name = cell_id[:cell_id.rfind("_")]
        if roi_name not in roi_lst:
            roi_lst.append(roi_name)
    # collect roi cell numbers
    roi_cell_dict = {}
    for cell_id in cell_ids:
        roi_name = cell_id[:cell_id.rfind("_")]
        if roi_name not in roi_cell_dict:
            roi_cell_dict[roi_name] = 1
        else:   
            roi_cell_dict[roi_name] += 1
    # collect ordered cell_phenotypes
    cell_phenotypes = []
    for roi_name in roi_lst:
        cell_num = roi_cell_dict[roi_name]
        cell_phenotypes.extend([cell_phenotype_dict[roi_name + "_" + str(num)] for num in range(1, cell_num + 1)])

    # load AnnData and add cell type
    steinbock_dir = os.path.join(roi_root_dir, args.steinbock_dir)
    anndata_path = os.path.join(steinbock_dir, args.object_name)
    imc_adata = sc.read_h5ad(anndata_path)
    imc_adata.obs["cell_type"] = pd.Categorical(cell_phenotypes) 
    import pdb; pdb.set_trace()  
    # imc_ic = sq.im.ImageContainer.from_adata(imc_adata)
    
     

    