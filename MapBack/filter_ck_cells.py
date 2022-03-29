# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse, pickle
import numpy as np
import pyreadr



def set_args():
    parser = argparse.ArgumentParser(description = "IMC Cell Filtering")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--batchcorrection_dir",    type=str,       default="BatchCorrection")
    parser.add_argument("--phenotype_dir",          type=str,       default="Phenotype")
    parser.add_argument("--fea_option",             type=str,       default="Transform", choices = ["Transform", "SelfCorrect", "ControlCorrect"])
    parser.add_argument("--seed",                   type=int,       default=1234)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    ck_cluster_ids = None # starting from 1
    # Transform
    if args.fea_option == "Transform":
        ck_cluster_ids = [7, 8, 15, 16, 24, 32]

    # cellfea root
    cellfea_dir = os.path.join(args.data_root, args.batchcorrection_dir, "RData")

    # load communites
    community_path = os.path.join(cellfea_dir, args.fea_option + "CommunitiesSOM.RData")
    community_rdata = pyreadr.read_r(community_path)
    communities = community_rdata["communities"]
    communities = communities.to_numpy().astype(int)
    communities = np.squeeze(communities).tolist()
    import pdb; pdb.set_trace()
    # cell id starting from 1
    ck_cell_ids = [cell_ind + 1 for (cell_ind, community_id) in enumerate(communities) if community_id in ck_cluster_ids]
    import pdb; pdb.set_trace()

    # load roi files and cell numbers
    metadata_path = os.path.join(cellfea_dir, "StudyMeta.RData")
    metadata_dict = pyreadr.read_r(metadata_path)
    fea_filenames = metadata_dict["fea_filenames"]["fea_filenames"].to_list()
    roi_nrows = metadata_dict["roi_nrows"]["roi_nrows"].to_list()
