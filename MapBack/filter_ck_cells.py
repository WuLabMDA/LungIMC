# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse, pickle
import numpy as np
import pyreadr



def set_args():
    parser = argparse.ArgumentParser(description = "IMC Cell Filtering")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--batchcorrection_dir",    type=str,       default="BatchCorrection")
    parser.add_argument("--vis_dir",                type=str,       default="Phenotype")
    parser.add_argument("--fea_option",             type=str,       default="Transform", choices = ["Transform", "SelfCorrect", "ControlCorrect"])
    parser.add_argument("--seed",                   type=int,       default=1234)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    key_antibody = "ck"
    key_cluster_ids = None # starting from 1
    # Transform
    if args.fea_option == "Transform":
        key_cluster_ids = [2, 1, 9, 17, 33, 42, 43, 27, 35, 18, 26, 20, 10, 19]
    elif args.fea_option == "SelfCorrect":
        key_cluster_ids = [57, 58, 59, 60, 50, 55]
    elif args.fea_option == "ControlCorrect":
        key_cluster_ids = [64, 62, 63, 39, 40, 48, 47, 56]
    else:
        print("No option {}".format(args.fea_option))
        exit()

    # cellfea root
    cellfea_dir = os.path.join(args.data_root, args.batchcorrection_dir, "RData")

    # load communites
    community_path = os.path.join(cellfea_dir, args.fea_option + "CommunitiesSOM.RData")
    community_rdata = pyreadr.read_r(community_path)
    communities = community_rdata["communities"]
    communities = communities.to_numpy().astype(int)
    communities = np.squeeze(communities).tolist()
    # cell id starting from 1
    key_cell_ids = [cell_ind + 1 for (cell_ind, community_id) in enumerate(communities) if community_id in key_cluster_ids]
    # load roi files and cell numbers
    metadata_path = os.path.join(cellfea_dir, "StudyMeta.RData")
    metadata_dict = pyreadr.read_r(metadata_path)
    fea_filenames = metadata_dict["fea_filenames"]["fea_filenames"].to_list()
    roi_nrows = metadata_dict["roi_nrows"]["roi_nrows"].to_list()

    # build roi vs cell ids mapping
    roi_cells_dict = {}
    cum_num = 0
    for roi_name, cell_num in zip(fea_filenames, roi_nrows):
        start_num = cum_num
        end_num = cum_num + cell_num
        cell_list = [cell_id - start_num for cell_id in key_cell_ids if cell_id > start_num and cell_id <= end_num]
        cum_num = cum_num + cell_num
        roi_cells_dict[roi_name] = cell_list

    # save roi vs cell ids mapping
    key_roi_cells_path = os.path.join(args.data_root, args.vis_dir, args.fea_option, "{}_roi_cells.pkl".format(key_antibody))
    with open(key_roi_cells_path, 'wb') as handle:
        pickle.dump(roi_cells_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
