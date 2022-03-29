# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse, pickle
import numpy as np


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
    ck_cluter_ids = None # starting from 1
    # Transform
    if args.fea_option == "Transform":
        ck_cluter_ids == [7, 8, 15, 16, 24, 32]

    # load communites
    cellfea_dir = os.path.join(args.data_root, args.batchcorrection_dir, "RData")
    community_path = os.path.join(cellfea_dir, args.fea_option + "CommunitiesSOM.RData")
    community_rdata = pyreadr.read_r(community_path)
    communities = community_rdata["communities"]
    import pdb; pdb.set_trace()
