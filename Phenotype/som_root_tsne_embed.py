# -*- coding: utf-8 -*-

import os, sys
import numpy as np
import argparse, shutil, pickle, pytz
from datetime import datetime
import pandas as pd
import phenograph
from sklearn.manifold import TSNE
import umap
import matplotlib.pyplot as plt
import seaborn as sns


def set_args():
    parser = argparse.ArgumentParser(description = "IMC Community Detection")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--batchcorrection_dir",    type=str,       default="BatchCorrection")
    parser.add_argument("--control_option",         type=str,       default="NoControl", choices = ["NoControl", "WithControl"])
    parser.add_argument("--fea_option",             type=str,       default="Corrected", choices = ["Transformed", "Corrected"])
    parser.add_argument("--sample",                 default=False,  action="store_true")
    parser.add_argument("--sample_cell_size",       type=int,       default=100000)
    parser.add_argument("--seed",                   type=int,       default=1234)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()
    # prepare directory
    phenotype_dir = os.path.join(args.data_root, args.batchcorrection_dir, args.control_option)
    fea_path = os.path.join(phenotype_dir, args.fea_option + "Feas.RData")
    community_path = os.path.join(phenotype_dir, args.fea_option + "CommunitiesSOM.RData")
    
