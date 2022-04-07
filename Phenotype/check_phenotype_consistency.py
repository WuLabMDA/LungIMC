# -*- coding: utf-8 -*-

import os, sys
import numpy as np
import argparse, shutil, pickle, pytz
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
import pyreadr


def set_args():
    parser = argparse.ArgumentParser(description = "IMC Community Detection Consistency Evaluation")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--batchcorrection_dir",    type=str,       default="BatchCorrection")
    parser.add_argument("--phenotype_dir",          type=str,       default="Phenotype")
    parser.add_argument("--seed",                   type=int,       default=1234)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    # prepare directory
    phenotype_dir = os.path.join(args.data_root, args.phenotype_dir)
    if not os.path.exists(phenotype_dir):
        os.makedirs(phenotype_dir)
    cellfea_dir = os.path.join(args.data_root, args.batchcorrection_dir, "RData")
    # load the communites from uncorrected data
    raw_community_path = os.path.join(cellfea_dir, "TransformCommunitiesSOM.RData")
    raw_community_rdata = pyreadr.read_r(raw_community_path)
    raw_communities = raw_community_rdata["communities"]
    raw_communities = raw_communities.to_numpy().astype(int)
    raw_communities = np.squeeze(raw_communities)
    # load the communties from corrected data
    correct_community_path = os.path.join(cellfea_dir, "SelfCorrectCommunitiesSOM.RData")
    correct_community_rdata = pyreadr.read_r(correct_community_path)
    correct_communities = correct_community_rdata["communities"]
    correct_communities = correct_communities.to_numpy().astype(int)
    correct_communities = np.squeeze(correct_communities)
    # draw the confusion matrix
    consistency_cm = confusion_matrix(raw_communities, correct_communities)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.matshow(consistency_cm)
    fig.colorbar(cax)
    consistency_path = os.path.join(phenotype_dir, "RawCorrectConsistency.png")
    plt.savefig(consistency_path, dpi=300)
    plt.close()
