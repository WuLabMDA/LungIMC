# -*- coding: utf-8 -*-

import os, sys
import shutil, argparse, random
import json, pytz
import pandas as pd
import numpy as np
from datetime import datetime
from skimage import io
from skimage.segmentation import slic, mark_boundaries
import cv2


def set_args():
    parser = argparse.ArgumentParser(description = "Visualize Cell Neighborhood")
    parser.add_argument("--data_root",              type=str,       default="/Data")
    parser.add_argument("--data_set",               type=str,       default="HumanWholeIMC", choices=["HumanWholeIMC", "HumanSampling35"])
    parser.add_argument("--data_type",              type=str,       default="LungROIProcessing")
    parser.add_argument("--cellphenotype_dir",      type=str,       default="CellPhenotyping")    
    parser.add_argument("--steinbock_dir",          type=str,       default="Steinbock")
    parser.add_argument("--feature_dir",            type=str,       default="FeatureAnalysis")     
    parser.add_argument("--result_dir",             type=str,       default="Results")    

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = set_args()

    root_data_dir = os.path.join(args.data_root, args.data_set)
    steinbock_dir = os.path.join(root_data_dir, args.data_type, args.steinbock_dir)
    feature_ana_dir = os.path.join(root_data_dir, args.feature_dir)    
    region_prop_dir = os.path.join(steinbock_dir, "regionprops")
    cell_type_dir = os.path.join(root_data_dir, args.cellphenotype_dir, "Cell_CT_CN_Morphs")

    voronoi_vis_dir = os.path.join(root_data_dir, args.result_dir, "VisCN-Voronoi")
    if os.path.exists(voronoi_vis_dir):
        shutil.rmtree(voronoi_vis_dir)
    os.makedirs(voronoi_vis_dir)    

    roi_meta_info_path = os.path.join(steinbock_dir, "images.csv")
    roi_meta_df = pd.read_csv(roi_meta_info_path)
    roi_lst = [os.path.splitext(ele)[0] for ele in roi_meta_df["image"].tolist()]

    # define color mapping
    cell_color_dict = {"1": (141, 211, 199), "2": (228, 26, 28), "3": (190, 186, 218), "4": (251, 128, 114), "5": (128, 177, 211),
                       "6": (253, 180, 98), "7": (51, 160, 44), "8": (252, 205, 229), "9": (166, 86, 40), "10": (106, 61, 154)}    

    # generate voronoi graph one-by-one
    for roi_ind, cur_roi in enumerate(roi_lst):
        print("CN Voronoi on {:4d}/{} ROI: {}".format(roi_ind+1, len(roi_lst), cur_roi))
        roi_meta_info = roi_meta_df[roi_meta_df["image"] == cur_roi + ".tiff"]
        roi_height = roi_meta_info["height_px"].tolist()[0]
        roi_width = roi_meta_info["width_px"].tolist()[0]  
        # load cell type
        roi_cell_type_path = os.path.join(cell_type_dir, cur_roi + ".csv")
        cell_id_type_df = pd.read_csv(roi_cell_type_path)    

        seg_ids, seg_types = cell_id_type_df["seg_id"].tolist(), cell_id_type_df["cell_cn"].tolist(), 
        cell_id_type_dict = {cell_id: cell_type for cell_id, cell_type in zip(seg_ids, seg_types)}

        # obtain cell id & coordinates
        cell_morph_path = os.path.join(region_prop_dir, cur_roi + ".csv")
        cell_morph_df = pd.read_csv(cell_morph_path)
        cell_ids = cell_morph_df["Object"].tolist()

        # build voronoi graph
        cell_heights = [int(np.floor(ele + 0.5)) for ele in cell_morph_df["centroid-0"].tolist()]
        cell_widths = [int(np.floor(ele + 0.5)) for ele in cell_morph_df["centroid-1"].tolist()]
        center_dict = {}
        for ind in np.arange(len(cell_ids)):
            center_dict[(cell_heights[ind], cell_widths[ind])] = cell_id_type_dict[cell_ids[ind]]
        # create an instance of Subdiv2D
        roi_rect = (0, 0, roi_width, roi_height)
        roi_subdiv = cv2.Subdiv2D(roi_rect)
        # insert cell center
        for cw, ch in zip(cell_widths, cell_heights):
            roi_subdiv.insert((cw, ch))
        (facets, centers) = roi_subdiv.getVoronoiFacetList([])    

        roi_voronoi = np.zeros((roi_height, roi_width, 3), dtype=np.uint8)
        max_face_area = 500
        for ind in np.arange(0, len(facets)) :
            ifacet_arr = []
            for f in facets[ind] :
                ifacet_arr.append(f)
            ifacet = np.array(ifacet_arr, np.int32)
            convex_hull = cv2.convexHull(ifacet)
            face_type = None
            for key in center_dict.keys():
                c_h, c_w = key[0], key[1] 
                if cv2.pointPolygonTest(convex_hull, (c_w, c_h), False) == 1:
                    face_type = center_dict[key]
                    break
            convex_hull_area = cv2.contourArea(convex_hull)
            if convex_hull_area <= max_face_area and face_type != None:
                cv2.fillConvexPoly(roi_voronoi, ifacet, cell_color_dict[str(face_type)])   
        # save voronoi graph
        voronoi_path = os.path.join(voronoi_vis_dir, cur_roi + ".png")
        io.imsave(voronoi_path, roi_voronoi)
        # # superpixel
        # segments_slic = slic(roi_voronoi, n_segments=1000, compactness=10, sigma=1, start_label=1)
        # roi_slic = mark_boundaries(roi_voronoi, segments_slic)                         