# -*- coding: utf-8 -*-

import os, sys
import xml.etree.ElementTree as ET
import numpy as np


def parse_imagescope_annotations(xml_path):
    xml_tree = ET.parse(xml_path)
    lesion_vertices = []
    roi_anno_dict = {}

    annotations_xml = xml_tree.findall('.//Annotation')
    anno_num = len(annotations_xml)
    for ianno in range(anno_num):
        annotation_xml = annotations_xml[ianno]
        annotation_name = annotation_xml.attrib["Name"]

        if annotation_name.startswith("ROI"):
            regions = annotation_xml.findall('.//Region')
            for idx in range(len(regions)):
                region_xml = regions[idx]
                vertices = []
                for vertex_xml in region_xml.findall('.//Vertex'):
                    attrib = vertex_xml.attrib
                    vertices.append([float(attrib['X']) + 0.5,
                                    float(attrib['Y']) + 0.5])
                vertices = np.asarray(vertices, dtype=np.int32)
                region_name = region_xml.attrib["Text"]
                roi_anno_dict[region_name] = vertices
        else:
            lesion_xml = annotation_xml.findall('.//Region')[0]
            for vertex_xml in lesion_xml.findall('.//Vertex'):
                attrib = vertex_xml.attrib
                lesion_vertices.append([float(attrib['X']) + 0.5,
                                float(attrib['Y']) + 0.5])
            lesion_vertices = np.asarray(lesion_vertices, dtype=np.int32)

    return lesion_vertices, roi_anno_dict