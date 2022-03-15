# -*- coding: utf-8 -*-

import os, sys
import numpy as np
import random, colorsys


# Stain list
antibody_names = ['B2M', 'B7_H3', 'CD11b', 'CD11c', 'CD14', 'CD163', 'CD19', 'CD31', 'CD33', 'CD3e',
                  'CD4', 'CD45', 'CD45RO', 'CD68', 'CD73', 'CD8a', 'CD94', 'CK', 'CTLA_4', 'FoxP3',
                  'GranzymeB', 'HLA_DR', 'ICOS', 'IDO_1', 'Ir191_193', 'Ki67', 'LAG3', 'MPO', 'NaKATPase', 'PD_1',
                  'PD_L1', 'TIGIT', 'TIM3', 'VISTA', 'aSMA']


####
def random_colors(N, bright=True):
    """
    Generate random colors.
    To get visually distinct colors, generate them in HSV space then
    convert to RGB.
    """
    brightness = 1.0 if bright else 0.7
    hsv = [(i / N, 1, brightness) for i in range(N)]
    colors = list(map(lambda c: colorsys.hsv_to_rgb(*c), hsv))
    random.shuffle(colors)

    return colors
