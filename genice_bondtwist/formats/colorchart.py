# coding: utf-8
"""
Bond twist II

Define twist from oxygen positions. HB is not referred.
界面では、水素結合の相手が水とは限らないので、ひねりの定義が難しい。

"""

import numpy as np
import math
import cmath
import networkx as nx
import logging
import colorsys

import pairlist as pl
import yaplotlib as yp

logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s %(message)s")


maxcir = 16
maxrad = 6
# 16x6=96 color variations
s = ""
for bin in range(maxcir):
    for rad in range(maxrad):
        angle = bin*360/maxcir
        hue = ((angle - 60 + 360) / 360) % 1.0
        bri = rad / (maxrad-1)
        sat = math.sin(angle*math.pi/180)**2
        logging.info((angle,sat,bri))
        r,g,b = colorsys.hsv_to_rgb(hue, sat, bri)
        n = bin*maxrad + rad + 3
        s += yp.SetPalette(n,int(r*255),int(g*255),int(b*255))
        s += yp.Color(n)
        s += yp.Circle([bin, rad, 0])

print(s)
