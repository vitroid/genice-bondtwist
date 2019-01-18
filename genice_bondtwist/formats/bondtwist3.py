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
from genice_bondtwist.formats.bondtwist  import hook0
from genice_bondtwist.formats.bondtwist2 import BondTwist2


class BondTwist3(BondTwist2):
    def _bondtwist(self, edge):
        #1. make vertices list
        logger = logging.getLogger()
        a,b = edge
        # 遠くても近い順に4つとってくる
        alist = [x for x in self.proximity[a]]
        # 遠いものは採用しない
        #alist = [x for x in self.proximity[a] if self.proximity[a][x]['distance']<3.0]
        assert len(alist) >= 4
        alist.sort(key=lambda x:self.proximity[a][x]['distance'])
        alist = alist[:4]
        #for x in alist:
        #    print(self.proximity[a][x]['distance'])
        #assert False
        aset = set(alist)

        blist = [x for x in self.proximity[b]]
        assert len(blist) >= 4
        blist.sort(key=lambda x:self.proximity[b][x]['distance'])
        blist = blist[:4]
        #for x in blist:
        #    print(self.proximity[b][x]['distance'])
        #assert False
        bset = set(blist)

        center = self.relcoord[a]
        # a,bが、互いの最近接分子の1つでなければならない。
        #これにより、水素結合の定義を全く使わず一意的に結合を定義できる。
        #ただし、気体には使えない。界面や水以外にも使えない。
        #その点では、最初の定義がいちばん使い道が広い。
        if a not in bset:
            return 0
        if b not in aset:
            return 0
        vertices = aset | bset
        #2. calculate the direction vectors
        vecs = dict()
        for vertex in vertices:
            d = self.relcoord[vertex] - center
            d -= np.floor(d+0.5)
            vecs[vertex] = np.dot(d, self.cell)
        pivot = vecs[b] - vecs[a]
        pivot /= np.linalg.norm(pivot)
        #direction vectors, perpendicular to the pivot vector
        for vertex in vertices:
            v = vecs[vertex] - np.dot(pivot, vecs[vertex])*pivot
            vecs[vertex] = v / np.linalg.norm(v)
        sum = 0.0
        n = 0
        for i in aset - set(edge):
            for j in bset - set(edge):
                if not i==j: #not triangle
                    sine = np.dot(np.cross(vecs[i],vecs[j]), pivot)
                    cosine = np.dot(vecs[i],vecs[j])
                    angle = math.atan2(sine,cosine)
                    sum += cmath.exp(angle*3j)
                    n += 1
        if n == 0:
            return 0
        op = sum / n
        return op
    

def hook1(lattice):
    lattice.logger.info("Hook1: Complex bond twists rev.3.")
    cell = lattice.repcell.mat 
    positions = lattice.reppositions
    twist = BondTwist3(positions, cell, Rcutoff=0.45)
    twist.graph = twist.proximity
    if lattice.bt_mode == "yaplot":
        print(twist.yaplot())
    elif lattice.bt_mode == "svg":
        print(twist.svg(lattice.bt_rotmat))
    else: # "file"
        print(lattice.repcell.serialize_BOX9(), end="")
        print(twist.serialize("@BTC3"), end="")
    lattice.logger.info("Hook1: end.")


    
hooks = {0:hook0, 1:hook1}
