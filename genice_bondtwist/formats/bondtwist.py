# coding: utf-8
"""
Bond twist

Matsumoto, M., Yagasaki, T. & Tanaka, H. Chiral Ordering in Supercooled Liquid Water and Amorphous Ice. Phys. Rev. Lett. 115, 197801 (2015).
"""

import math
import cmath
import colorsys
import logging

import numpy as np
import networkx as nx
import yaplotlib as yp

class BondTwist():
    def __init__(self, graph, relcoord, cell):
        self.graph = graph
        self.relcoord = relcoord
        self.cell = cell

    def _bondtwist(self, edge):
        logger = logging.getLogger()
        #1. make vertices list
        a,b = edge
        aset = set([x for x in self.graph[a]])
        bset = set([x for x in self.graph[b]])
        center = self.relcoord[a]
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
            assert -1e-10 < np.dot(v,pivot) < +1e-10
            vecs[vertex] = v / np.linalg.norm(v)
        sum = 0.0
        n = 0
        for i in aset - set(edge):
            for j in bset - set(edge):
                if not i==j: #not triangle
                    # print(np.dot(vectors[i],pivot))
                    # print(np.dot(vectors[j],pivot))
                    sine = np.dot(np.cross(vecs[i],vecs[j]), pivot)
                    cosine = np.dot(vecs[i],vecs[j])
                    angle = math.atan2(sine,cosine)# * 180 / math.pi
                    degree = angle*180/math.pi
                    if 29 < degree < 31:
                        logger.debug("Angle, norm, cosine: {0} {1} {2}".format(degree, (sine**2+cosine**2)**0.5, cosine))
                        logger.debug((i,a,b,j))
                        d = self.relcoord[a] - self.relcoord[i]
                        d-= np.floor(d+0.5)
                        logger.debug(np.linalg.norm(d))
                        d = self.relcoord[a] - self.relcoord[b]
                        d-= np.floor(d+0.5)
                        logger.debug(np.linalg.norm(d))
                        d = self.relcoord[j] - self.relcoord[b]
                        d-= np.floor(d+0.5)
                        logger.debug(np.linalg.norm(d))
                        
                    sum += cmath.exp(angle*3j)
                    n += 1
        if n == 0:
            return 0
        op = sum / n
        return op

    def iter(self):
        for edge in self.graph.edges():
            twist = self._bondtwist(edge)
            a,b = edge
            d = self.relcoord[b] - self.relcoord[a]
            d -= np.floor(d + 0.5)
            center = self.relcoord[a] + d/2
            center -= np.floor(center)
            center = np.dot( center, self.cell )
            yield a,b,center,twist

    def serialize(self, tag):
        s = "# i, j, center, of, bond, bond-twist"\n
        s += "{0}\n".format(tag)
        s += self.relcoord.shape[0].__str__()+"\n"
        for a,b,center,twist in self.iter():
            s += "{0} {1} {2:.4f} {3:.4f} {4:.4f} {5:.4f} {6:.4f}\n".format(a,b,*center*10, twist.real, twist.imag)
        s += "-1 -1 0 0 0 0 0\n"
        return s

    def yaplot(self, Rcutoff=0.3):
        maxcir = 16
        maxrad = 6
        # 16x6=96 color variations
        s = ""
        for cir in range(maxcir):
            for rad in range(maxrad):
                angle = cir*360/maxcir
                hue = ((angle - 60 + 360) / 360) % 1.0
                bri = rad / (maxrad-1)
                sat = math.sin(angle*math.pi/180)**2
                logging.debug((angle,sat,bri))
                r,g,b = colorsys.hsv_to_rgb(hue, sat, bri)
                n = cir*maxrad + rad + 3
                s += yp.SetPalette(n,int(r*255),int(g*255),int(b*255))
        for a,b in self.graph.edges():
            twist = self._bondtwist((a,b))
            if twist == 0:
                # not an appropriate pair
                continue
            d = self.relcoord[b] - self.relcoord[a]
            d -= np.floor(d + 0.5)
            apos = np.dot(self.relcoord[a], self.cell)
            bpos = apos + np.dot(d, self.cell)
            cosine = twist.real
            sine   = twist.imag
            angle = math.atan2(sine,cosine) * 180 / math.pi
            if angle < 0:
                angle += 360
            cir = int(angle*maxcir/360+0.5)
            # rad is squared.
            rad = int(abs(twist)**2 * maxrad)
            if cir > maxcir-1:
                cir -= maxcir
            if rad > maxrad-1:
                rad = maxrad-1
            palette = cir*maxrad + rad + 3
            # logging.info((abs(twist),rad,cir,palette))
            s += yp.Color(palette)
            s += yp.Line(apos,bpos)
            s += "# angle {0} rad {1} cir {2} rad {3}\n".format(angle, abs(twist), hue, rad)
        return s


def hook2(lattice):
    lattice.logger.info("Hook2: Complex bond twists.")
    cell = lattice.repcell 
    print(cell.serialize_BOX9(), end="")

    positions = lattice.reppositions
    graph = nx.Graph(lattice.graph) #undirected

    print("# i, j, center, of, bond, bond-twist")
    twist = BondTwist(graph, positions, cell.mat)
    if lattice.bt_yaplot:
        print(twist.yaplot())
    else:
        print(twist.serialize("@BTWC"), end="")
    lattice.logger.info("Hook2: end.")



def hook0(lattice, arg):
    lattice.logger.info("Hook0: ArgParser.")
    lattice.bt_yaplot = False
    if arg == "":
        pass
        #This is default.  No reshaping applied.
    else:
        args = arg.split(":")
        for a in args:
            if a.find("=") >= 0:
                key, value = a.split("=")
                lattice.logger.info("Option with arguments: {0} := {1}".format(key,value))
            else:
                lattice.logger.info("Flags: {0}".format(a))
                if a == "yaplot":
                    lattice.bt_yaplot = True
                else:
                    assert False, "Wrong options."
    lattice.logger.info("Hook0: end.")

    
hooks = {0:hook0, 2:hook2}
