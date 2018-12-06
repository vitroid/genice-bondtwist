# coding: utf-8
"""
Bond twist

Matsumoto, M., Yagasaki, T. & Tanaka, H. Chiral Ordering in Supercooled Liquid Water and Amorphous Ice. Phys. Rev. Lett. 115, 197801 (2015).
"""

import numpy as np
import math
import cmath
import networkx as nx

class BondTwist():
    def __init__(self, graph, relcoord, cell):
        self.graph = graph
        self.relcoord = relcoord
        self.cell = cell

    def _bondtwist(self, edge):
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
                    # print(angle*180/math.pi)
                    #sum += math.sin(angle*3)
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

    def serialize_BTWC(self):
        s = "@BTWC\n"
        s += self.relcoord.shape[0].__str__()+"\n"
        for a,b,center,twist in self.iter():
            s += "{0} {1} {2:.4f} {3:.4f} {4:.4f} {5:.4f} {6:.4f}\n".format(a,b,*center*10, twist.real, twist.imag)
        s += "-1 -1 0 0 0 0 0\n"
        return s


def hook4(lattice):
    lattice.logger.info("Hook4: Complex bond twists.")
    cell = lattice.repcell 
    print(cell.serialize_BOX9(), end="")

    positions = lattice.reppositions
    graph = nx.Graph(lattice.spacegraph) #undirected

    print("# i, j, center, of, bond, bond-twist")
    twist = BondTwist(graph, positions, cell.mat)
    print(twist.serialize_BTWC(), end="")
    lattice.logger.info("Hook4: end.")


hooks = {4:hook4}
