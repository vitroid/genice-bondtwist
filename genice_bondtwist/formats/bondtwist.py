# coding: utf-8
"""
Bond twist

Matsumoto, M., Yagasaki, T. & Tanaka, H. Chiral Ordering in Supercooled Liquid Water and Amorphous Ice. Phys. Rev. Lett. 115, 197801 (2015).
"""

import numpy as np
import math
import networkx as nx

def BondTwist(edge, graph, relcoord, cell):
    #1. make vertices list
    a,b = edge
    aset = set([x for x in graph[a]])
    bset = set([x for x in graph[b]])
    center = relcoord[a]
    vertices = aset | bset
    #2. calculate the direction vectors
    vecs = dict()
    for vertex in vertices:
        d = relcoord[vertex] - center
        d -= np.floor(d+0.5)
        vecs[vertex] = np.dot(d, cell)
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
                sum += math.sin(angle*3)
                n += 1
    if n == 0:
        return 0.0
    op = sum / n
    return op



def hook4(lattice):
    lattice.logger.info("Hook4: Ring test.")
    cellmat = lattice.repcell.mat * 10 # in AA
    positions = lattice.reppositions
    graph = nx.Graph(lattice.spacegraph) #undirected
    print("@BOX9")
    for d in range(3):
        print("{0:.4f} {1:.4f} {2:.4f}".format(*cellmat[d]))
    print("# i, j, center, of, bond, bond-twist")
    print("@BTWI")
    print(positions.shape[0])
    for edge in graph.edges():
        twist = BondTwist(edge, graph, positions, cellmat)
        a,b = edge
        d = positions[b] - positions[a]
        d -= np.floor(d + 0.5)
        center = np.dot( positions[a] + d/2, cellmat )
        print("{0} {1} {2:.4f} {3:.4f} {4:.4f} {5:.4f}".format(a,b,*center, twist))
    print("-1 -1 0 0 0 0")
    lattice.logger.info("Hook4: end.")


hooks = {4:hook4}
