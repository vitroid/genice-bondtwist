# coding: utf-8
"""
Bond twist

Matsumoto, M., Yagasaki, T. & Tanaka, H. Chiral Ordering in Supercooled Liquid Water and Amorphous Ice. Phys. Rev. Lett. 115, 197801 (2015).
"""

from math import atan2, sin, cos, pi
import cmath
import colorsys
import logging
import io
import sys

import numpy as np
import networkx as nx
import yaplotlib as yp

from genice_svg.formats.svg import Render, draw_cell
from png import Render as pRender

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
                    angle = atan2(sine,cosine)# * 180 / pi
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
        s = "# i, j, center, of, bond, bond-twist\n"
        s += "{0}\n".format(tag)
        s += self.relcoord.shape[0].__str__()+"\n"
        for a,b,center,twist in self.iter():
            s += "{0} {1} {2:.4f} {3:.4f} {4:.4f} {5:.4f} {6:.4f}\n".format(a,b,*center*10, twist.real, twist.imag)
        s += "-1 -1 0 0 0 0 0\n"
        return s

    def yaplot(self):
        maxcir = 16
        maxrad = 6
        # 16x6=96 color variations
        s = ""
        for cir in range(maxcir):
            for rad in range(maxrad):
                angle = cir*360/maxcir
                hue = ((angle - 60 + 360) / 360) % 1.0
                bri = rad / (maxrad-1)
                sat = sin(angle*pi/180)**2
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
            angle = atan2(sine,cosine) * 180 / pi
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


    def svg(self, rotmat, render=Render):
        Rsphere = 0.03  # nm
        Rcyl    = 0.02  # nm
        RR      = (Rsphere**2 - Rcyl**2)**0.5
        prims = []
        proj = np.dot(self.cell, rotmat)
        xmin, xmax, ymin, ymax = draw_cell(prims, proj)
        for a,b in self.graph.edges():
            twist = self._bondtwist((a,b))
            if twist == 0:
                # not an appropriate pair
                continue
            d = self.relcoord[b] - self.relcoord[a]
            d -= np.floor(d + 0.5)
            apos = np.dot(self.relcoord[a], proj)
            dp = np.dot(d, proj)
            bpos = apos + dp
            center = apos + dp/2
            o = dp / np.linalg.norm(dp)
            o *= RR
            
            #Color setting
            cosine = twist.real
            sine   = twist.imag
            angle = atan2(sine,cosine) * 180 / pi
            rad = sine**2 + cosine**2
            hue = ((angle - 60 + 360) / 360) % 1.0
            bri = rad
            sat = sin(angle*pi/180)**2
            r,g,b = colorsys.hsv_to_rgb(hue, sat, bri)
            colorcode = "#{0:02x}{1:02x}{2:02x}".format(int(r*255),int(g*255),int(b*255))
            prims.append([center, "L", apos+o, bpos-o,Rcyl, {"fill":colorcode}])
        for v in self.relcoord:
            prims.append([np.dot(v, proj),"C",Rsphere, {"fill":"#fff"}]) #circle
        return render(prims, Rsphere, shadow=False,
                   topleft=np.array((xmin,ymin)),
                   size=(xmax-xmin, ymax-ymin))



    def svg2(self, rotmat, phasefiles, render=Render):
        #
        # Twist order parameter to distinguish phases.
        #
        logger = logging.getLogger()

        pxF = dict()
        phases = list(phasefiles)
        logger.info(phases)
        X = []
        for phase in phases:
            pxF[phase] = np.load(open(phasefiles[phase], "rb"))
            X.append(pxF[phase].reshape([1600,]))
        X = np.array(X).T

        chirs = []
        chiis = []
        for a,b in self.graph.edges():
            twist = self._bondtwist((a,b))
            if twist == 0:
                continue
            chirs.append(twist.real)
            chiis.append(twist.imag)

        #pX: p(chi)
        pX = np.histogram2d(chirs, chiis,
                            bins=(40,40),
                            range=[[-1.,1.],[-1.,1.]],
                            normed=True)
        Y = pX[0].reshape([1600,])
        from sklearn import linear_model
        reg = linear_model.LinearRegression(fit_intercept=False)
        reg.fit(X,Y)

        # pF : p(F)

        pF = dict()
        for i, phase in enumerate(phases):
            pF[phase] = reg.coef_[i]

        pFx = dict()
        for phase in phases:
            pFx[phase] = pxF[phase]*pF[phase] / pX[0]
        pFx["ice"] = pFx["1c"] + pFx["1h"]

        Rsphere = 0.03  # nm
        Rcyl    = 0.02  # nm
        RR      = (Rsphere**2 - Rcyl**2)**0.5
        prims = []
        proj = np.dot(self.cell, rotmat)
        xmin, xmax, ymin, ymax = draw_cell(prims, proj)
        for a,b in self.graph.edges():
            twist = self._bondtwist((a,b))
            if twist == 0:
                # not an appropriate pair
                continue
            d = self.relcoord[b] - self.relcoord[a]
            d -= np.floor(d + 0.5)
            apos = np.dot(self.relcoord[a], proj)
            dp = np.dot(d, proj)
            bpos = apos + dp
            center = apos + dp/2
            o = dp / np.linalg.norm(dp)
            o *= RR
            
            #Color setting
            bin = int(twist.real*20+20), int(twist.imag*20+20)
            green = pFx["HDL"][bin]
            blue  = pFx["LDL"][bin]
            red   = pFx["ice"][bin]
            logger.debug((red,green,blue))
            if green < 0:
                green = 0
            if green > 1:
                green = 1
            if red < 0:
                red = 0
            if red > 1:
                red = 1
            if blue < 0:
                blue = 0
            if blue > 1:
                blue = 1
            colorcode = "#{0:02x}{1:02x}{2:02x}".format(int(red*255),int(green*255),int(blue*255))
            prims.append([center, "L", apos+o, bpos-o,Rcyl, {"fill":colorcode}])
        for v in self.relcoord:
            prims.append([np.dot(v, proj),"C",Rsphere, {"fill":"#fff"}]) #circle
        return render(prims, Rsphere, shadow=False,
                   topleft=np.array((xmin,ymin)),
                   size=(xmax-xmin, ymax-ymin))


    def svg3(self, rotmat, render=Render):
        #
        # Twist order parameter to distinguish phases.
        # Simple criteria.
        #
        logger = logging.getLogger()
        anglerange = 20
        
        Rsphere = 0.03  # nm
        Rcyl    = 0.02  # nm
        RR      = (Rsphere**2 - Rcyl**2)**0.5
        prims = []
        proj = np.dot(self.cell, rotmat)
        xmin, xmax, ymin, ymax = draw_cell(prims, proj)
        for a,b in self.graph.edges():
            twist = self._bondtwist((a,b))
            if twist == 0:
                # not an appropriate pair
                continue
            d = self.relcoord[b] - self.relcoord[a]
            d -= np.floor(d + 0.5)
            apos = np.dot(self.relcoord[a], proj)
            dp = np.dot(d, proj)
            bpos = apos + dp
            center = apos + dp/2
            o = dp / np.linalg.norm(dp)
            o *= RR
            
            #Color setting
            cosine, sine  = twist.real, twist.imag
            angle = atan2(sine,cosine) * 180 / pi
            radius = (cosine**2 + sine**2)**0.5
            red = 0
            green = 0
            blue = 0
            if radius > 0.7:
                if (-anglerange < angle <anglerange ) or (angle < -180+anglerange) or (180-anglerange < angle):
                    red = 1
                else:
                    blue = 1
            else:
                green = 1
            colorcode = "#{0:02x}{1:02x}{2:02x}".format(int(red*255),int(green*255),int(blue*255))
            prims.append([center, "L", apos+o, bpos-o,Rcyl, {"fill":colorcode}])
        for v in self.relcoord:
            prims.append([np.dot(v, proj),"C",Rsphere, {"fill":"#fff"}]) #circle
        return render(prims, Rsphere, shadow=False,
                   topleft=np.array((xmin,ymin)),
                   size=(xmax-xmin, ymax-ymin))


    
    def svg4(self, rotmat, phasefiles, render=Render):
        #
        # Twist order parameter to distinguish phases.
        # Color only when likelihood > 0.5.
        #
        logger = logging.getLogger()
        logger.info("Calculating OP...")
        pxF = dict()
        phases = list(phasefiles)
        logger.info(phases)
        X = []
        for phase in phases:
            pxF[phase] = np.load(open(phasefiles[phase], "rb"))
            X.append(pxF[phase].reshape([1600,]))
        X = np.array(X).T

        chirs = []
        chiis = []
        for a,b in self.graph.edges():
            twist = self._bondtwist((a,b))
            if twist == 0:
                continue
            chirs.append(twist.real)
            chiis.append(twist.imag)

        #pX: p(chi)
        pX = np.histogram2d(chirs, chiis,
                            bins=(40,40),
                            range=[[-1.,1.],[-1.,1.]],
                            normed=True)
        Y = pX[0].reshape([1600,])
        from sklearn import linear_model
        reg = linear_model.LinearRegression(fit_intercept=False)
        reg.fit(X,Y)

        # pF : p(F)

        pF = dict()
        for i, phase in enumerate(phases):
            pF[phase] = reg.coef_[i]

        pFx = dict()
        for phase in phases:
            pFx[phase] = pxF[phase]*pF[phase] / pX[0]
        pFx["ice"] = pFx["1c"] + pFx["1h"]
        logger.info("Done.")

        Rsphere = 0.03  # nm
        Rcyl    = 0.02  # nm
        RR      = (Rsphere**2 - Rcyl**2)**0.5
        prims = []
        proj = np.dot(self.cell, rotmat)
        xmin, xmax, ymin, ymax = draw_cell(prims, proj)
        for a,b in self.graph.edges():
            twist = self._bondtwist((a,b))
            if twist == 0:
                # not an appropriate pair
                continue
            d = self.relcoord[b] - self.relcoord[a]
            d -= np.floor(d + 0.5)
            apos = np.dot(self.relcoord[a], proj)
            dp = np.dot(d, proj)
            bpos = apos + dp
            center = apos + dp/2
            o = dp / np.linalg.norm(dp)
            o *= RR
            
            #Color setting
            bin = int(twist.real*20+20), int(twist.imag*20+20)
            green = pFx["HDL"][bin]
            blue  = pFx["LDL"][bin]
            red   = pFx["ice"][bin]
            
            if green > 0.5:
                green = 1.0
                blue = 0
                red = 0
            elif blue > 0.5:
                blue = 1
                red = 0.0
                green = 0.3
            elif red > 0.5:
                red = 1
                blue = 0
                green = 0
            else:
                green = 0
                blue = 0
                red = 0
            colorcode = "#{0:02x}{1:02x}{2:02x}".format(int(red*255),int(green*255),int(blue*255))
            prims.append([center, "L", apos+o, bpos-o,Rcyl, {"fill":colorcode}])
        for v in self.relcoord:
            prims.append([np.dot(v, proj),"C",Rsphere, {"fill":"#fff"}]) #circle
        return render(prims, Rsphere, shadow=False,
                   topleft=np.array((xmin,ymin)),
                   size=(xmax-xmin, ymax-ymin))
    

def hook2(lattice):
    lattice.logger.info("Hook2: Complex bond twists.")
    cell = lattice.repcell 

    positions = lattice.reppositions
    graph = nx.Graph(lattice.graph) #undirected

    twist = BondTwist(graph, positions, cell.mat)
    if lattice.bt_mode == "yaplot":
        print(twist.yaplot())
    elif lattice.bt_mode == "svg":
        print(twist.svg(lattice.bt_rotmat))
    elif lattice.bt_mode == "svg2":
        print(twist.svg2(lattice.bt_rotmat, lattice.bt_phases))
    elif lattice.bt_mode == "svg3":
        print(twist.svg3(lattice.bt_rotmat))
    elif lattice.bt_mode == "svg4":
        print(twist.svg4(lattice.bt_rotmat, lattice.bt_phases))
    elif "png" in lattice.bt_mode:
        if lattice.bt_mode == "png":
            image = twist.svg(lattice.bt_rotmat, render=pRender)
        elif lattice.bt_mode == "png2":
            image = twist.svg2(lattice.bt_rotmat, lattice.bt_phases, render=pRender)
        elif lattice.bt_mode == "png3":
            image = twist.svg3(lattice.bt_rotmat, render=pRender)
        elif lattice.bt_mode == "png4":
            image = twist.svg4(lattice.bt_rotmat, lattice.bt_phases, render=pRender)
        imgByteArr = io.BytesIO()
        image.save(imgByteArr, format='PNG')
        imgByteArr = imgByteArr.getvalue()
        sys.stdout.buffer.write(imgByteArr)
    else: # "file"
        print(cell.serialize_BOX9(), end="")
        print(twist.serialize("@BTWC"), end="")
    lattice.logger.info("Hook2: end.")



def hook0(lattice, arg):
    lattice.logger.info("Hook0: ArgParser.")
    lattice.bt_mode = "file"
    lattice.bt_rotmat = np.array([[1., 0, 0], [0, 1, 0], [0, 0, 1]])
    lattice.bt_phases = {"1h": "1h.4P2K.250K.pkl",
                         "1c": "1c.4P2K.250K.pkl",
                         "LDL": "LDL.ST2.240K0.88.pkl",
                         "HDL": "HDL.ST2.240K1.04.pkl"}

    if arg == "":
        pass
        #This is default.  No reshaping applied.
    else:
        args = arg.split(":")
        for a in args:
            if a.find("=") >= 0:
                key, value = a.split("=")
                lattice.logger.info("Option with arguments: {0} := {1}".format(key,value))
                if key == "rotmat":
                    value = re.search(r"\[([-0-9,.]+)\]", value).group(1)
                    lattice.bt_rotmat = np.array([float(x) for x in value.split(",")]).reshape(3,3)
                elif key == "rotatex":
                    value = float(value)*pi/180
                    cosx = cos(value)
                    sinx = sin(value)
                    R = np.array([[1, 0, 0], [0, cosx, sinx], [0,-sinx, cosx]])
                    lattice.bt_rotmat = np.dot(lattice.bt_rotmat, R)
                elif key == "rotatey":
                    value = float(value)*pi/180
                    cosx = cos(value)
                    sinx = sin(value)
                    R = np.array([[cosx, 0, -sinx], [0, 1, 0], [sinx, 0, cosx]])
                    lattice.bt_rotmat = np.dot(lattice.bt_rotmat, R)
                elif key == "rotatez":
                    value = float(value)*pi/180
                    cosx = cos(value)
                    sinx = sin(value)
                    R = np.array([[cosx, sinx, 0], [-sinx, cosx, 0], [0, 0, 1]])
                    lattice.bt_rotmat = np.dot(lattice.bt_rotmat, R)
                else:
                    #may be the file names for phases. just record them.
                    lattice.bt_phases[key] = value
                    if value == "":
                        del lattice.bt_phases[key]
            else:
                lattice.logger.info("Flags: {0}".format(a))
                if a == "yaplot":
                    lattice.bt_mode = "yaplot"
                elif a in ("svg", "svg2", "svg3", "svg4", "png", "png2", "png3", "png4"):
                    lattice.bt_mode = a
                else:
                    assert False, "Wrong options."
    lattice.logger.info("Hook0: end.")

    
hooks = {0:hook0, 2:hook2}