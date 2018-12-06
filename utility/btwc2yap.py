#!/usr/bin/env python3

# show the bond twists of the solvation water

from fileloaders import gro, btwc, pdb
import yaplotlib as yap
import numpy as np
import pairlist as pl


def main():
    import sys
    groname, twiname = sys.argv[1:3]

    with open(twiname) as twifile:
        with open(groname) as grofile:
            while True:
                if groname[-4:] == ".pdb":
                    gr = pdb.PDB(grofile, "O", "H", resname="WAT")
                else:
                    # assume it is gromacs
                    gr = gro.Gromacs(grofile, "Ow", "Hw", "Mw")
                btw = btwc.load_BTWC(twifile)
                if btw is None or gr is None:
                    break
                celli = np.linalg.inv(gr.cell)
                # relative coord
                rs = np.array([np.dot(x, celli) for x in gr.solutes])
                ro = np.array([np.dot(x, celli) for x in gr.waters])
                print(len(rs))
                print(len(ro))
                grid = pl.determine_grid(gr.cell, 0.245)
                shell1 = np.vstack(pl.pairs_fine_hetero(rs, ro, 0.4, gr.cell, grid, distance=False))
                # some solute atoms are isolated from water.
                shell1s=set(shell1[:,0])
                shell1w=set(shell1[:,1])
                #show all bonds in yaplot
                s = ""
                for i in range(-10, 0):
                    s += yap.SetPalette(i+10+3, 255*(-i)//10, 0, 0)
                for i in range(11):
                    s += yap.SetPalette(i+10+3, 0,0,255*i//10)
                for b in btw:
                    i,j = b[0]
                    if i in shell1w or j in shell1w:
                        s += yap.Layer(1)
                    else:
                        s += yap.Layer(2)
                    r = abs(b[2])
                    sine = b[2].imag
                    s += yap.Size(r)
                    s += yap.Color(int(sine*10)+10+3)
                    s += yap.Circle(b[1])
                print(s)



if __name__ == "__main__":
    main()

