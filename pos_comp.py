
# This program reproduces the geometry of the cylindrical aggregates studied in the paper
# "Macroscopic coherence as an emergent property in molecular nanotubes". Parameters can
# be tuned by the user to generate all the models analyzed in the work.

import pandas as pd
import numpy as np
import math


def positions_components(N, n1, R, eps, alpha, beta, phi, h, csi):

    n2 = N // n1 # number of molecules/ring (should be integer)

    positions = []
    components = []

    for i in range(n1):
        for j in range(n2):

            theta = i*csi + j*phi

            rx = R*math.cos(theta)
            ry = R*math.sin(theta)
            rz = h*i

            xx = -math.sin(beta)*math.sin(theta + alpha*(-1)**j)
            yy = math.sin(beta)*math.cos(theta + alpha*(-1)**j)
            zz = math.cos(beta)

            positions.append([rx, ry, rz])
            components.append([xx, yy, zz])

    positions = pd.DataFrame(positions, columns = ["r$_x$", "r$_y$", "r$_z$"])
    components = pd.DataFrame(components, columns = ["$x$", "$y$", "$z$"])

    return positions, components



if __name__ == '__main__':

    n_mols = 6000 # total number of molecules in the cylinder
    n_rings = 100 # total number of rings in the cylinder
    n2 = n_mols/n_rings # number of molecules/ring

    Phi = 2 * math.pi / n2 # azimuthal angle between two consecutive dipoles in the same ring
    H = 8.3  # vertical distance between two consecutive rings, measured in Angstrom
    r = 60 # cylinder radius, measured in Angstrom
#    epsilon = 8*math.pi/45 # vertical displacement between two consecutive dipoles on neighbour rings
#    Alpha = math.pi/45 # tilt between the xy component of the dipole and the tangent to the ring
#    Beta = 11*math.pi/90 # angle between the dipole and the cylindrical axis
    epsilon, Alpha, Beta = 0, 0, 0
    xi = H*math.tan(epsilon)/r # twist angle between two molecules on consecutive rings


    mol_pos = positions_components(n_mols, n_rings, r, epsilon, Alpha, Beta, Phi, H, xi)[0].to_csv("positions.csv")
    mol_comp = positions_components(n_mols, n_rings, r, epsilon, Alpha, Beta, Phi, H, xi)[1].to_csv("components.csv")





