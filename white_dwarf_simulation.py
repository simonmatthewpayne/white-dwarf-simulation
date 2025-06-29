"""
@author: Simon Matthew Payne
"""
import numpy as np
import matplotlib.pyplot as plt
from math import *

def coupled4thorderrungekutta(differential1, differential2, x_0, y_0, z_0, d_x):
    i1 = equationsofstate(differential2, x_0, y_0, z_0)
    j1 = equationsofstate(differential1, x_0, y_0, z_0)
    i2 = equationsofstate(differential2, x_0 + d_x / 2, y_0 + d_x / 2 * i1, z_0 + d_x / 2 * j1)
    j2 = equationsofstate(differential1, x_0 + d_x / 2, y_0 + d_x / 2 * i1, z_0 + d_x / 2 * j1)
    i3 = equationsofstate(differential2, x_0 + d_x / 2, y_0 + d_x / 2 * i2, z_0 + d_x / 2 * j2)
    j3 = equationsofstate(differential1, x_0 + d_x / 2, y_0 + d_x / 2 * i2, z_0 + d_x / 2 * j2)
    i4 = equationsofstate(differential2, x_0 + d_x, y_0 + d_x * i3, z_0 + d_x * j3)
    j4 = equationsofstate(differential1, x_0 + d_x, y_0 + d_x * i3, z_0 + d_x * j3)

    y = y_0 + (d_x / 6) * (i1 + 2 * i2 + 2 * i3 + i4)
    z = z_0 + (d_x / 6) * (j1 + 2 * j2 + 2 * j3 + j4)
    x = x_0 + d_x

    return x, y, z   

def equationsofstate(differential, x, y, z):
    if y > 0:
        array = [
            y * x**2,
            -3 * z / x**2 * y**(1/3) * sqrt(1 + y**(2/3))
        ]
        return array[differential]
    else:
        return 0

def chandrasekharlimitgraph():
    Ye = 0.5  # Y_e = 0.5 for carbon star, Y_e = 0.464 for iron star
    r0 = 100**-10  # Initial radius

    Masses = []
    Radii = []
    densities = []

    for a in range(3):  # Increase 'a' for finer step sizes (slower, but more accurate)
        d_x = 1e-3 * 10**-a
        fr = []  # Final radius
        fd = []  # Final density
        fm = []  # Final mass

        for b in range(60):
            pCentral = 1e-4 * 1.5**b
            p01 = pCentral - (pCentral * r0**2) / (6 * Ye)
            m01 = (1/3) * p01 * r0**3
            r1, p1, m1 = constants(r0, p01, m01, d_x)

            fd.append(p1[-1] * 9.79e8 / Ye)
            fm.append(m1[-1] * 5.67e30 * Ye**2 * (1 / 1.98e-30))
            fr.append(r1[-1] * 7.72e6 * Ye * (1 / 6.95e-8))

        densities.append(fd)
        Masses.append(fm)
        Radii.append(fr)

        # Plotting
        plt.plot(fm, fr)
        plt.xlabel("Solar Mass")
        plt.ylabel("Solar Radius")
        plt.title(f"Chandrasekhar Limit Curve (dx = {d_x:.1e})")
        plt.grid(True)

        # Show mass-radius data
        print(np.column_stack((fm, fr)))  # Display column of Solar Mass, Solar Radius
        plt.show()

def constants(r0, p0, m0, d_x):
    r = [r0]
    p = [p0]
    m = [m0]

    while p0 > 0:   
        r0, p0, m0 = coupled4thorderrungekutta(0, 1, r0, p0, m0, d_x)
        r.append(r0)
        p.append(p0)
        m.append(m0)

    return r, p, m 

# Run the function
chandrasekharlimitgraph()

input("Press Enter to exit...")
