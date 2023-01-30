# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 20:57:22 2023

@author: Utente
"""

import numpy as np
import os

L = 6
g = 0.2
lam = 1.0
kappa = 1.0
h = 0.2
p = 0.03
a = 1
b = 1
tmax = 2 * L
tm = 0.01
steptm = 2
with open(f"jobs{L}L.sh", "w") as f:
    for i in range(10000):
        if i%2 == 0:
            seed = 100
        else:
            seed = 200
        f.write(f"nohup ./central L {L} g {g} lambda {lam} kappa {kappa} h {h} p {p} a {a} b {b} tmax {tmax} seed {seed} tm {tm} steptm {steptm}\n")