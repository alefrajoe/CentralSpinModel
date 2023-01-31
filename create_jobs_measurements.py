# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 20:57:22 2023

@author: Utente
"""

import numpy as np
import os

L = 6
g = 0.0
lam = 1.0
kappa = 1.0
h = 1.0
p = 1.0
a = 1
b = 1
tmax = 1 * L
tm = 0.1
seed = 500
steptm = 20
with open(f"jobs{seed}seed{L}L.sh", "w") as f:
    for i in range(10000):
        if i%2 == 0:
            seed += 100
        else:
            seed -= 100
        f.write(f"nohup ./central L {L} g {g} lambda {lam} kappa {kappa} h {h} p {p} a {a} b {b} tmax {tmax} seed {seed} tm {tm} steptm {steptm}\n")