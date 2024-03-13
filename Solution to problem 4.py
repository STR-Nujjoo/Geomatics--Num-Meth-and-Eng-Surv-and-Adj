#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GEOMATICS ENGINEERING NUMERICAL METHODS PROBLEM 
Created on Mon Mar 11 16:23:28 2024
@author: Syed Tanweer Raza Nujjoo
Credits: Geomatics department in the University of Cape Town (UCT)
"""

# Importing Relevant libraries
import pandas as pd
import sympy as sp
import numpy as np

# Import excel spreadsheet ------------------------------------------------
df = pd.read_excel('Levelling network info.xlsx')
print()
print(df)

# Known heights -----------------------------------------------------------
FH1 = 100
FH2 = 99.729
hAp = FH1+df.deltah_1A[0]
hBp = hAp+df.deltah_AB[0]

# define symbols
h1,h2,dh1,dh2,dh3,hAo,hBo = sp.symbols('h1 h2 dh1 dh2 dh3 hAo hBo')

# Defining Design/Jacobian matrix with observation equations --------------
DH1 = hAo-h1 # observation equation 1
DH2 = hBo-hAo # observation equation 2
DH3 = h2-hBo # observation equation 3

A_mat=sp.Matrix([[DH1],[DH2],[DH3]]) # insert observation equation in a matrix
var=[hAo, hBo] # jacobian w.r.t dhA and dhB
A=A_mat.jacobian(var) # design matrix

# Defining the weight matrix ----------------------------------------------
weights = [1/df.sd_1A[0]**2,1/df.sd_AB[0]**2,1/df.sd_B2[0]**2]
P=np.diagflat(weights)

# Defining misclosure matrix ----------------------------------------------
l1 = dh1-(hAo-h1) # misclosure equation 1
l2 = dh2-(hBo-hAo) # misclosure equation 2
l3 = dh3-(h2-hBo) # misclosure equation 3

# Misclosure matrix
Mis=sp.Matrix([[l1],[l2],[l3]])
l = Mis.subs([(dh1, df.deltah_1A[0]),(hAo,hAp),(h1,FH1),(dh2,df.deltah_AB[0]),(hBo,hBp),(hAo,hAp),(dh3,df.deltah_B2[0]),(h2,FH2)])

# Solving for unknowns using parametric least squares formula
x = (A.T * P * A).inv() * (A.T * P * l)

# Adjusted heights
hA_adj = hAp + x[0] # adjusted height at A
hB_adj = hBp + x[1] # adjusted height at B
print()
print("Adjusted hA: ",hA_adj)
print("Adjusted hB: ",hB_adj)

# The end -----------------------------------------------------------------