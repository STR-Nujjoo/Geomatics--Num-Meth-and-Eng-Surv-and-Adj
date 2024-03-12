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
print(df)

# Known heights -----------------------------------------------------------
h1 = 100
h2 = 99.729
hAp = h1+df.deltah_1A[0]
hBp = hAp+df.deltah_AB[0]

# define symbols
dh1A,dhAB,dhB2,dh1,dhA,dhB,dh2,hAo,hBo = sp.symbols('dh1A dhAB dhB2 dh1 dhA dhB dh2 hAo hBo')

# Defining Design/Jacobian matrix with observation equations --------------

A1 = -dh1+dhA # observation equation 1
A2 = -dhA+dhB # observation equation 2
A3 = -dhB+dh2 # observation equation 3

A_mat=sp.Matrix([[A1],[A2],[A3]]) # insert observation equation in a matrix
var=[dhA, dhB] # jacobian w.r.t dhA and dhB
A=A_mat.jacobian(var) # design matrix

# Defining the weight matrix ----------------------------------------------
weights = [1/df.sd_1A[0]**2,1/df.sd_AB[0]**2,1/df.sd_B2[0]**2]
P=np.diagflat(weights)

# Defining misclosure matrix ----------------------------------------------
l1 = dh1A-(hAo-h1) # misclosure equation 1
l2 = dhAB-(hBo-hAo) # misclosure equation 2
l3 = dhB2-(h2-hBo) # misclosure equation 3

# Misclosure matrix
Mis=sp.Matrix([[l1],[l2],[l3]])
l = Mis.subs([(dh1A, df.deltah_1A[0]),(hAo,hAp),(h1,h1),(dhAB,df.deltah_AB[0]),(hBo,hBp),(hAo,hAo),(dhB2,df.deltah_B2[0]),(h2,h2)])

# Solving for unknowns using parametric least squares formula

x = (A.T * P * A).inv() * (A.T * P * l)

# Adjusted heights

hA_adj = hAp + x[0] # adjusted height at A
hB_adj = hBp + x[1] # adjusted height at B

# The end -----------------------------------------------------------------




