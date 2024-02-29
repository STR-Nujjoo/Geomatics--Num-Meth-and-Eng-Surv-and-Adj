#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GEOMATICS ENGINERRING NUMERICAL METHODS PROBLEM 
Created on Thu Feb 29 12:37:43 2024
@author: Syed Tanweer Raza Nujjoo
Credits: Geomatics department in the University of Cape Town (UCT)
"""

# Importing Relevant libraries
import pandas as pd
import sympy as sp
import numpy as np


# (i) Simulated Data ------------------------------------------------------
# importing excel spreadsheet
df = pd.read_excel('Sequential polar obs.xlsx')
print(df)

# (ii) Computing Coordinates for Remaining Points -------------------------

# Polar function
def Polar_Coords (YCoords, XCoords, distance, direction_in_decimaldegrees):
    
    Y = YCoords + (distance * np.sin(np.radians(direction_in_decimaldegrees)))
    X = XCoords + (distance * np.cos(np.radians(direction_in_decimaldegrees)))
    
    return [round(Y,4), round(X,4)]

# Assumption:
P1Y = 1000 # let x coordinates of P1 be 1000
P1X = 1000 # let y coordinates be P1 be 1000

# calculating sequential polars using the polar function
# a loop will not be relevant in this case as sequential polar is solely subjective to the network design itself
P2_coords = Polar_Coords(P1Y, P1X, df.Dist[0], df.Direc[0])
P3_coords = Polar_Coords(P1Y, P1X, df.Dist[1], df.Direc[1])
P4_coords = Polar_Coords(P3_coords[0], P3_coords[1], df.Dist[2], df.Direc[2])
P5_coords = Polar_Coords(P3_coords[0], P3_coords[1], df.Dist[3], df.Direc[3])
P6_coords = Polar_Coords(P5_coords[0], P5_coords[1], df.Dist[4], df.Direc[4])
P7_coords = Polar_Coords(P5_coords[0], P5_coords[1], df.Dist[5], df.Direc[5])
P8_coords = Polar_Coords(P6_coords[0], P6_coords[1], df.Dist[6], df.Direc[6])
P9_coords = Polar_Coords(P6_coords[0], P6_coords[1], df.Dist[7], df.Direc[7])

# stack all coordinates
all_coords = np.vstack([[P1Y,P1X],
           P2_coords,
           P3_coords,
           P4_coords,
           P5_coords,
           P6_coords,
           P7_coords,
           P8_coords,
           P9_coords])

P_coords_df = pd.DataFrame(all_coords, index=['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9'])
print()

# rename the columns
P_coords_df.columns = ['Y(m)', 'X(m)']
print(P_coords_df)

# (iii) Covariance matrix of P2, P3, ..., P9 ------------------------------

# defining the symbols
Y, X = sp.symbols('Y X')
d12, d13, d34, d35, d56, d57, d68, d69= sp.symbols('d12 d13 d34 d35 d56 d57 d68 d69')
a12, a13, a34, a35, a56, a57, a68, a69= sp.symbols('a12 a13 a34 a35 a56 a57 a68 a69')

# defining observation equations
P2Y = Y+d12*sp.sin(a12)
P2X = X+d12*sp.cos(a12)

P3Y = Y+d13*sp.sin(a13)
P3X = X+d13*sp.cos(a13)

P4Y = Y+d34*sp.sin(a34)
P4X = X+d34*sp.cos(a34)

P5Y = Y+d35*sp.sin(a35)
P5X = X+d35*sp.cos(a35)

P6Y = Y+d56*sp.sin(a56)
P6X = X+d56*sp.cos(a56)

P7Y = Y+d57*sp.sin(a57)
P7X = X+d57*sp.cos(a57)

P8Y = Y+d68*sp.sin(a68)
P8X = X+d68*sp.cos(a68)

P9Y = Y+d69*sp.sin(a69)
P9X = X+d69*sp.cos(a69)

# defining the jacobian matrix
J = sp.Matrix([[P2Y],[P2X],
               [P3Y],[P3X],
               [P4Y],[P4X],
               [P5Y],[P5X],
               [P6Y],[P6X],
               [P7Y],[P7X],
               [P8Y],[P8X],
               [P9Y],[P9X]])

# defining the variables to which partial derivatives should be executed
var = [d12,a12,d13,a13,d34,a34,d35,a35,d56,a56,d57,a57,d68,a68,d69,a69]

# evaluting the jacobian matrix
jac = J.jacobian(var)
Jacs = jac.subs([(d12,df.Dist[0]),(a12,np.radians(df.Direc[0])),
                  (d13,df.Dist[1]),(a13,np.radians(df.Direc[1])),
                  (d34,df.Dist[2]),(a34,np.radians(df.Direc[2])),
                  (d35,df.Dist[3]),(a35,np.radians(df.Direc[3])),
                  (d56,df.Dist[4]),(a56,np.radians(df.Direc[4])),
                  (d57,df.Dist[5]),(a57,np.radians(df.Direc[5])),
                  (d68,df.Dist[6]),(a68,np.radians(df.Direc[6])),
                  (d69,df.Dist[7]),(a69,np.radians(df.Direc[7]))])

# defining the covariance matrix of observations
C_l = []
for i in range(len(df)):
    std_d = df.Sdist[i]
    std_a = np.radians(df.Sdirec_in_seconds[i]/3600)
    C_l.append(std_d)
    C_l.append(std_a)

# calculating the covariance matrix of parameters/points
cov_mat = Jacs * np.diagflat(C_l)**2 * Jacs.T
print()
print(cov_mat)






