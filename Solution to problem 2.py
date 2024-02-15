#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GEOMATICS ENGINERRING NUMERICAL METHODS PROBLEM 
Created on Thu Feb 15 14:50:46 2024
@author: Syed Tanweer Raza Nujjoo
Credits: Geomatics department in the University of Cape Town (UCT)
"""

# Importing relevant libraries
import numpy as np
import sympy as sp
import math as m

# Defining Known Parameters and Equation ----------------------------------
obs_alp = np.radians(60) # angle alpha in radians
std_alp = np.radians(30/3600) # standard deviation of angle alpha in radians

obs_beta = np.radians(60) # angle beta in radians
std_beta = np.radians(2/60) # standard deviation of angle beta in radians

obs_C = 100 # side C in meters
std_C = 3/100 # standard deviation of side C in meters

alp, beta, C = sp.symbols('alp beta C') # defining symbols to use in the function below
f = C * (sp.sin(beta)/ (sp.sin(m.pi - (alp + beta)))) # given formula derived from area of triangle ABC

# (i) Calculating an actual value of b ------------------------------------
b = f.subs([(alp,obs_alp), (beta,obs_beta), (C,obs_C)]) # evaluating expression on known values defined above

# (ii) Computing jacobian matrix for expression relating to b -------------
mat = sp.Matrix([f]) # convert function into a matrix
var = [alp, beta, C] # define variables to which partial derivatives should be executed
J = mat.jacobian(var) # Jacobian w.r.t variables alpha, beta and C
jacs = J.subs([(alp,obs_alp), (beta,obs_beta), (C,obs_C)]) # jacobian matrix for the function defined as f- related to side b's expression

# (iii) Applying the law of propagation of variances to find std b --------
std_mat = np.matrix([std_alp, std_beta, std_C]) # display standard deviations in a matrix
std_b = m.sqrt(np.square(jacs) * np.square(std_mat).T) # formula derived in supporting document as Equation 3

# (iv) Finding highest error contributor from the error propagation -------
alp_comp = jacs[0]**2 * std_mat[0,0]**2 # evaluating the alpha component (or component 1) from the error propagation
beta_comp = jacs[1]**2 * std_mat[0,1]**2  # evaluating the beta component (or component 2) from the error propagation
C_comp = jacs[2]**2 * std_mat[0,2]**2 # evaluating the C component (or component 3) from the error propagation

list_comp = [alp_comp, beta_comp, C_comp] # put the above in a list for comparison

# comparison to declare largest error contributor
print()
if list_comp.index(max(list_comp)) == 0: # extracting the index of the maximum error
    print(f"(iv) From the comparison, {alp_comp} is the largest value, therefore angle α affects the propagation the most.")
elif list_comp.index(max(list_comp)) == 1: # extracting the index of the maximum error
    print(f"(iv) From the comparison, {beta_comp} is the largest value, therefore angle ß affects the propagation the most.")
else: # extracting the index of the maximum error
    print(f"(iv) From the comparison, {C_comp} is the largest value, therefore side c affects the propagation the most.")

# (v) How to reduce propagated variance -----------------------------------
print()
print("(v) The logical way to improve the result is to adopt a more precise procedure of observing the direction using circle right and circle left. Also, the standard deviations of both angles are very different and that is quite unusual unless the observer used 2 different instruments to measure the angle. In addition, this should be avoided but if the observer chose to observe with 2 instruments then he should consider adopting weightings in his calculations.")

# The end -----------------------------------------------------------------