#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GEOMATICS ENGINERRING NUMERICAL METHODS PROBLEM 
Created on Sun Feb 11 13:51:05 2024
@author: Syed Tanweer Raza Nujjoo
Credits: Geomatics department in the University of Cape Town (UCT)
"""

# Importing relevant libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import statistics as st

"""
X = [72.5, 83.4, 69.7, 69.5, 73.8, 82.6, 70.4, 80.8, 77.5, 80.4]
Y = [1.72, 1.84, 1.68, 1.65, 1.72, 1.81, 1.66, 1.77, 1.75, 1.80]
"""

# (i) Prompting User for Values and calculate specific statistics ---------
N = int(input("Enter number of distances/heights: ")) # prompt user to input no. of distances/heights 
print()

if N>=20:
    print("Suggestion: Number of values are too much to input by hand. Rather read the measurements from an excel sheet!")
else:
    X = [] # empty list to append distances X
    Y = [] # empty list to append heights Y
    
    for i in range(N):
        dist = float(input(f"Enter Distance {i+1}: ")) # prompt user to input distances
        X.append(dist) # append inputs to distance variable X
        height = float(input(f"Enter Height {i+1}: ")) # prompt user to input heights
        Y.append(height) # append inputs to height variable Y
        print()
  
    
    df = pd.DataFrame({'Distance': X, 'Height': Y}) # insert inputs in a dataframe
    
    mean_X = st.mean(df['Distance']) # sample mean for distance variable
    mean_Y = st.mean(df['Height']) # sample mean for height variable
    
    var_X = st.variance(df['Distance']) # sample variance for distance variable
    var_Y = st.variance(df['Height']) # sample variance for height variable
    
    std_X = np.sqrt(var_X) # standard deviation for the distance variable 
    # st.stdev(df['Distance']) # another alternative to calculate above
    std_Y = np.sqrt(var_Y) # standard deviation for the height variable
    # st.stdev(df['Height']) # another alternative to calculate above
    
    # (ii) Covariance of XY ---------------------------------------------------
    covmat_XY = np.cov(df['Distance'], df['Height']) # covariance calculation for distance and height variables
    cov_XY = covmat_XY[0][1] # extract covariance XY from covariance matrix above
    
    # (iii) Display Observations on separate barplots -------------------------
    plt.bar(x= list(range(1,N+1)), height = df['Distance'], width = 1, edgecolor = 'blue')
    plt.ylabel('Distance')
    plt.xlabel('N Observation')
    plt.show()
    
    plt.bar(x = list(range(1,N+1)), height = df['Height'], width = 1, edgecolor = 'blue')
    plt.xlabel('N Observation')
    plt.ylabel('Height')
    plt.show()
    
    # (iv) Substituting observations to dF ------------------------------------
    
    Xi,Yi = sp.symbols('Xi Yi') # defining symbols to use in the function below
    fX = Xi**3 + 2*Yi # defining function as an expression
    
    dfX = sp.diff(fX, Xi) # partial derivative w.r.t Xi
    dfX_sol = [dfX.subs(Xi, x) for x in df['Distance']] # evaluating the partial derivative w.r.t. Xi at each distance observation

# The end -----------------------------------------------------------------
        
        
        
        
        
































































