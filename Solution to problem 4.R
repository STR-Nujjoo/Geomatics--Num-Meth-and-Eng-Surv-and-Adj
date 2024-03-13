# GEOMATICS ENGINERRING NUMERICAL METHODS PROBLEM 
# Author: Syed Tanweer Raza Nujjoo
# Date: 11/03/2024
# Credits: Geomatics department in the University of Cape Town (UCT)

# Source showing various package to help calculating jacobians in R: https://ms.mcmaster.ca/~bolker/misc/jacobian.html

#rm(list = ls()) # Clear global environment

# Importing Relevant libraries

library(Deriv)
library(readxl)

# Import excel spreadsheet ------------------------------------------------
df <- read_excel('Levelling network info.xlsx')

# Known heights -----------------------------------------------------------
h1 <- 100
h2 <- 99.729
hAo <- h1+df$deltah_1A
hBo <- hAo+df$deltah_AB

# Defining Design/Jacobian matrix with observation equations --------------
J <- Deriv(~c(hAo-h1, # observation equation 1
              hBo-hAo, # observation equation 2
              h2-hBo), # observation equation 3
           c('hAo', 'hBo')) # jacobian w.r.t hA and hB

A <- matrix(J, 3, 2) # design matrix

# Defining the weight matrix ----------------------------------------------
weights <- c(1/(df$sd_1A)^2, 1/(df$sd_AB)^2, 1/(df$sd_B2)^2) # turning standard deviations into weights
P <- diag(weights) # weight matrixs

# Defining misclosure matrix ----------------------------------------------
misclosure1 <- expression(dh1-(hAo-h1)) # misclosure equation 1
misclosure2 <- expression(dh2-(hBo-hAo)) # misclosure equation 2
misclosure3 <- expression(dh3-(h2-hBo)) # misclosure equation 3

# defining terms
dh1 <- df$deltah_1A 
dh2 <- df$deltah_AB
dh3 <- df$deltah_B2

# evaluating misclosures
l1 <- eval(misclosure1) 
l2 <- eval(misclosure2)
l3 <- eval(misclosure3)

l <- matrix(c(l1,l2,l3), 3, 1) # misclosure matrix

# Solving for unknowns using parametric least squares formula
x <- solve(t(A) %*% P %*% A) %*% (t(A) %*% P %*% l)

# Adjusted heights
hA_adj <- hAo + x[1,] # adjusted height at A
hB_adj <- hBo + x[2,] # adjusted height at B

# The end -----------------------------------------------------------------