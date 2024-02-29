# GEOMATICS ENGINERRING NUMERICAL METHODS PROBLEM 
# Author: Syed Tanweer Raza Nujjoo
# Date: 24/02/2024
# Credits: Geomatics department in the University of Cape Town (UCT)

# Source showing various package to help calculating jacobians in R: https://ms.mcmaster.ca/~bolker/misc/jacobian.html

# Importing Relevant libraries
library(Deriv)
library(readxl)

# (i) Simulated Data ------------------------------------------------------
# importing excel spreadsheet
df <- read_excel('Sequential polar obs.xlsx')

# (ii) Computing Coordinates for Remaining Points -------------------------
# Degrees to radians function
rad <- function(angle){ # creating a function to convert degrees to radians
  angle * (pi/180)
}

# Polar function
Polar_Coords <- function(YCoords,XCoords,distance,direction_in_decimaldegrees){
  
  Y <- YCoords + (distance * sin(rad(direction_in_decimaldegrees)))
  X <- XCoords + (distance * cos(rad(direction_in_decimaldegrees)))
  return(c(newY = Y, newX = X))
}

# Assumption:
P1Y <- 1000 # let x coordinates of P1 be 1000
P1X <- 1000 # let y coordinates be P1 be 1000

# calculating sequential polars using the polar function
# a loop will not be relevant in this case as sequential polar is solely subjective to the network design itself
P2_coords <- Polar_Coords(P1Y, P1X, df$Dist[1], df$Direc[1])
P3_coords <- Polar_Coords(P1Y, P1X, df$Dist[2], df$Direc[2])
P4_coords <- Polar_Coords(P3_coords[1], P3_coords[2], df$Dist[3], df$Direc[3])
P5_coords <- Polar_Coords(P3_coords[1], P3_coords[2], df$Dist[4], df$Direc[4])
P6_coords <- Polar_Coords(P5_coords[1], P5_coords[2], df$Dist[5], df$Direc[5])
P7_coords <- Polar_Coords(P5_coords[1], P5_coords[2], df$Dist[6], df$Direc[6])
P8_coords <- Polar_Coords(P6_coords[1], P6_coords[2], df$Dist[7], df$Direc[7])
P9_coords <- Polar_Coords(P6_coords[1], P6_coords[2], df$Dist[8], df$Direc[8])


all_coords <- rbind(c(P1Y, P1X), P2_coords, P3_coords, P4_coords, P5_coords, P6_coords, 
      P7_coords, P8_coords, P9_coords) # put the coordinates together

P_coords_df <- data.frame(all_coords, row.names = c('P1','P2','P3','P4',
                                                    'P5','P6','P7','P8','P9')) # put the coordinates in a dataframe format

colnames(P_coords_df) <- c('Y(m)', 'X(m)') # rename the columns 


# (iii) Covariance matrix of P2, P3, ..., P9 ------------------------------
# Defining jacobian matrix
J <- Deriv(~c(P1Y+d12*sin(a12), P1X+d12*cos(a12), # P2 observation equation for the Y and X coordinates
              P1Y+d13*sin(a13), P1X+d13*cos(a13), # P3 observation equation for the Y and X coordinates
              P3Y+d34*sin(a34), P3X+d34*cos(a34), # P4 observation equation for the Y and X coordinates
              P3Y+d35*sin(a35), P3X+d35*cos(a35), # P5 observation equation for the Y and X coordinates
              P5Y+d56*sin(a56), P5X+d56*cos(a56), # P6 observation equation for the Y and X coordinates
              P5Y+d57*sin(a57), P5X+d57*cos(a57), # P7 observation equation for the Y and X coordinates
              P6Y+d68*sin(a68), P6X+d68*cos(a68), # P8 observation equation for the Y and X coordinates
              P6Y+d69*sin(a69), P6X+d69*cos(a69)), # P9 observation equation for the Y and X coordinates
          # distance and angle observations- order in which partial derivatives will takes place
            c('d12','a12', 
                                                     'd13','a13',
                                                     'd34','a34',
                                                     'd35','a35',
                                                     'd56','a56',
                                                     'd57','a57',
                                                     'd68','a68',
                                                     'd69','a69'))

d_var <- c('d12','d13','d34','d35','d56','d57','d68','d69') # define distance variables in the same order in which observations occurred
a_var <- c('a12','a13','a34','a35','a56','a57','a68','a69') # define direction variables in the same order in which observations occurred
d_obs <- df$Dist # distance vector
a_obs <- rad(df$Direc) # direction vector

# Assign values from the vector to variables
for (i in 1:nrow(df)) {
  assign(d_var[i], d_obs[i])
  assign(a_var[i], a_obs[i])
}

jac <- eval(J) # evaluating jacobian matrix on known observations defined
Jacs <- matrix(jac,
               nrow=length(d_var) + length(a_var),
               ncol=length(d_var) + length(a_var)) # convert jacobian evaluation in a matrix

C_l <- diag(c(rbind(df$Sdist^2,rad(df$Sdirec_in_seconds/3600)^2))) # covariance matrix of observations
cov_mat <- Jacs %*% C_l %*% t(Jacs) # covariance matrix of parameters/points

