# GEOMATICS ENGINERRING NUMERICAL METHODS PROBLEM 
# Author: Syed Tanweer Raza Nujjoo
# Date: 14/02/2024
# Credits: Geomatics department in the University of Cape Town (UCT)

# Source showing various package to help calculating jacobians in R: https://ms.mcmaster.ca/~bolker/misc/jacobian.html

# Importing relevant libraries
library(Deriv)

# Defining Known Parameters and Equation ----------------------------------
rad <- function(angle){ # creating a function to convert degrees to radians
  angle * (pi/180)
}

alp <-  rad(60) # angle alpha in radians
std_alp <-  rad((30/60)/60) # standard deviation of angle alpha in radians

beta <- rad(60) # angle beta in radians
std_beta <-  rad(2/60) # standard deviation of angle beta in radians

C <-  100 # side C in meters
std_C <-  3/100 # standard deviation of side C in meters

f <-  expression(C * (sin(beta)/(sin(pi-(alp+beta))))) # given formula derived from area of triangle ABC

# (i) Calculating an actual value of b ------------------------------------
b <-  eval(f) # evaluating expression on known values defined above

# (ii) Computing jacobian matrix for expression relating to b -------------
J <- Deriv(f, c("alp", "beta", "C")) # Jacobian w.r.t variables alpha, beta and C
jacs <- matrix(eval(J)) # jacobian matrix for the function defined as f- related to side b's expression

# (iii) Applying the law of propagation of variances to find std b --------
std_b = sqrt(t(jacs^2) %*% matrix(c(std_alp^2, std_beta^2, std_C^2))) # formula derived in supporting document as Equation 3


# (iv) Finding highest error contributor from the error propagation -------
alp_comp <- jacs[1,]^2 * std_alp^2 # evaluating the alpha component (or component 1) from the error propagation
beta_comp <- jacs[2,]^2 * std_beta^2 # evaluating the beta component (or component 2) from the error propagation
C_comp <- jacs[3,]^2 * std_C^2 # evaluating the C component (or component 3) from the error propagation

if (which.max(c(alp_comp, beta_comp, C_comp)) == 1){
  cat("From the comparison,", alp_comp, "is the largest value, therefore angle α affects the propagation the most.")
}else if (which.max(c(alp_comp, beta_comp, C_comp)) == 2){
  cat("From the comparison,", beta_comp, "is the largest value, therefore angle ß affects the propagation the most.")
}else {
  cat("From the comparison,", C_comp, "is the largest value, therefore side c affects the propagation the most.")
} # comparison to declare largest error contributor


# (v) How to reduce propagated variance -----------------------------------
cat("The logical way to improve the result is to adopt a more precise procedure of observing the direction using circle right and circle left. Also, the standard deviations of both angles are very different and that is quite unusual unless the observer used 2 different instruments to measure the angle. In addition, this should be avoided but if the observer chose to observe with 2 instruments then he should consider adopting weightings in his calculations.")

# The end -----------------------------------------------------------------