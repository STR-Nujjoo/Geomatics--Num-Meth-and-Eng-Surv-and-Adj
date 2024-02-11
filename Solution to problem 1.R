# GEOMATICS ENGINERRING NUMERICAL METHODS PROBLEM 
# Author: Syed Tanweer Raza Nujjoo
# Date: 11/02/2024
# Credits: Geomatics department in the University of Cape Town (UCT)


# (i) Prompting User for Values and calculate specific statistics ---------
N <- readline("Enter number of distances/heights: ") # prompt user to input no. of distances/heights
N <- as.numeric(N) # convert the latter to the numeric class

if (N>=20){
  readline("Suggestion: Number of values are too much to input by hand. Rather read the measurements from an excel sheet!")
}

X <- NULL # empty list to append distances X
Y <- NULL # empty list to append heights Y
for (i in 1:N){
  dist <- readline(paste0("Enter distance ",i, ": ")) # prompt user to input distances
  X[i] <- as.numeric(dist) # append inputs to distance variable X
  height <- readline(paste0("Enter height ",i, ": ")) # prompt user to input heights
  cat("\n")
  Y[i] <- as.numeric(height) # append inputs to height variable Y
}

# X <- c(72.5, 83.4, 69.7, 69.5, 73.8, 82.6, 70.4, 80.8, 77.5, 80.4)
# Y <- c(1.72, 1.84, 1.68, 1.65, 1.72, 1.81, 1.66, 1.77, 1.75, 1.80)

df <- data.frame(Distance = X, Height = Y) # insert inputs in a dataframe
# summary(df)

mean_X <- mean(df$Distance) # sample mean for distance variable
mean_Y <- mean(df$Height) # sample mean for height variable

var_X <- var(df$Distance) # sample variance for distance variable
var_Y <- var(df$Height) # sample variance for height variable

std_X <- sqrt(var_X) # standard deviation for the distance variable 
std_Y <- sqrt(var_Y) # standard deviation for the height variable


# (ii) Covariance of XY ---------------------------------------------------
cov_XY <- cov(df$Distance, df$Height) # covariance calculation for distance and height variables


# (iii) Display Observations on separate barplots -------------------------
barplot(df$Distance, 
        col = "steelblue",
        xlab = "N Observation",
        ylab = "Distance",
        names.arg = 1:10,
        space = 0)

barplot(df$Height, 
        col = "steelblue",
        xlab = "N Observation",
        ylab = "Height",
        names.arg = 1:10,
        space = 0)


# (iv) Substituting observations to dF ------------------------------------

fX <- expression(Xi^3+ 2*Yi) # defining function as an expression
dfX <- D(fX,'Xi') # partial derivative w.r.t Xi
Xi <- df$Distance # defining values for Xi
eval(dfX) # evaluating the partial derivative w.r.t. Xi at each distance observation

# The end -----------------------------------------------------------------



