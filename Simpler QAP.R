# Calculates the correlation between sharing X and sharing Y from a dataframe 
#
# Takes:
# df = A dataframe
# X  = A string with the name of the variable X (factor) to compute adjacencies
# Y  = A string with the name of the variable Y (factor) to compute adjacencies
# 
# Returns:
# The Person r between sharing attribute X and sharing attribute Y
lcor <- function(df, X, Y){
  cor <-lcov(df, X, Y)/(lcov(df, X, X)*lcov(df, Y, Y))^0.5
  cor
}

# Calculates the number of adjacencies from the dataframes variable X, Y and the combination of X and Y
#
# Takes:
# df = A dataframe
# X  = A string with the name of the variable X (factor) to compute adjacencies
# Y  = A string with the name of the variable Y (factor) to compute adjacencies
#
# Returns a "D" object that contains:
# dyads_tot      = The total number of possible dyadic relationships in df
# dyads_x        = The number of dyads that share X in df
# dyads_y        = The number of dyads that share Y in df
# dyads_xY       = The number of dyads that share both X and Y in df
# dyads_xyblanck = The number of dyads that do not share X or Y
calc_dyads <- function(df, X, Y){
  x <- summary(df[,X])
  y <- summary(df[,Y])
  xy <-summary(as.factor(paste(df[,X], df[,Y])))
  D <- list()
  D$dyads_tot <- nrow(df)^2-nrow(df)
  D$dyads_x <- sum(x^2-x)
  D$dyads_y <- sum(y^2-y)
  D$dyads_xy <- sum(xy^2-xy)
  D$dyads_xyblanck <- (D$dyads_tot-D$dyads_xy-(D$dyads_x-D$dyads_xy)-(D$dyads_y-D$dyads_xy))
  D
}

# Calculates the covariance from the information on the adjacencies
#
# Takes a "D" object that contains:
# dyads_tot      = The total number of possible dyadic relationships in df
# dyads_x        = The number of dyads that share X in df
# dyads_y        = The number of dyads that share Y in df
# dyads_xY       = The number of dyads that share both X and Y in df
# dyads_xyblanck = The number of dyads that do not share X or Y
#
# Returns:
# The covariance between sharing attribute X and sharing attribute Y
cov_from_D <- function(D){
  x_mean <- D$dyads_x/D$dyads_tot
  y_mean <- D$dyads_y/D$dyads_tot
  cov_xy_positive <- D$dyads_xy*(1-x_mean)*(1-y_mean)
  cov_xy_negativ <- D$dyads_xyblanck*-x_mean*-y_mean
  cov_only_x <- (D$dyads_x-D$dyads_xy)*(1-x_mean)*(-y_mean)
  cov_only_y <- (D$dyads_y-D$dyads_xy)*(-x_mean)*(1-y_mean)
  cov <- (cov_xy_positive+cov_xy_negativ+cov_only_x+cov_only_y)/(D$dyads_tot-1)
  cov
}

# Calculates the covariance between sharing X and sharing Y from a dataframe 
#
# Takes:
# df = A dataframe
# X  = A string with the name of the variable X (factor) to compute adjacencies
# Y  = A string with the name of the variable Y (factor) to compute adjacencies
#
# Returns:
# The covariance between sharing attribute X and sharing attribute Y
lcov <- function(df, X, Y){
  D <- calc_dyads(df, X, Y)
  cov <- cov_from_D(D)
  cov
}

# Calculates the correlation between sharing two different attributes and test the null hypothesis.
# The correlation is calculated without generating matrices, but with dyad information inferred from the transitivity. 
# Opposed to the quadratic assignment procedure, the reordering happens on the dependent attributes. And null hypothesis coefficients is calculated from this reordering as on the real value.
#
# Takes:
# df = A dataframe
# X  = A string with the name of the variable X (factor) to compute adjacencies
# Y  = A string with the name of the variable Y (factor) to compute adjacencies
#
# Returns a object similar to that of qaptest (sna v2.4):
# testval = Persons r between adjacencies of X and adjacencies of Y 
# dist    = A vector containing the Monte Carlo draws
# pgreq   = The proportion of draws which were greater than or equal to the observed value
# pleeq   = The proportion of draws which were less than or equal to the observed value
lap <- function(df, X, Y, reps=1000){
  testval <- lcor(df, X, Y)
  dist <- numeric()
  for(p in 1:reps){
    Perm <- data.frame(df[,X], sample(df[,Y]))
    dist[p] <- lcor(Perm, 1, 2)
  }
  pgreq <- sum(dist >= testval)/reps
  pleeq <- sum(dist <= testval)/reps
  out <- list()
  out$testval <- testval; out$dist <- dist; out$pgreq <- pgreq; out$pleeq <- pleeq
  out
}
