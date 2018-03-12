lcor <- function(df, X, Y, subgroup=NULL){
  cor <-lcov(df, X, Y, subgroup)/(lcov(df, X, X, subgroup)*lcov(df, Y, Y, subgroup))^0.5
  cor
}

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

lcov <- function(df, X, Y, subgroup=NULL){
  D <- calc_dyads(df, X, Y)
  if (!is.null(subgroup)){
    D <- remove_subgroup(df, X, Y, subgroup)
  }
  cov <- cov_from_D(D)
  cov
}

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