library(lars)
library(tidyverse)
library(latex2exp)
# library(glmnet)
# `glmnet' has convergence issue when p>>n or p is large
# `lars' is more stable and as fast


get_xi <- function(w, h, n){
  # given omega, n and hkk, return xi
  ans = n*(1-w)/(n-1+w-n*(1-w)*h)
  abs(ans)
}


get_w <- function(xi, h, n){
  (n-xi*(n-1-n*h))/(n+(1+n*h)*xi)
}


centralize <- function(x, r = F){
  # centralize x with mean 0 and sd 1/sqrt(n-1)
  n = dim(x)[1]
  meanx <- drop(rep(1,n) %*% x)/n
  x = scale(x, meanx, FALSE)
  normx = sqrt(drop(rep(1,n) %*% (x^2)))
  if (r){
    return(list(m = meanx, d = scale(x,FALSE,normx), v = normx)) 
  }
  scale(x, FALSE, normx) 
}