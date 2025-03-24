
library(lars)
library(tidyverse)
library(latex2exp)
library(Matrix)
library(MASS)
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
    return(list(d = scale(x,FALSE,normx), m = meanx, v = normx))
  }
  scale(x, FALSE, normx)
}


qr_delete = function(Q,R,col){
  # function test
  # A = matrix(c(1, 2,3, 4,
  #              0, 4, 1, 5,
  #              0, 2, 2,9,
  #              1,  0,0,11),
  #            nrow = 4, byrow = TRUE)
  # QR = qr(A)
  # PP = qr_delete(qr.Q(QR),qr.R(QR),2)
  # PP$Q%*%PP$R
  if (dim(R)[1]<col){
    warning("no such feature included", call. = FALSE)
  }
  p = dim(R)[1]
  R = R[,-col,drop=F]
  if (col==p){
    return(list(Q=Q[,-p,drop=F],R=R[-p,,drop=F]))
  }
  if (col <= (p-2)){
    for (i in col:(p-2)){
      a = R[i,i]
      b = R[i+1,i]
      rot = matrix(c(a/sqrt(a**2+b**2), -b/sqrt(a**2+b**2),
                     b/sqrt(a**2+b**2), a/sqrt(a**2+b**2)), nrow=2)
      Q[,i:(i+1)] = Q[,i:(i+1)]%*%t(rot)
      R[i:(i+1),] = rot%*%R[i:(i+1),]
    }
  }
  temp = sqrt(R[p,p-1]**2 + R[p-1,p-1]**2)
  Q[,(p-1)] = (R[p-1,p-1] * Q[,(p-1)] + R[p,p-1] * Q[,p])/temp
  R[p-1,p-1] = temp
  return(list(Q=Q[,-p,drop=F],R=R[-p,,drop=F]))
}


qr_insert = function(Q,R,col,v){
  # test code
  # A = matrix(c(1, 2,3,
  #              0, 4, 1,
  #              0, 2, 2,
  #              1,  0,0),
  #            nrow = 4, byrow = TRUE)
  # QR = qr(A)
  # PP = qr_insert(qr.Q(QR),qr.R(QR),1,c(1,2,3,4))
  # PP$Q%*%PP$R
  p = dim(R)[1]
  if (is.null(p)){
    nor = sqrt(sum(v**2))
    return(list(Q= as.matrix(v/nor,length(v),1), R = as.matrix(nor,1,1)))
  }
  if (col>(1+p)){
    warning("there is not enough features before", call. = FALSE)
  }
  if (col == (1+p)){
    coor = t(Q)%*%v
    res = v - Q%*%coor
    nor = sqrt(sum(res**2))
    return(list(Q=cbind(Q, res/nor), R = cbind(rbind(R,0), c(coor, nor))))
  }
  new_coor = t(Q)%*%v
  new_q = v - Q%*%new_coor
  nor = sqrt(sum(new_q**2))
  new_q = new_q/nor
  Q = cbind(Q, new_q)
  if (col == 1){
    R = rbind(cbind(new_coor,R), c(nor,rep(0,p)))
  }
  else{
    R = rbind(cbind(R[,1:(col-1)], new_coor, R[,col:p]), c(rep(0, col-1), nor, rep(0,p-col+1)))
  }

  for (i in (p+1):(col+1)){
    a = R[i-1,col]
    b = R[i,col]
    rot = matrix(c(a/sqrt(a**2+b**2), -b/sqrt(a**2+b**2),
                   b/sqrt(a**2+b**2), a/sqrt(a**2+b**2)), nrow=2)
    Q[,(i-1):i] = Q[,(i-1):i]%*%t(rot)
    R[(i-1):i,] = rot%*%R[(i-1):i,]
  }
  return(list(Q = Q, R = R))
}

