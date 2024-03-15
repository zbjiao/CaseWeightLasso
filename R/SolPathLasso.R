source('../CaseWeightLasso/R/common.R')


#' obtain the solution path of observation k when lambda is fixed. 
#' @param X           matrix n by p      design matrix
#' @param y           vector n by 1      response vector
#' @param k           integer            observation of interest
#' @param s           a value, or vector of values, indexing the penalty levels to consider.
#'                    Its values depends on the mode= argument. By default (mode="fraction"), 
#'                    s should take on values between 0 and 1. 
#' @param mode        Mode="fraction", then s should be a number between 0 and 1, and it refers 
#'                    to the ratio of the L1 norm of the coefficient vector, relative to the 
#'                    norm at the full LS solution if n>p or the minimum L1 norm linear 
#'                    interpolation solution if n<=p. Mode="norm" means s refers to the L1 norm 
#'                    of the coefficient vector. Mode="lambda" uses the lasso regularization 
#'                    parameter for s. Abbreviations allowed.
#' @return \item{w_path}{a vector of breakpoints}
#'         \item{hkk_path}{leverages at each breakpoint}
#'         \item{beta_path}{beta estimate at each breakpoint}
#'         \item{s_path}{beta hat's sign at each breakpoint}
#'         \item{beta0_path}{beta0 at each breakpoint}
#'         \item{l1norm}{l1 norm of beta at each breakpoint}
#'         \item{normx}{standard deviation of X by column}
#'         \item{meanx}{mean of X by column}
#'         \item{obs_of_interest}{k}
#'         \item{special}{A boolean recording if penalty is 0 or not}
#' @description
#' Only predictors that have sign change are included in the solution path plot
#' @seealso predict, plot, coef methods
#' @examples
#' library(lars)
#' data(diabetes)
#' attach(diabetes)
#' obj = SolPathLooLasso(x,y,k = 182, 5, 'lambda')
#' 
#' detach(diabetes)
#' 
#' set.seed(100)
#' x = matrix(rnorm(50*200),nrow=50)
#' y = x[,1:5]%*%c(5,4,3,2,1) + rnorm(50)
#' ans = SolPathLooLasso(x,y,k = 1, s = 0.7, mode = "fraction")
#' plot(ans)
#' predict(ans)
#' coef(ans)
#' @export
SolPathLooLasso <- function(X, y, k = 1, s, mode = c("fraction", "norm", "lambda")){
  
  mode <- match.arg(mode)
  X_ = centralize(X, r = T)
  X = X_$d
  meanx = X_$m
  normx = X_$v
  

  n = dim(X)[1]
  p = dim(X)[2]
  
  if (k > n) {
    stop("index k out of bounds", call. = FALSE)
  }
  if (is.null(k)){
    stop("please specify the index of interest", call. = FALSE)
  }
  if (n!=length(y)){
    stop("X and y don't match", call. = FALSE)
  }
  if(is.null(colnames(X))){
    colnames(X) = 1:p
  }
  if (( (s==0 & mode=="lambda") | (s==1 & mode == "fraction") ) & n<p){
    stop("it is beyond the scope of this function", call. = FALSE)
  }
  
  # obtain lasso solution
  # fit = glmnet(X,y,family="gaussian",lambda=lambda/n,standardize=FALSE,thresh=1e-16,intercept=T)
  # beta_hat = as.vector(fit$beta)
  obj = lars(X,y,type='lasso',use.Gram = !(n<p | p>500))
  # beta_hat = coef(obj,s=lambda,mode='lambda')
  beta_hat = as.vector(predict.lars(obj,mode=mode, s=s, type = 'coefficients')$coefficients)
  
  if (mode == 'lambda'){
    lambda = s}
  else{
    nbeta1 = drop(abs(obj$beta) %*% rep(1, p))
    nbeta2 = drop(rep(1, p) %*% abs(beta_hat))
    lambda = approx(nbeta1, c(obj$lambda,0), xout = nbeta2)$y
  }
  
  
  ybar = mean(y)
  # record the sign of each covariate 
  s = sign(beta_hat)
  # active set 
  A = s!=0
  
  # derivative of f in terms of covariates in nonactive set
  d_hat = t(X[,!A])%*%y - t(X[,!A])%*%X%*%beta_hat
  if (sum(A)==0){
    XX = matrix(0,nrow=p,ncol=p)
    hk = as.matrix(rep(0,n),nrow=n)
  }
  else{
    XX = solve(crossprod(X[,A]))
    hk = X[,A]%*%XX%*%X[k,A]
  }
  
  # current predictor of yk 
  yk_hat = drop(ybar + X[k,] %*% beta_hat)
  
  # eazy case
  if (lambda == 0 & n>p){
    beta1 = as.vector(coef(lm(y[-k]~X[-k,])))
    beta10 = beta1[1]
    beta1 = beta1[2:(p+1)]
    ans = list(w_path = c(1,0),hkk_path = c(hk[k], hk[k]), 
               beta_path = rbind(beta_hat, beta1), s_path = rbind(s,s = sign(beta1)), 
               beta0_path = c(ybar, beta10),l1norm = sum(abs(beta_hat)),
               normx = normx, meanx = meanx, obs_of_interest = k, special=T)
    class(ans) <- "SolPathLooLasso"
    return(invisible(ans))
  } 

  hk_path = c()
  # beta path records beta hat's value at each breakpoint
  beta_path = c(beta_hat)
  # so does beta0_path records intercept
  beta0_path = c(ybar)
  # and sign change
  s_path = c(s)
  # and omega value at each breakpoint
  w_path = c()
  hkk_path = c()
  w = 1
  while (T){
    hk_path = cbind(hk_path,hk)
    w_path = c(w_path, w)
    hkk_path = c(hkk_path, hk[k])
    xi = get_xi(w, hk[k], n)
    bias = yk_hat - y[k]
    if (sum(A) == 0){
      xi_cand1 = c()
    }
    else{
      slope_beta = XX%*%X[k,A]*bias
      xi_cand1 = -XX%*%(t(X[,A])%*%y-lambda*s[A])/slope_beta
    }
    if (sum(A) == (n-1)){
      slope_d = rep(0, p-sum(A))
      xi_cand2p = rep(Inf, p-sum(A))
      xi_cand2m = rep(Inf, p-sum(A))
    }
    else{
      slope_d = (X[k,!A] - t(X[,!A])%*%hk)*bias
      xi_cand2p = (-lambda-d_hat)/slope_d
      xi_cand2m = (lambda-d_hat)/slope_d
    }
    # xi candidates
    # xi_cand1 = -beta_hat[A]/slope_beta
    xi_cand = c(xi_cand1, xi_cand2p, xi_cand2m)
    xi_cand0 = min(xi_cand[xi_cand>(xi+0.00001)],Inf)
    ind = which(xi_cand == xi_cand0)
    # update beta
    if (sum(A) > 0){
      beta_hat[A] = beta_hat[A] + min(get_xi(0,hk[k],n),xi_cand0)*slope_beta
    }
    beta_path = rbind(beta_path, beta_hat)
    beta_hat0 = ybar + min(get_xi(0,hk[k],n),xi_cand0) * bias/n
    beta0_path = c(beta0_path, beta_hat0)
    
    # if the xi is off the bound, stop the algorithm
    if (xi_cand0 > get_xi(0,hk[k],n)){
      w_path = c(w_path, 0)
      hkk_path = c(hkk_path, hk[k])
      s_path = rbind(s_path, s)
      hk_path = cbind(hk_path, hk)
      break
    }
    # if not, locate the covariate (go to or leave the active set)
    else if(ind<=sum(A)){
      coef = which(cumsum(A)==ind)[1]
      A[coef] = F
      s[coef] = 0
      # check if active set is empty 
      if (sum(A) == 0){
        hk_path = cbind(hk_path, 0,0)
        w_path = c(w_path, get_w(xi_cand0, hk[k], n), 0)
        hkk_path = c(hkk_path, 0, 0)
        s_path = rbind(s_path, s, s)
        beta_path = rbind(beta_path, beta_hat)
        beta0_path = c(beta0_path, ybar + get_xi(0,0,n) * (ybar-y[k])/n)
        break
      }
    }
    else if(ind>p){
      coef = which(cumsum(1-A)==(ind-p))[1]
      A[coef] = T
      s[coef] = 1
    }
    else{
      coef = which(cumsum(1-A)==(ind-sum(A)))[1]
      A[coef] = T
      s[coef] = -1
    }
    
    s_path = rbind(s_path, s)
    # update omega, XX, beta_hat, d_hat
    w = get_w(xi_cand0, hk[k], n)
    XX = solve(t(X[,A])%*%X[,A])
    beta_hat[A] = XX%*%(t(X[,A])%*%y-lambda*s[A])
    beta_hat[!A] = 0
    d_hat = t(X[,!A])%*%y - t(X[,!A])%*%X%*%beta_hat
    # y_hat = ybar + X%*%beta_hat
    hk = X[,A]%*%XX%*%X[k,A]
    yk_hat = drop(ybar + X[k,] %*% beta_hat)
  }
  rownames(beta_path) <- NULL
  rownames(s_path) <- NULL
  
  ans = list(w_path = w_path, hkk_path = hkk_path, beta_path = beta_path, 
              s_path = s_path, beta0_path = beta0_path,
              l1norm = as.vector(apply(abs(beta_path), 1, sum)),
             normx = normx, meanx = meanx, obs_of_interest = k,
             special = F)
  class(ans) <- "SolPathLooLasso"
  invisible(ans)
}

#' extract coefficients or estimates from a SolPathLooLasso object 
#' @param obj         fitted SolPathLooLasso object 
#' @param newx        Matrix of new values for x at which predictions are to be made. 
#'                    Must be a matrix.
#'                    This argument is not used for type="coefficients"
#' @param type        Type of prediction required. 
#'                    Type "coefficients" computes the LOO Lasso estimate 
#'                    Type "fit" computes the y estimates of given newx
#' @return The object returned depends on type.
#' @seealso plot, coef methods
#' @examples
#' library(lars)
#' data(diabetes)
#' attach(diabetes)
#' obj = SolPathLooLasso(x,y,k = 182, 5, 'lambda')
#' 
#' detach(diabetes)
#' 
#' set.seed(100)
#' x = matrix(rnorm(50*200),nrow=50)
#' y = x[,1:5]%*%c(5,4,3,2,1) + rnorm(50)
#' ans = SolPathLooLasso(x,y,k = 1, s = 0.7, mode = "fraction")
#' plot(ans)
#' predict(ans)
#' coef(ans)
#' @export
## S3 method for class 'SolPathLooLasso'
predict.SolPathLooLasso <- function(obj, newx, type = c("fit", "coefficients")){
  type <- match.arg(type)
  if(missing(newx) & type == "fit") {
    warning("Type=fit with no newx argument; type switched to coefficients"
    )
    type <- "coefficients"
  }
  L = length(obj$w_path)
  ans <- switch(type,
                    coefficients = list(obs_of_interest = obj$obs_of_interest,
                                        loobeta = drop(obj$beta_path[L,]),
                                        loobeta0 = drop(obj$beta0_path[L])),
                    fit = list(obs_of_interest = obj$obs_of_interest,
                                fit = drop(scale(newx, obj$meanx, obj$normx) %*% 
                                             obj$beta_path[L,]) + obj$beta0_path[L]))
  ans
}

#' plot the solution path from a SolPathLooLasso object 
#' @param object      fitted SolPathLooLasso object 
#' @param plot        If plot=1, this function plots the linear approximated solution path.
#'                    If plot=2, plot exact solution path.
#' @return The solution path plot
#' @seealso predict, coef methods
#' @examples
#' library(lars)
#' data(diabetes)
#' attach(diabetes)
#' obj = SolPathLooLasso(x,y,k = 182, 5, 'lambda')
#' 
#' detach(diabetes)
#' 
#' set.seed(100)
#' x = matrix(rnorm(50*200),nrow=50)
#' y = x[,1:5]%*%c(5,4,3,2,1) + rnorm(50)
#' ans = SolPathLooLasso(x,y,k = 1, s = 0.7, mode = "fraction")
#' plot(ans)
#' predict(ans)
#' coef(ans)
#' @export
## S3 method for class 'SolPathLooLasso'
plot.SolPathLooLasso <- function(object, plot = 2){
  k = object$obs_of_interest
  plot_helper <-function(x, df){
    i = findInterval(-x, -w_path, rightmost.closed = T)
    beta1 = df[i]
    beta2 = df[i+1]
    w1 = w_path[i]
    w2 = w_path[i+1]
    hkk = hkk_path[i]
    beta1+(beta2-beta1)*(get_xi(x, hkk, n) - get_xi(w1, hkk, n))/(get_xi(w2, hkk, n) - get_xi(w1, hkk, n))
  }
  beta_path = object$beta_path
  w_path = object$w_path
  hkk_path = object$hkk_path
  
  coln = 1:length(beta_path[1,])
  
  # ind = c(which((beta_hat>0) & (beta1<0)),  which((beta_hat<0) & (beta1>0)))
  num_z = apply(beta_path, 2, function(c) sum(abs(c)< 1e-10))
  ind = which(num_z>0 & num_z<length(w_path))
  
  if(object$special){
    if (length(ind)>0){
      if (length(ind) == 1){
        xx = data.frame(coef = coln[ind])
        print(ggplot()+xlim(0,1) + xlab(TeX("Case weight $\\omega$")) + 
                ylab(TeX("Coefficient estimate $\\hat{\\beta}$")) + 
                geom_function(data = xx, fun = plot_helper0, args = list(ind),aes(color = coef))+ 
                ggtitle(paste("Solution path for Case",toString(k)))+
                theme(plot.title = element_text(hjust = 0.5),
                      panel.background = element_rect(fill = "white"),
                      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)))
      }
      else{
        fig = ggplot()+xlim(0,1)
        for (j in 1:length(ind)){
          xx = data.frame(coef = coln[ind[j]])
          fig = fig + geom_function(data = xx, fun = plot_helper0, args = list(ind[j]),
                                    aes(color = coef)) 
        }
        print(fig + xlab(TeX("Case weight $\\omega$")) + ylab(TeX("Coefficient estimate $\\hat{\\beta}$")) + ggtitle(paste("Solution path for Case",toString(k))) +
                theme(plot.title = element_text(hjust = 0.5),
                      panel.background = element_rect(fill = "white"),
                      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)))
      }
    }
    else{
      print(paste('there is no sign change for case',toString(k)))
    }
  }
  else{
    if (length(ind)>0){
      if(plot ==1){
        df = cbind(beta_path[,ind], w_path)
        colnames(df) = c(coln[ind], 'w')
        df = as_tibble(df) %>% gather('coef', 'val', -w)
        print(ggplot(df, aes(x = w, y = val, group=coef, color = coef))+
                geom_line()+ggtitle(paste("Approx Solution path for Case",toString(k)))+
                ylab(TeX("Coefficient estimate $\\hat{\\beta}$")) +
                xlab(TeX("Case weight $\\omega$")) +
                theme(plot.title = element_text(hjust = 0.5),
                      panel.background = element_rect(fill = "white"),
                      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)))
      }
      else{
        if (length(ind) == 1){
          xx = data.frame(coef = coln[ind])
          print(ggplot()+xlim(0,1) + xlab(TeX("Case weight $\\omega$")) + 
                  ylab(TeX("Coefficient estimate $\\hat{\\beta}$")) + 
                  geom_function(data = xx, fun = plot_helper, args = list(beta_path[,ind]), aes(color = coef))+ 
                  ggtitle(paste("Solution path for Observation",toString(k)))+
                  theme(plot.title = element_text(hjust = 0.5),
                        panel.background = element_rect(fill = "white"),
                        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)))
        }
        else{
          df = beta_path[,ind]
          colnames(df) = coln[ind]
          fig = ggplot()+xlim(0,1)
          for (j in 1:length(ind)){
            xx = data.frame(coef = colnames(df)[j])
            fig = fig + geom_function(data = xx, fun = plot_helper, 
                                      args = list(df[,j]), aes(color = coef)) 
          }
          fig = fig + labs(color = "Feature")
          print(fig + xlab(TeX("Case weight $\\omega$")) + ylab(TeX("Coefficient estimate $\\hat{\\beta}$")) + 
                  theme(plot.title = element_text(hjust = 0.5),
                        panel.background = element_rect(fill = "white"),
                        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)))
        }
      }
    }
    else{
      print(paste('there is no sign change for case',toString(k)))
    }
  }
}
  

#' extract coefficients from a SolPathLooLasso object 
#' @param object      fitted SolPathLooLasso object 
#' @param plot        If plot=1, this function plots the linear approximated solution path.
#'                    If plot=2, plot exact solution path.
#' @return The solution path plot
#' @seealso predict, plot methods
#' @examples
#' library(lars)
#' data(diabetes)
#' attach(diabetes)
#' obj = SolPathLooLasso(x,y,k = 182, 5, 'lambda')
#' 
#' detach(diabetes)
#' 
#' set.seed(100)
#' x = matrix(rnorm(50*200),nrow=50)
#' y = x[,1:5]%*%c(5,4,3,2,1) + rnorm(50)
#' ans = SolPathLooLasso(x,y,k = 1, s = 0.7, mode = "fraction")
#' plot(ans)
#' predict(ans)
#' coef(ans)
#' @export
## S3 method for class 'SolPathLooLasso'
coef.SolPathLooLasso <- function(object){ans = predict(object, type = "coefficient")
    data.frame(coef = c(ans$loobeta0,ans$loobeta),row.names = c('intercept',1:length(ans$loobeta)))}


