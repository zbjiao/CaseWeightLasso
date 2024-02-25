source('../CaseWeightLasso/R/common.R')


#' Output the solution path of observation k when lambda is fixed. 
#' Only covariates that have sign change are included in the solution path plot
#' @param X           matrix n by p      design matrix
#' @param y           vector n by 1      response vector
#' @param k           integer            observation of interest
#' @param lambda      float > 0          penalty parameter
#' @param plot        0                  no plot 
#'                    1                  plot approximate 
#'                    2                  plot exact solution path
#' @param lb          float < 1          the lower bound of omega 
#' @return w_path          vector b by 1      a vector of breakpoints
#'         hkk_path        vector b by 1      leverages at each breakpoint
#'         beta_path       matrix b by p      beta estimate at each breakpoint
#'         s_path          matrix b by p      beta hat's sign at each breakpoint
#'         beta0_path      vector b by 1      beta0 at each breakpoint 
#'         l1norm          vector b by 1      l1 norm of beta at each breakpoint
#'         the solution path visualization 
#' @examples
#' library(lars)
#' data(diabetes)
#' attach(diabetes)
#' SolPathLasso(x,y,k = 182, lambda = 5, plot = 2)
#' detach(diabetes)
#' 
#' set.seed(100)
# x = matrix(rnorm(50*200),nrow=50)
# y = x[,1:5]%*%c(5,4,3,2,1) + rnorm(50)
# SolPathLasso(x,y,k = 1, lambda = 1, plot = 2)
#' @export
SolPathLasso <- function(X,y,k=1,lambda=50, plot=0, lb = 0){
  
  X = centralize(X)

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
  
  # obtain lasso solution
  # fit = glmnet(X,y,family="gaussian",lambda=lambda/n,standardize=FALSE,thresh=1e-16,intercept=T)
  # beta_hat = as.vector(fit$beta)
  obj = lars(X,y,type='lasso',use.Gram = !(n<p | p>500))
  beta_hat = coef(obj,s=lambda,mode='lambda')
  # beta_hat = as.vector(predict.lars(obj,mode='lambda', s=lambda, type = 'coefficients')$coefficients)
  
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
    ind = c(which((beta_hat>0) & (beta1<0)),  which((beta_hat<0) & (beta1>0)))
    if(plot){
      coln = colnames(X)
      plot_helper0 <- function(w,ind){
        A = get_xi(w,hk[k],n)%*%t(XX%*%X[k,]*(yk_hat - y[k]))
        as.vector(t(A[,ind]) + (XX%*%(t(X)%*%y))[ind])
      }
      if (length(ind)>0){
        if (length(ind) == 1){
          xx = data.frame(coef = coln[ind])
          print(ggplot()+xlim(lb,1) + xlab(TeX("Case weight $\\omega$")) + 
                  ylab(TeX("Coefficient estimate $\\hat{\\beta}$")) + 
                  geom_function(data = xx, fun = plot_helper0, args = list(ind),aes(color = coef))+ 
                  ggtitle(paste("Solution path for Case",toString(k)))+
                  theme(plot.title = element_text(hjust = 0.5),
                        panel.background = element_rect(fill = "white"),
                        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)))
        }
        else{
          fig = ggplot()+xlim(lb,1)
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
      return(list(w_path = c(1,0),hkk_path = c(hk[k], hk[k]), 
                  beta_path = rbind(beta_hat, beta1), s_path = rbind(s,s = sign(beta1)), 
                  beta0_path = c(ybar, beta10),l1norm = sum(abs(beta_hat))))
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
      beta_hat[A] = beta_hat[A] + min(get_xi(lb,hk[k],n),xi_cand0)*slope_beta
    }
    beta_path = rbind(beta_path, beta_hat)
    beta_hat0 = ybar + min(get_xi(lb,hk[k],n),xi_cand0) * bias/n
    beta0_path = c(beta0_path, beta_hat0)
    
    # if the xi is off the bound, stop the algorithm
    if (xi_cand0 > get_xi(lb,hk[k],n)){
      w_path = c(w_path, lb)
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
        w_path = c(w_path, get_w(xi_cand0, hk[k], n), lb)
        hkk_path = c(hkk_path, 0, 0)
        s_path = rbind(s_path, s, s)
        beta_path = rbind(beta_path, beta_hat)
        beta0_path = c(beta0_path, ybar + get_xi(lb,0,n) * (ybar-y[k])/n)
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
  if(plot){
    plot_helper <-function(x, df){
      i = findInterval(-x, -w_path, rightmost.closed = T)
      beta1 = df[i]
      beta2 = df[i+1]
      w1 = w_path[i]
      w2 = w_path[i+1]
      hkk = hkk_path[i]
      beta1+(beta2-beta1)*(get_xi(x, hkk, n) - get_xi(w1, hkk, n))/(get_xi(w2, hkk, n) - get_xi(w1, hkk, n))
    }
    if (is.null(colnames(X))){
      coln = 1:p
    }
    else{
      coln = colnames(X)
    }
    num_z = apply(beta_path, 2, function(c) sum(abs(c)< 1e-10))
    print(num_z)
    ind = which(num_z>0 & num_z<length(w_path))
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
          print(ggplot()+xlim(lb,1) + xlab(TeX("Case weight $\\omega$")) + 
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
          fig = ggplot()+xlim(lb,1)
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
  rownames(beta_path) <- NULL
  rownames(s_path) <- NULL
  
  return(list(w_path = w_path, hkk_path = hkk_path, beta_path = beta_path, 
              s_path = s_path, beta0_path = beta0_path,
              l1norm = as.vector(apply(beta_path, 1, function(x) sum(abs(x))))))
}