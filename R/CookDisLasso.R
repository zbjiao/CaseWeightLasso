source('../CaseWeightLasso/R/common.R')


#' Calculate case influence for the lasso and plot Case influence graph
#' @param X           matrix n by p      design matrix
#' @param y           vector n by 1      response vector
#' @param k           number or vector   if one or a subset of observations are of interest
#' @param fineness    float <=0.05       how fine the lambda candidate set is
#' @param lambda      number or vector   if only one or a subset of penalty parameters are of interest
#' @param threshold   Boolean            whether to calculate the threshold for influential points
#' @param plot        Boolean            whether to plot the case influence graph
#' @return CD_Mat          matrix n by f      matrix of cook's distance at each fraction
#'         Lambda_list     vector f by 1      a vector of lambdas at each fraction
#'         fraction        vector f by 1      a vector of fractions based on required fineness
#'         threshold_table matrix f by 2      threshold at each fraction 
#'         beta_hat_table  matrix p by f      a matrix of betahats at each fraction
#'         fig        the case influence graph
#' @description
#' Using case-weight lasso model and the solution path algorithm, this function calculates 
#' the case influence for the lasso for all observations(for specified subset), all fractions(for specified lambdas) and also gives 
#' case influence graph. 
#' @details
#' 1. threshold needs case influence of all observations to calculate. So threshold should be set to False if the functions only calculate 
#' case influence for a subset of observations,i.e. k is not NULL.\\
#' 2. fineness automatically chooses a list of lambda candidates from 0 to maximum lambda that penalizes all beta to 0. If lambda candidate 
#' is pre-specified, i.e. lambda is not NULL, fineness should be set to NULL.
#' @examples
#' library(lars)
#' data(diabetes)
#' attach(diabetes)
#' CookDisLasso(x,y)
#' detach(diabetes)
#' 
#' set.seed(1)
#' x = matrix(rnorm(50*200),nrow=50)
#' y = x[,1:5]%*%c(5,4,3,2,1) + rnorm(50)
#' CookDisLasso(x,y, k=1:10, lambda=c(1,2), threshold = FALSE, plot=FALSE)
#' @export
CookDisLasso <- function(X, y, k = NULL, fineness = NULL, lambda = NULL, threshold = TRUE, plot = TRUE) {
  if (!is.null(k) && threshold) {
    warning("Unable to calculate threshold without all observations.", call. = FALSE)
    threshold = FALSE
  }
  if (!is.null(lambda) && !is.null(fineness)) {
    warning("'lambda' and 'fineness' cannot be specified together.", call. = FALSE)
    fineness = NULL
  }
  # if (!is.null(lambda) && plot) {
  #  warning("plot is not suggested if not using fineness.", call. = FALSE)
  # }
  if (is.null(fineness) && is.null(lambda)) {
    fineness = 100
  }
  
  X = centralize(X)
  n = dim(X)[1]
  p = dim(X)[2]

  # record all gram inverses we come across during the calculation 
  XX_recorder = list()
  XX_recorder[[paste(rep(0,p),collapse = '')]] = matrix(0,nrow=p,ncol=p)
  
  lambda_max = max(abs(t(X)%*%y))
  
  obj = lars(X,y,type='lasso',use.Gram = !(n<p | p>500))
  nbeta1 = drop(abs(obj$beta) %*% rep(1, p))
  ybar = mean(y)
  beta_hat_table = c()
  if (!is.null(fineness)){
    fraction = seq(0.99,0,length.out=fineness)
    for (f in fraction){
      beta_hat = as.vector(predict.lars(obj,mode='fraction', s=f, type = 'coefficients')$coefficients)
      beta_hat_table = cbind(beta_hat_table, beta_hat)
    }
    nbeta2 = drop(rep(1, p) %*% abs(beta_hat_table))
    l_list = approx(nbeta1, c(obj$lambda,0), xout = nbeta2)$y
  } else
  {    
    l_list = sort(lambda[lambda<lambda_max])
    for (l in l_list){
      beta_hat = as.vector(predict.lars(obj,mode='lambda', s=l, type = 'coefficients')$coefficients)
      beta_hat_table = cbind(beta_hat_table, beta_hat)
    }
    fraction = approx(c(obj$lambda,0), nbeta1/nbeta1[length(nbeta1)], xout = l_list)$y
  }
  
  CD_Mat = c()
  
  if (is.null(k)){
    indices_vec = 1:n
  }
  else{
    indices_vec = k
  }
  
  if (n<=p){
    denom = 1
  } else{
    denom = sum(lm(y~X)$residual**2)/(n-p-1)*(p+1)
  }
  
  if (0 %in% l_list){
    starter = 2
    if (n<=p){
      cd0 = rep(NA,length(indices_vec))
    } else{
      cd0 = cooks.distance(lm(y~X))[indices_vec]
    }
  } else{
    starter = 1
  }
  
  
  for (i in starter:length(fraction)){
    lambda = l_list[i]
    
    beta_hat_backup = beta_hat_table[,i]
    s_backup = sign(beta_hat_backup)
    A_backup = s_backup!=0
    # derivative of f in terms of covariates in nonactive set
    d_hat_backup = t(X[,!A_backup])%*%y - t(X[,!A_backup])%*%X[,A_backup]%*%beta_hat_backup[A_backup]
    
    A_id = paste(A_backup*1,collapse = '')
    if (A_id %in% names(XX_recorder)){
      XX_backup = XX_recorder[[A_id]]
    }
    else{
      XX_backup = solve(t(X[,A_backup])%*%X[,A_backup])
      XX_recorder[[A_id]] = XX_backup
    }
    
    if (sum(A_backup)==0){
      h_backup = matrix(0,ncol=n,nrow=n)
    }
    else{
      h_backup = X[,A_backup]%*%XX_backup%*%t(X[,A_backup])
    }
    
    CD_list = c()
    
    # print(c(k,n,K))
    A = s_backup!=0
    if (sum(A)==1){
      y_tilde = X[,A]*beta_hat_backup[A] + ybar
    } 
    else{
      y_tilde = X[,A]%*%beta_hat_backup[A] + ybar
    }
    
    for (k in indices_vec){
      beta_hat = beta_hat_backup
      s = s_backup
      A = A_backup
      
      d_hat = d_hat_backup
      XX = XX_backup
      hk = h_backup[,k]
      
      # current predictor of yk 
      yk_hat = drop(ybar + X[k,A] %*% beta_hat[A])
      
      # beta path records beta hat's value at each breakpoint
      # beta_path = c(beta_hat)
      # and sign change
      # s_path = c(s)
      # and omega value at each breakpoint
      # hkk_path = c()
      w = 1
      while (T){
        # hkk_path = c(hkk_path, hk[k])
        xi = get_xi(w, hk[k], n)
        bias = yk_hat - y[k]
        if (sum(A) == 0){
          xi_cand1 = c()
        }
        else{
          slope_beta = XX%*%X[k,A]*bias
          xi_cand1 = -beta_hat[A]/slope_beta
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
        # slope_d = (X[k,!A] - t(X[,!A])%*%hk)*bias
        # xi candidates
        # xi_cand1 = -beta_hat[A]/slope_beta
        xi_cand2p = (-lambda-d_hat)/slope_d
        xi_cand2m = (lambda-d_hat)/slope_d
        xi_cand = c(xi_cand1, xi_cand2p, xi_cand2m)
        xi_cand0 = min(xi_cand[xi_cand>(xi+0.00001)],Inf)
        # print(c(xi_cand0, xi))
        ind = which(xi_cand == xi_cand0)
        
        # update beta
        if (sum(A) > 0){
          beta_hat[A] = beta_hat[A] + min(get_xi(0,hk[k],n),xi_cand0)*slope_beta
        }
        # beta_path = rbind(beta_path, beta_hat)
        beta_hat0 = ybar + min(get_xi(0,hk[k],n),xi_cand0) * bias/n
        
        # if the xi is off the bound, stop the algorithm
        # print(get_xi(0,hk[k],n))
        # print(beta_hat[A])
        
        if (xi_cand0 > get_xi(0,hk[k],n)){
          # hkk_path = c(hkk_path, hk[k])
          # s_path = rbind(s_path, s)
          break
        }
        # if not, locate the covariate (go to or leave the active set)
        else if(ind<=sum(A)){
          coef = which(cumsum(A)==ind)[1]
          # print(c('-',coef))
          A[coef] = F
          s[coef] = 0
          # check if active set is empty 
          if (sum(A) == 0){
            hk[k] = 0
            # hkk_path = c(hkk_path, 0, 0)
            # s_path = rbind(s_path, s, s)
            # beta_path = rbind(beta_path, beta_hat)
            break
          }
        }
        else if(ind>p){
          coef = which(cumsum(1-A)==(ind-p))[1]
          # print(c('+',coef))
          A[coef] = T
          s[coef] = 1
        }
        else{
          coef = which(cumsum(1-A)==(ind-sum(A)))[1]
          # print(c('+',coef))
          A[coef] = T
          s[coef] = -1
        }
        
        # s_path = rbind(s_path, s)
        # update omega, XX, beta_hat, d_hat
        w = get_w(xi_cand0, hk[k], n)
        A_id = paste(A*1,collapse = '')
        if (A_id %in% names(XX_recorder)){
          XX = XX_recorder[[A_id]]
        }
        else{
          XX = solve(crossprod(X[,A]))
          # XX = solve(t(X[,A])%*%X[,A])
          XX_recorder[[A_id]] = XX
        }
        beta_hat[A] = XX%*%(t(X[,A])%*%y-lambda*s[A])
        beta_hat[!A] = 0
        
        # d_hat_backup = t(X[,!A_backup])%*%y - t(X[,!A_backup])%*%X[,A_backup]%*%beta_hat_backup[A_backup]
        d_hat = t(X[,!A])%*%y - t(X[,!A])%*%X[,A]%*%beta_hat[A]
        hk = X[,A]%*%XX%*%X[k,A]
        yk_hat = drop(ybar + X[k,A] %*% beta_hat[A])
        # print(c(yk_hat,w))
      }
      
      #--------------------
    
      A = s!=0
      A_id = paste(A*1,collapse = '')
      if (sum(A) == 0){
        y_last = ybar + get_xi(0, hk[k], n)*1/n*(ybar - y[k])
      }
      else{
        xxx =  X[,A]%*%XX_recorder[[A_id]]
        y_hat =  xxx%*%(t(X[,A])%*%y-lambda*s[A]) + ybar
        y_last = y_hat + get_xi(0, hk[k], n)*(xxx%*%X[k,A]+1/n)*(y_hat[k] - y[k])
      }
      
      CD_list = c(CD_list, sum((y_last - y_tilde)**2))
    }
    
    CD_Mat = cbind(CD_Mat, CD_list)
  }
  
  CD_Mat = CD_Mat/denom
  if (starter == 2){
    CD_Mat = cbind(cd0,CD_Mat)
  }
  
  colnames(CD_Mat) <- NULL
  colnames(beta_hat_table) <- NULL
  rownames(threshold) <- NULL
  ans = list(CD_Mat = CD_Mat, Lambda_list = l_list, fraction = fraction, beta_table = beta_hat_table)

  if (threshold){
    value = sqrt(apply(CD_Mat,2,var)/2)*qchisq(0.95,1)
    threshold_table = cbind(fraction,value)
    ans$threshold = threshold_table
  }
  
  if (plot){
    plot_table = t(rbind(fraction, CD_Mat))
    colnames(plot_table) = c('fraction',1:length(indices_vec))
    df = as_tibble(plot_table)%>%gather('obs','val',-fraction)
    fig = ggplot(data=df,aes(x=fraction,y=val,group=obs,color=obs)) +
      geom_line() + xlim(c(0,NA))+ labs(x = "|coef|/max|coef|")+
      theme(plot.title = element_text(hjust = 0.5),legend.position = "none") +labs(y="Cook's distance under Lasso")+
      theme(panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(colour = "grey",linewidth = 0.5),
            panel.grid.minor = element_line(colour = "white",linewidth = 0.5),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
    if (threshold){
      colnames(threshold_table)=c('fraction','val')
      df2 = as_tibble(threshold_table)%>%mutate(obs='threshold')
      fig = fig + geom_line(data=df2,aes(x=fraction, y=val),linetype = "dotted", color = 'red1')
    }
    ans$fig = fig
  }
  ans
}

