source('../CaseWeightLasso/R/common.R')


#' Plot the residuals and leverages of all observations against their case influence when penalty is fixed 
#' @param X           matrix n by p      design matrix
#' @param y           vector n by 1      response vector
#' @param lambda      number or vector   if only one or a subset of penalty parameters are of interest
#' @param studentize  Boolean            whether to studentize the residual 
#' @return A plot includes all observations. X axis: leverage(when w=1), Y axis: residuals
#' @description Each point in the plot represents a observation and its size indicates its case influence for the lasso.
#'      Observations that are identified as influential points (having case influence greater than the threshold) is marked in red. 
#' @examples
#' library(lars)
#' data(diabetes)
#' attach(diabetes)
#' ResidLevPlot(x,y,1, FALSE)
#' detach(diabetes)
#' 
#' x = matrix(rnorm(50*200),nrow=50)
#' y = x[,1:5]%*%c(5,4,3,2,1) + rnorm(50)
#' ResidLevPlot(x,y,3)
#' @export
ResidLevPlot <- function(X, y, lambda, studentize = TRUE){
  # result = CD_one_fraction(X,y,lambda)
  X = centralize(X)
  n = dim(X)[1]
  p = dim(X)[2]

  result = CookDisLasso(X, y, lambda = lambda, threshold = TRUE, plot = FALSE)
  residy = y - mean(y) - (X%*%result$beta_table)
  X_ = X[,result$beta!=0]
  lev = diag(X_%*%solve(crossprod(X_))%*%t(X_))

  if (n<=p){
    denom = 1
  } else{
    denom = sum(lm(y~X)$residual**2)/(n-p-1)*(p+1)
  }

  if (studentize){
    df_ = data.frame(index = 1:n, lev = lev, res = residy/sqrt(denom*(1-lev)), Case_Influence = drop(result$CD_Mat))
    ylabel = 'studentized residual'
  }else{
    df_ = data.frame(index = 1:n, lev = lev, res = residy, Case_Influence = drop(result$CD_Mat))
    ylabel = 'residual'
  }
  ggplot(data=df_, aes(x=lev, y=res, size=Case_Influence)) +
    geom_point(aes(color = Case_Influence > result$threshold[1,2])) +
    scale_color_manual(values = c("TRUE" = "red", "F" = "black"), guide = FALSE) +
    xlab(TeX('leverage')) +
    ylab(TeX(ylabel)) +
    geom_hline(yintercept=0, linetype="dashed", color = "gray") +
    theme(panel.background = element_rect(fill = "white"),
          legend.position = 'None',
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
  return(list(influential=which(df_$Case_Influence > result$threshold[1,2])))
}

