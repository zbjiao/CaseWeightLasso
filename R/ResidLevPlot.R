source('../CaseWeightLasso/R/common.R')


#' Plot the residuals and leverages of all observations against their case influence when penalty is fixed 
#' @param X           matrix n by p      design matrix
#' @param y           vector n by 1      response vector
#' @param s           a value indexing the penalty level to consider.
#'                    Its values depends on the mode= argument. By default (mode="fraction"), 
#'                    s should take on values between 0 and 1. 
#' @param mode        Mode="fraction", then s should be a number between 0 and 1, and it refers 
#'                    to the ratio of the L1 norm of the coefficient vector, relative to the 
#'                    norm at the full LS solution if n>p or the minimum L1 norm linear 
#'                    interpolation solution if n<=p. Mode="norm" means s refers to the L1 norm 
#'                    of the coefficient vector. Mode="lambda" uses the lasso regularization 
#'                    parameter for s. Abbreviations allowed.
#' @param studentize  If true, ResidLevPlot studentizes the residual. 
#' @return A plot includes all observations. X axis: leverage(when w=1), Y axis: residuals
#' @description Each point in the plot represents a observation and its size indicates its case influence for the lasso.
#'      Observations that are identified as influential points (having case influence greater than the threshold) is marked in red. 
#' @examples
#' library(lars)
#' data(diabetes)
#' attach(diabetes)
#' ResidLevPlot(x,y,1,'lambda', FALSE)
#' detach(diabetes)
#' 
#' x = matrix(rnorm(50*200),nrow=50)
#' y = x[,1:5]%*%c(5,4,3,2,1) + rnorm(50)
#' ResidLevPlot(x,y,0.7)
#' @export
ResidLevPlot <- function(X, y, s, mode = c("fraction", "norm", "lambda"), studentize = TRUE){
  mode <- match.arg(mode)
  X = centralize(X)
  n = dim(X)[1]
  p = dim(X)[2]

  result = CookDisLasso(X, y, s = s, mode=mode, threshold = TRUE)
  residy = y - mean(y) - (X%*%result$beta_table)
  X_ = X[,result$beta!=0]
  lev = result$lev_history

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
  fig = ggplot(data=df_, aes(x=lev, y=res, size=Case_Influence)) +
    geom_point(aes(color = Case_Influence > result$threshold[1,2])) +
    scale_color_manual(values = c("TRUE" = "red", "F" = "black")) +
    xlab(TeX('leverage')) +
    ylab(TeX(ylabel)) +
    geom_hline(yintercept=0, linetype="dashed", color = "gray") +
    theme(panel.background = element_rect(fill = "white"),
          legend.position = 'None',
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
  return(list(influential=which(df_$Case_Influence > result$threshold[1,2]),
              fig = fig))
}

