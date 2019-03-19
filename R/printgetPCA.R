##' Print PCA
##' 
##' @title Print of the main results of the PCA performed.
##' @param x an object of class getPCA
##' @param ... other parameters
##' @export 

print.getPCA <- function(x, ...){
  cat("**Results for the Principal Component Analysis using method '", x$method, "' **\n", sep="")
  cat("\n")
  mvar = matrix(x$var,nrow=1)
  mpercvar = matrix(x$percvar,nrow=1)
  colnames(mvar) <- paste("PC",1:ncol(mvar),sep="")
  colnames(mpercvar) <- paste("PC",1:ncol(mpercvar),sep="")
  cat("Variance of the components\n")
  print(mvar)
  cat("\n")
  cat("Percentatge of variance explained for each component\n")
  print(mpercvar*100)
}
