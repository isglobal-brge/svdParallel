##' Print PCA
##' 
##' @title Print of the main results of the PCA performed.
##' @param x an object of class getPCA
##' @param ... other parameters
##' @export 

print.getPCA <- function(x, ...){
  cat("**Results for the Principal Component Analysis using method '", x$method, "' **\n")
  cat("\n")
  cat("Variance of the components\n")
  print(x$var)
  cat("\n")
  cat("Percentatge of variance explained for each component\n")
  print(x$percvar)
}
