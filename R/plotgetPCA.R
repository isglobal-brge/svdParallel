##' Plots PCA
##' 
##' @title Plot of the variables and individuals in PCA.
##' @param x an object of class getPCA
##' @param type variables for the variables plot or individuals for the individuals plot
##' @param comps components to be plot. Default c(1,2)
##' @param ... other parameters of the plot
##' @export 


plot.getPCA <- function(x, type="variables", comps=c(1,2), ...){

  if(!inherits(x, "getPCA"))
    stop("x must be an object of class pcaJacobi")
  
  type.plot <- charmatch(type, c("variables", "individuals"))
  
  if (is.na(type.plot))
    stop("type must be 'variables' or 'individuals'")
  
  if( max(comps)>ncol(x$varcoord) || min(comps)<1)
    stop("incorrect number of components")
  
  if (type.plot==1) {
    plot(x$varcoord[,comps], xlim=c(-1,1), ylim=c(-1,1), 
         main="Variables", type="n", 
         xlab=paste("PC",comps[1],"(",round(x$percvar[comps[1]]*100,2),"%)"),
         ylab=paste("PC",comps[2],"(",round(x$percvar[comps[2]]*100,2),"%)"), ...)
    abline(h=0, lty=2)
    abline(v=0, lty=2)
    lab.var <- abbreviate(rownames(x$varcoord))
    text(x$varcoord[, comps], adj=c(1,-0.5), pos=c(1,2), labels=lab.var, ...)
    arrows(rep(0,10), rep(0,10), x$varcoord[,comps[1]], x$varcoord[,comps[2]],
         length=0.1,angle=20)
    curve(sqrt(1-x^2),-1,1,add=T)
    curve(-sqrt(1-x^2),-1,1,add=T)
  }

  else{
    plot(x$Y[,1:2], main="Individuals",
         type="n", xlab=paste("PC",comps[1],"(",round(x$percvar[comps[1]]*100,2),"%)"),
         ylab=paste("PC",comps[2],"(",round(x$percvar[comps[2]]*100,2),"%)"), ...)
    abline(h=0)
    abline(v=0)
    points(x$Y[,1], x$Y[,2], pch=16, ...)
  }
}
