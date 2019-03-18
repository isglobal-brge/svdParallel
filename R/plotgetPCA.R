

plot.getPCA <- function(x, type="variables", comps=c(1,2), ...){

  if(!inherits(x, "getPCA"))
    stop("x must be an object of class pcaJacobi")
  
  type.plot <- charmatch(type, c("variables", "individuals"))
  
  if (is.na(type.plot))
    stop("type must be 'variables' or 'individuals'")
  
  if (type.plot==1) {
    plot(var.coord[,comps], xlim=c(-1,1), ylim=c(-1,1), 
         main="Variables", type="n", 
         xlab=paste("PC",comps[1],"(",round(variance[comps[1]]*100,2),"%)"),
         ylab=paste("PC",comps[2],"(",round(variance[comps[2]]*100,2),"%)"), ...)
    abline(h=0, lty=2)
    abline(v=0, lty=2)
    lab.var <- abbreviate(rownames(var.coord))
    text(var.coord[, comps], adj=c(1,-1), labels=lab.var, ...)
    arrows(rep(0,10), rep(0,10), var.coord[,comps[1]], var.coord[,comps[2]],
         length=0.1,angle=20)
    curve(sqrt(1-x^2),-1,1,add=T)
    curve(-sqrt(1-x^2),-1,1,add=T)
  }

  else{
    plot(Y[,1:2], main="Individuals",
         type="n", xlab=paste("PC",comps[1],"(",round(variance[comps[1]]*100,2),"%)"),
         ylab=paste("PC",comps[2],"(",round(variance[comps[2]]*100,2),"%)"), ...)
    abline(h=0)
    abline(v=0)
    points(Y[,1], Y[,2], pch=16, ...)
  }
}