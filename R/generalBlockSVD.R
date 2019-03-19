##' SVD using an incremental SVD algorithm.
##'
##' Singular values and left singular vectors of a real nxp matrix 
##' @title SVD using an incremental SVD algorithm.
##' @param x a real nxp matrix
##' @param k number of local SVDs to concatenate at each level 
##' @param q number of levels
##' @export generalBlockSVD
##' @examples
##' (V <- (matrix(1:30, nrow=5, ncol=6)))
##' generalBlockSVD(V,2,3)
##' all.equal(generalBlockSVD(V,2,3)$u, base::svd(V)$u)
##' @return a list of two components with the singular values and left singular vectors of the matrix


generalBlockSVD <- function(A, k=2, q=1){
	n <- nrow(A)
	p <- ncol(A)
	M <- k^q
	
	if(M>p)
	  stop("k^q must not be greater than the number of columns of the matrix")
	Ai <- list()
	p <- ncol(A)
	cols <- seq(1,p,1)
	index <- split(cols, cut(cols,breaks=M))
	Ai <- lapply(index, function(v,index) v[,index], v=A) 
	
	for(j in 1:q){
		svdj <- lapply(Ai, svd)
		Ai <- list()
		for(i in 1:(length(svdj)/k)){
			if(length(svdj[[((i-1)*k+1)]]$d)>1){
				Ai[[i]] <- crossprod(t(svdj[[(i-1)*k+1]]$u),diag(svdj[[(i-1)*k+1]]$d))
				for(l in 1:(k-1)){
				  if(length(svdj[[((i-1)*k+1+l)]]$d)>1){
					  Ai[[i]] <- cbind(Ai[[i]],crossprod(t(svdj[[(i-1)*k+1+l]]$u),diag(svdj[[(i-1)*k+1+l]]$d)))
				  }
				  else{
				    Ai[[i]] <- cbind(Ai[[i]],crossprod(t(svdj[[(i-1)*k+1+l]]$u),svdj[[(i-1)*k+1+l]]$d))
				  }
				}
			}
			else{
				Ai[[i]] <- crossprod(t(svdj[[(i-1)*k+1]]$u),svdj[[(i-1)*k+1]]$d)
				for(l in 1:(k-1)){
				  if(length(svdj[[((i-1)*k+1+l)]]$d)>1){
				    Ai[[i]] <- cbind(Ai[[i]],crossprod(t(svdj[[(i-1)*k+1+l]]$u),diag(svdj[[(i-1)*k+1+l]]$d)))
				  }
				  else{
				    Ai[[i]] <- cbind(Ai[[i]],crossprod(t(svdj[[(i-1)*k+1+l]]$u),svdj[[(i-1)*k+1+l]]$d))
				  }
				}
			}
						 	
		}	
	}
	
	if(length(Ai)>1){
	  Ai <- Reduce(cbind,Ai)
	}
	
	svdA <- svd(Ai[[1]])
	v <- crossprod(A,svdA$u)
	v <- sweep(v,2,svdA$d,FUN="/")
	ans <- list(d=svdA$d, u=svdA$u, v=v)
	return(ans)
}
