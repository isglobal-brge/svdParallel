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


generalBlockSVD <- function(A, k, q){
	n <- nrow(A)
	p <- ncol(A)
	M <- k^q
	Ai <- list()
	for(i in 1:M){
		Ai[[i]] <- A[,(p/M*(i-1)+1):(p/M*i)]
	}
	
	for(j in 1:q){
		svdj <- lapply(Ai, svd)
		Ai <- list()
		for(i in 1:(length(svdj)/k)){
			if(length(svdj[[((i-1)*k+1)]]$d)>1){
				Ai[[i]] <- crossprod(t(svdj[[(i-1)*k+1]]$u),diag(svdj[[(i-1)*k+1]]$d))
				for(l in 1:(k-1)){
					Ai[[i]] <- cbind(Ai[[i]],crossprod(t(svdj[[(i-1)*k+1+l]]$u),diag(svdj[[(i-1)*k+1+l]]$d)))
				}
			}
			else{
				Ai[[i]] <- crossprod(t(svdj[[(i-1)*k+1]]$u),svdj[[(i-1)*k+1]]$d)
				for(l in 1:(k-1)){
					Ai[[i]] <- cbind(Ai[[i]],crossprod(t(svdj[[(i-1)*k+1+l]]$u),svdj[[(i-1)*k+1+l]]$d))
				}
			}
						 	
		}	
	}
	svdA <- svd(Ai[[1]])
	ans <- list(d=svdA$d, u=svdA$u)
	return(ans)
}
