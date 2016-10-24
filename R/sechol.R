####################################################################################################
###Schnabel/Eskow Cholesky Factorization                                                         ###
###Source: http://artsci.wustl.edu/~jgill/papers/sechol.s                                        ###
###http://dash.harvard.edu/bitstream/handle/1/4214881/King_HessianNotInvertible.pdf?sequence=2   ###
####################################################################################################

#################################################################################
# This is the generalized cholesky routine.  Reference:  	                #
#										#
# Gill, Jeff and Gary King. ``What to do When Your Hessian is Not Invertible:   #
# Alternatives to Model Respecification in Nonlinear Estimation,'' Sociological #
# Methods and Research, Vol. 32, No. 1 (2004): Pp. 54--87. 			#
#										#
# Abstract									#
# 										#
# What should a researcher do when statistical analysis software terminates 	#
# before completion with a message that the Hessian is not invertable? The 	#
# standard textbook advice is to respecify the model, but this is another way 	#
# of saying that the researcher should change the question being asked. 	#
# Obviously, however, computer programs should not be in the business of 	#
# deciding what questions are worthy of study. Although noninvertable 		#
# Hessians are sometimes signals of poorly posed questions, nonsensical 	#
# models, or inappropriate estimators, they also frequently occur when 		#
# information about the quantities of interest exists in the data, through 	#
# the likelihood function. We explain the problem in some detail and lay out 	#
# two preliminary proposals for ways of dealing with noninvertable Hessians 	#
# without changing the question asked. 						#
# 										#
# Also available is the software to implement the procedure described in this 	#
# paper in Gauss format.							#
#################################################################################

#' @export
sechol <- function(A)  {
  n <- nrow(A)
  L <- matrix(rep(0,n*n),ncol=ncol(A))
  macheps <- 2.23e-16; tau <- macheps^(1/3)  # made to match gauss
  gamm <- max(A)
  deltaprev <- 0
  Pprod <- diag(n)
  if (n > 2)  {
    for (k in 1:(n-2))  {
      if( (min(diag(A[(k+1):n,(k+1):n]) - A[k,(k+1):n]^2/A[k,k]) < tau*gamm)
          && (min(svd(A[(k+1):n,(k+1):n])$d)) < 0) {
        dmax <- order(diag(A[k:n,k:n]))[(n-(k-1))]
        if (A[(k+dmax-1),(k+dmax-1)] > A[k,k])  {
          print(paste("iteration:",k,"pivot on:",dmax,"with absolute:",(k+dmax-1)))
          P <- diag(n)
          Ptemp <-  P[k,]; P[k,] <- P[(k+dmax-1),]; P[(k+dmax-1),] = Ptemp
          A <- P%*%A%*%P
          L <- P%*%L%*%P
          Pprod <- P%*%Pprod
        }
        g <- rep(0,length=(n-(k-1)))
        for (i in k:n)  {
          if (i == 1) sum1 <- 0
          else sum1 <- sum(abs(A[i,k:(i-1)]))
          if (i == n) sum2 <- 0
          else sum2 <- sum(abs(A[(i+1):n,i]))
          g[i-(k-1)] <- A[i,i] - sum1 - sum2
        }
        gmax <- order(g)[length(g)]
        if (gmax != k)  {
          print(paste("iteration:",k,"gerschgorin pivot on:",gmax,"with absolute:",(k+gmax-1)))
          P <- diag(ncol(A))
          Ptemp <-  P[k,]; P[k,] <- P[(k+dmax-1),]; P[(k+dmax-1),] = Ptemp
          A <- P%*%A%*%P
          L <- P%*%L%*%P
          Pprod <- P%*%Pprod
        }
        normj <- sum(abs(A[(k+1):n,k]))
        delta <- max(0,deltaprev,-A[k,k]+max(normj,tau*gamm))
        if (delta > 0)  {
          A[k,k] <- A[k,k] + delta
          deltaprev <- delta
        }
      }
      L[k,k] <- A[k,k] <- sqrt(A[k,k])
      for (i in (k+1):n)  {
        L[i,k] <- A[i,k] <- A[i,k]/L[k,k]
        A[i,(k+1):i] <- A[i,(k+1):i] - L[i,k]*L[(k+1):i,k]
        if(A[i,i] < 0) A[i,i] <- 0
      }
    }
  }
  A[(n-1),n] <- A[n,(n-1)]
  eigvals <- eigen(A[(n-1):n,(n-1):n])$values
  delta <- max(0,deltaprev,
               -min(eigvals)+tau*max((1/(1-tau))*(max(eigvals)-min(eigvals)),gamm))
  if (delta > 0)  {
    print(paste("delta:",delta))
    A[(n-1),(n-1)] <- A[(n-1),(n-1)] + delta
    A[n,n] <- A[n,n] + delta
    deltaprev <- delta
  }
  L[(n-1),(n-1)] <- A[(n-1),(n-1)] <- sqrt(A[(n-1),(n-1)])
  L[n,(n-1)] <- A[n,(n-1)] <- A[n,(n-1)]/L[(n-1),(n-1)]
  L[n,n] <- A[n,n] <- sqrt(A[n,n] - L[n,(n-1)]^2)
  return(t(Pprod)%*%t(L)%*%t(Pprod))
}
